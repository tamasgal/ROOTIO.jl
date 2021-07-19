struct ROOTDirectory
    name::AbstractString
    header::ROOTDirectoryHeader
    keys::Vector{TKey}
end

struct ROOTFile
    filename::AbstractString
    format_version::Int32
    header::FileHeader
    fobj::IOStream
    tkey::TKey
    streamers::Streamers
    directory::ROOTDirectory
    customstructs::Dict{String, Type}
    lk::ReentrantLock
end
lock(f::ROOTFile) = lock(f.lk)
unlock(f::ROOTFile) = unlock(f.lk)
function Base.hash(rf::ROOTFile, h::UInt)
    hash(rf.fobj, h)
end

"""
    ROOTFile(filename::AbstractString; customstructs = Dict("TLorentzVector" => LorentzVector{Float64}))

`ROOTFile`'s constructor from a file. The `customstructs` dictionary can be used to pass user-defined
struct as value and its corresponding `fClassName` (in Branch) as key such that `UnROOT` will know
to intepret them, see [`interped_data`](@ref).

See also: [`LazyTree`](@ref), [`LazyBranch`](@ref)

# Example
```julia
julia> f = ROOTFile("test/samples/NanoAODv5_sample.root")
ROOTFile with 2 entries and 21 streamers.
test/samples/NanoAODv5_sample.root
└─ Events
   ├─ "run"
   ├─ "luminosityBlock"
   ├─ "event"
   ├─ "HTXS_Higgs_pt"
   ├─ "HTXS_Higgs_y"
   └─ "⋮"
```
"""
function ROOTFile(filename::AbstractString; customstructs = Dict("TLorentzVector" => LorentzVector{Float64}))
    fobj = Base.open(filename)
    preamble = unpack(fobj, FilePreamble)
    String(preamble.identifier) == "root" || error("Not a ROOT file!")
    format_version = preamble.fVersion

    if format_version < 1000000
        @debug "32bit ROOT file"
        header = unpack(fobj, FileHeader32)
    else
        @debug "64bit ROOT file"
        header = unpack(fobj, FileHeader64)
    end

    # Streamers
    if header.fSeekInfo != 0
        @debug "Reading streamer info."
        seek(fobj, header.fSeekInfo)
        streamers = Streamers(fobj)
    else
        @debug "No streamer info present, skipping."
    end

    seek(fobj, header.fBEGIN)
    tkey = unpack(fobj, TKey)

    # Reading the header key for the top ROOT directory
    seek(fobj, header.fBEGIN + header.fNbytesName)
    dir_header = unpack(fobj, ROOTDirectoryHeader)

    seek(fobj, dir_header.fSeekKeys)
    header_key = unpack(fobj, TKey)

    n_keys = readtype(fobj, Int32)
    keys = [unpack(fobj, TKey) for _ in 1:n_keys]

    directory = ROOTDirectory(tkey.fName, dir_header, keys)

    ROOTFile(filename, format_version, header, fobj, tkey, streamers, directory, customstructs, ReentrantLock())
end

function Base.show(io::IO, f::ROOTFile)
    n_entries = length(f.directory.keys)
    entries_suffix = n_entries == 1 ? "entry" : "entries"
    n_streamers = length(f.streamers)
    streamers_suffix = n_streamers == 1 ? "streamer" : "streamers"
    print(io, typeof(f))
    print(io, " with $n_entries $entries_suffix ")
    println(io, "and $n_streamers $streamers_suffix.")
    print_tree(io, f)
end


UUID(f::ROOTFile) = f.header.fUUID


function streamerfor(f::ROOTFile, name::AbstractString)
    for e in f.streamers.elements
        if e.streamer.fName == name
            return e
        end
    end
    error("No streamer found for $name.")
end


function Base.getindex(f::ROOTFile, s::AbstractString)
    S = _getindex(f, s)
    if S isa Union{TBranch, TBranchElement}
        try # if we can't construct LazyBranch, just give up (maybe due to custom class)
            return LazyBranch(f, S)
        catch
            @warn "Can't automatically create LazyBranch for branch $s. Returning a branch object"
        end
    end
    S
end

@memoize LRU(maxsize = 2000) function _getindex(f::ROOTFile, s)
# function _getindex(f::ROOTFile, s)
    if '/' ∈ s
        @debug "Splitting path '$s' and getting items recursively"
        paths = split(s, '/')
        return f[first(paths)][join(paths[2:end], "/")]
    end
    tkey = f.directory.keys[findfirst(isequal(s), keys(f))]
    @debug "Retrieving $s ('$(tkey.fClassName)')"
    streamer = getfield(@__MODULE__, Symbol(tkey.fClassName))
    lock(f)
    try
        S = streamer(f.fobj, tkey, f.streamers.refs)
        return S
    catch
    finally
        unlock(f)
    end
end

function Base.keys(f::ROOTFile)
    keys(f.directory)
end

function Base.keys(d::ROOTDirectory)
    [key.fName for key in d.keys]
end

Base.keys(t::TTree) = [b.fName for b in t.fBranches.elements]

function Base.getindex(t::T, s::AbstractString) where {T<:Union{TTree, TBranchElement}}
    if '/' ∈ s
        @debug "Splitting path '$s' and getting branches recursively"
        paths = split(s, '/')
        return t[first(paths)][join(paths[2:end], "/")]
    end
    @debug "Searching for branch '$s' in $(length(t.fBranches.elements)) branches."
    for branch in t.fBranches.elements
        if branch.fName == s
            return branch
        end
    end
    missing
end
function Base.getindex(t::TTree, s::Vector{T}) where {T<:AbstractString}
    [t[n] for n in s]
end

reinterpret(vt::Type{Vector{T}}, data::AbstractVector{UInt8}) where T <: Union{AbstractFloat, Integer} = reinterpret(T, data)

"""
    interped_data(rawdata, rawoffsets, ::Type{T}, ::Type{J}) where {T, J<:JaggType}

The function thats interpret raw bytes (from a basket) into corresponding Julia data, based
on type `T` and jagg type `J`.

In order to retrieve data from custom branches, user should defined more speialized
method of this function with specific `T` and `J`. See `TLorentzVector` example.
"""
function interped_data(rawdata, rawoffsets, ::Type{T}, ::Type{J}) where {T, J<:JaggType}
    # there are two possibility, one is the leaf is just normal leaf but the title has "[...]" in it
    # magic offsets, seems to be common for a lot of types, see auto.py in uproot3
    # only needs when the jaggedness comes from TLeafElements, not needed when
    # the jaggedness comes from having "[]" in TLeaf's title
    # the other is where we need to auto detector T bsaed on class name
    # we want the fundamental type as `reinterpret` will create vector
    if J !== Nojagg
        jagg_offset = J===Offsetjagg ? 10 : 0

        # for each "event", the index range is `offsets[i] + jagg_offset + 1` to `offsets[i+1]`
        # this is why we need to append `rawoffsets` in the `readbranchraw()` call
        # when you use this range to index `rawdata`, you will get raw bytes belong to each event
        # Say your real data is Int32 and you see 8 bytes after indexing, then this event has [num1, num2] as real data
        @views [
                ntoh.(reinterpret(
                                  T, rawdata[ (rawoffsets[i]+jagg_offset+1):rawoffsets[i+1] ]
                                 )) for i in 1:(length(rawoffsets) - 1)
               ]
    else # the branch is not jagged
        return ntoh.(reinterpret(T, rawdata))
    end

end

function _normalize_ftype(fType)
    # Taken from uproot4; thanks Jim ;)
    if Const.kOffsetL < fType < Const.kOffsetP
        fType - kOffsetP
    else
        fType
    end
end

const _leaftypeconstlookup = Dict(
                             Const.kBool   => Bool  ,
                             Const.kChar   => Int8  ,
                             Const.kUChar  => UInt8 ,
                             Const.kShort  => Int16 ,
                             Const.kUShort => UInt16,
                             Const.kInt    => Int32 ,
                             Const.kBits   => UInt32,
                             Const.kUInt   => UInt32,
                             Const.kCounter=>UInt32 ,
                             Const.kLong => Int64   ,
                             Const.kLong64 =>  Int64,
                             Const.kULong =>  UInt64,
                             Const.kULong64 =>UInt64,
                             Const.kDouble32 => Float32,
                             Const.kDouble =>   Float64,
                             Const.kFloat => Float32,
                            )

"""
    auto_T_JaggT(branch; customstructs::Dict{String, Type})

Given a branch, automatically return (eltype, Jaggtype). This function is aware of custom structs that
are carried with the parent `ROOTFile`.

This is also where you may want to "redirect" classname -> Julia struct name,
for example `"TLorentzVector" => LorentzVector` here and you can focus on `LorentzVectors.LorentzVector`
methods from here on.

See also: [`ROOTFile`](@ref), [`interped_data`](@ref)
"""
function auto_T_JaggT(branch; customstructs::Dict{String, Type})
# TODO Why is this broken on 1.8?
# @memoize LRU(;maxsize=10^3) function auto_T_JaggT(branch; customstructs::Dict{String, Type})
    leaf = first(branch.fLeaves.elements)
    _type = Nothing
    _jaggtype = JaggType(leaf)
    if hasproperty(branch, :fClassName)
        classname = branch.fClassName # the C++ class name, such as "vector<int>"
        try
            # this will call a customize routine if defined by user
            # see custom.jl
            _custom = customstructs[classname]
            return _custom, _jaggtype
        catch
        end
        m = match(r"vector<(.*)>", classname)
        if m!==nothing
            elname = m[1]
            try
                _custom = customstructs[elname]
                return Vector{_custom}, _jaggtype
            catch
            end
            elname = endswith(elname, "_t") ? lowercase(chop(elname; tail=2)) : elname  # Double_t -> double
            try
                _type = elname == "bool" ?          Bool : _type #Cbool doesn't exist
                _type = elname == "unsigned int" ?  UInt32 : _type #Cunsigned doesn't exist
                _type = elname == "unsigned char" ? Char   : _type
                _type = getfield(Base, Symbol(:C, elname))

                # we know it's a vector because we saw vector<>
                _type = Vector{_type}
            catch
                error("Cannot convert element of $elname to a native Julia type")
            end
        # Try to interpret by leaf type
        else
            leaftype = _normalize_ftype(leaf.fType)
            _type = get(_leaftypeconstlookup, leaftype, nothing)
            isnothing(_type) && error("Cannot interpret type.")
            if branch.fType == Const.kSubbranchSTLCollection
                _type = Vector{_type}
            end
        end
    else
        _type = primitivetype(leaf)
        _type = _jaggtype === Nojagg ? _type : Vector{_type}
    end

    return _type, _jaggtype
end


# read all bytes of DATA and OFFSET from a branch
function readbranchraw(f::ROOTFile, branch)
    nbytes = branch.fBasketBytes
    datas = sizehint!(Vector{UInt8}(), sum(nbytes)) # maximum length if all data are UInt8
    offsets = sizehint!(zeros(Int32, 1), branch.fEntries+1) # this is always Int32
    position = 0
    foreach(branch.fBasketSeek) do seek
        seek==0 && return
        data, offset = readbasketseek(f, branch, seek)
        append!(datas, data)
        # FIXME: assuming offset has always 0 or at least 2 elements ;)
        append!(offsets, (@view offset[2:end]) .+ position)
        if length(offset) > 0
            position = offset[end]
        end
    end
    datas, offsets
end

# Thanks Jim and Philippe
# https://groups.google.com/forum/#!topic/polyglot-root-io/yeC0mAizQcA
# The offsets start at fKeylen - fLast + 4. A singe basket of data looks like this:
#                                           4 bytes          4 bytes
# ┌─────────┬────────────────────────────────┬───┬────────────┬───┐
# │ TKey    │ content                        │ X │ offsets    │ x │
# └─────────┴────────────────────────────────┴───┴────────────┴───┘
#           │←        fLast - fKeylen       →│(l1)                │
#           │                                                     │
#           │←                       fObjlen                     →│
# 3GB cache for baskets
"""
    readbasket(f::ROOTFile, branch, ith)
    readbasketseek(f::ROOTFile, branch::Union{TBranch, TBranchElement}, seek_pos::Int)

The fundamental building block of reading read data from a .root file. Read read one
basket's raw bytes and offsets at a time. These raw bytes and offsets then (potentially) get
processed by [`interped_data`](@ref).

See also: [`auto_T_JaggT`](@ref), [`basketarray`](@ref)
"""
readbasket(f::ROOTFile, branch, ith) = readbasketseek(f, branch, branch.fBasketSeek[ith])

@memoize LRU(; maxsize=3 * 1024^3, by=x -> sum(sizeof, x)) function readbasketseek(
# function readbasketseek(
f::ROOTFile, branch::Union{TBranch, TBranchElement}, seek_pos::Int
)::Tuple{Vector{UInt8},Vector{Int32},Int32}  # just being extra careful
    lock(f)
    seek(f.fobj, seek_pos)
    basketkey = unpack(f.fobj, TBasketKey)
    compressedbytes = compressed_datastream(f.fobj, basketkey)
    unlock(f)

    basketrawbytes = decompress_datastreambytes(compressedbytes, basketkey)

    @debug begin
        ibasket = findfirst(==(seek_pos), branch.fBasketSeek)
        mbcompressed = length(compressedbytes)/1024^2
        mbuncompressed = length(basketrawbytes)/1024^2
        "Read branch $(branch.fName), basket $(ibasket), $(mbcompressed) MB compressed, $(mbuncompressed) MB uncompressed"
    end

    Keylen = basketkey.fKeylen
    contentsize = Int32(basketkey.fLast - Keylen)

    offsetbytesize = basketkey.fObjlen - contentsize - 8

    data = @view basketrawbytes[1:contentsize]
    if offsetbytesize > 0

        #indexing is inclusive on both ends
        offbytes = @view basketrawbytes[(contentsize + 4 + 1):(end - 4)]

        # offsets starts at -fKeylen, same as the `local_offset` we pass in in the loop
        offset = ntoh.(reinterpret(Int32, offbytes)) .- Keylen
        push!(offset, contentsize)
        data, offset, contentsize
    else
        data, Int32[], contentsize
    end
end
