<img style="height:5em;" alt="UnROOT.jl" src="docs/src/assets/unroot.svg"/>

<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-5-orange.svg?style=flat-square)](#contributors-)
<!-- ALL-CONTRIBUTORS-BADGE:END -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tamasgal.github.io/UnROOT.jl/dev)
[![Build Status](https://github.com/tamasgal/UnROOT.jl/workflows/CI/badge.svg)](https://github.com/tamasgal/UnROOT.jl/actions)
[![Codecov](https://codecov.io/gh/tamasgal/UnROOT.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/tamasgal/UnROOT.jl)

UnROOT.jl is a (WIP) reader for the [CERN ROOT](https://root.cern) file format
written entirely in pure Julia, without any dependence on ROOT or Python.

While the ROOT documentation does not contain a detailed description of the
binary structure, the format can be triangulated by other packages like

- [uproot3](https://github.com/scikit-hep/uproot) (Python), see also [UpROOT.jl](https://github.com/JuliaHEP/UpROOT.jl/)
- [groot](https://godoc.org/go-hep.org/x/hep/groot#hdr-File_layout) (Go)
- [root-io](https://github.com/cbourjau/alice-rs/tree/master/root-io) (Rust)
- [Laurelin](https://github.com/spark-root/laurelin) (Java)

Here's a detailed [from-scratch walk through](https://jiling.web.cern.ch/jiling/dump/ROOT_Fileformat.pdf) 
on reading a jagged branch from a ROOT file, recommended for first time contributors or those who just want to learn
about ROOT file format.

Three's also a [discussion](https://github.com/scikit-hep/uproot/issues/401) reagarding the ROOT binary format
documentation on uproot's issue page.

## Status
We support reading all scalar and jagged branches of "basic" types, provide
indexing and iteration interface with a "per branch" basket-cache. There is a low level
API to provide interpretation functionalities for custom types and classes.
As a metric, UnROOT can read all branches (~1800) of CMS NanoAOD including jagged `TLorentzVector` branch.

## Quick Start
The easiest way to access data is through `LazyTree`, which is `<: AbstractDataFrame` and
a thin-wrap around `TypedTable` under the hood. It supports most accessing patterns from
the loved `DataFrames` eco-system.
```julia
julia> using UnROOT

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

julia> mytree = LazyTree(f, "Events", ["Electron_dxy", "nMuon", r"Muon_(pt|eta)$"])
 Row │ Electron_dxy     nMuon   Muon_eta         Muon_pt
     │ Vector{Float32}  UInt32  Vector{Float32}  Vector{Float32}
─────┼───────────────────────────────────────────────────────────
 1   │ [0.000371]       0       []               []
 2   │ [-0.00982]       2       [0.53, 0.229]    [19.9, 15.3]
 3   │ []               0       []               []
 4   │ [-0.00157]       0       []               []
 ⋮   │     ⋮            ⋮             ⋮                ⋮
 
 
julia> mytree[1:3, :nMuon]
3-element Vector{UInt32}:
 0x00000000
 0x00000002
 0x00000000
```

You can iterate through a `LazyTree`:
```julia
julia> for event in mytree
           @show event.Electron_dxy
           break
       end
event.Electron_dxy = Float32[0.00037050247]
```

Only one basket per branch will be cached so you don't have to worry about running out of RAM.
At the same time, `event` inside the for-loop is not materialized until a field is accessed. If your event
is fairly small or you need all of them anyway, you can `collect(event)` first inside the loop.

## Branch of custom struct

We provide an experimental interface for hooking up UnROOT with your custom types
that only takes 2 steps, as explained [here](https://github.com/tamasgal/UnROOT.jl/wiki/CustomBranch).
As a show case for this functionality, the `TLorentzVector` support in UnROOT is implemented
with the said plug-in system.

Alternatively, reading raw data is also possible
using the `UnROOT.array(f::ROOTFile, path; raw=true)` method. The output can
be then reinterpreted using a custom type with the method
`UnROOT.splitup(data, offsets, T::Type; skipbytes=0, jagged=true)`. This provides more fine grain control in case
your branch is highly irregular. You can then define suitable Julia `type` and `readtype` method for parsing these data.
Alternatively, you can of course parse the `data` and `offsets` entirely manually.
Here is it in action, with the help of the `type`s from `custom.jl`, and some data from the KM3NeT experiment:
``` julia
julia> using UnROOT

julia> f = ROOTFile("test/samples/km3net_online.root")
ROOTFile("test/samples/km3net_online.root") with 10 entries and 41 streamers.

julia> data, offsets = array(f, "KM3NET_EVENT/KM3NET_EVENT/snapshotHits"; raw=true)
2058-element Array{UInt8,1}:
 0x00
 0x03
   ⋮
   
julia> UnROOT.splitup(data, offsets, UnROOT.KM3NETDAQHit)
4-element Vector{Vector{UnROOT.KM3NETDAQHit}}:
 [UnROOT.KM3NETDAQHit(1073742790, 0x00, 9, 0x60)......
```

## Main challenges

- ROOT data is generally stored as big endian and is a
  self-descriptive format, i.e. so-called streamers are stored in the files
  which describe the actual structure of the data in the corresponding branches.
  These streamers are read during runtime and need to be used to generate
  Julia structs and `unpack` methods on the fly.
- Performance is very important for a low level I/O library.


## Low hanging fruits

Pick one ;)

- [x] Parsing the file header
- [x] Read the `TKey`s of the top level dictionary
- [x] Reading the available trees
- [x] Reading the available streamers
- [x] Reading a simple dataset with primitive streamers
- [x] Reading of raw basket bytes for debugging
- [ ] Automatically generate streamer logic
- [x] Prettier `show` for `Lazy*`s
- [ ] Clean up `Cursor` use
- [x] Reading `TNtuple` #27
- [x] Reading histograms (`TH1D`, `TH1F`, `TH2D`, `TH2F`, etc.) #48
- [ ] Clean up the `readtype`, `unpack`, `stream!` and `readobjany` construct
- [ ] Refactor the code and add more docs
- [ ] Class name detection of sub-branches
- [ ] High-level histogram interface

## Acknowledgements

Special thanks to Jim Pivarski ([@jpivarski](https://github.com/jpivarski))
from the [Scikit-HEP](https://github.com/scikit-hep) project, who is the
main author of [uproot](https://github.com/scikit-hep/uproot), a native
Python library to read and write ROOT files, which was and is a great source
of inspiration and information for reverse engineering the ROOT binary
structures.

## Behind the scene
<details><summary>Some additional debug output: </summary>
<p>


``` julia
julia> using UnROOT

julia> f = ROOTFile("test/samples/tree_with_histos.root")
Compressed stream at 1509
ROOTFile("test/samples/tree_with_histos.root") with 1 entry and 4 streamers.

julia> keys(f)
1-element Array{String,1}:
 "t1"

julia> keys(f["t1"])
Compressed datastream of 1317 bytes at 1509 (TKey 't1' (TTree))
2-element Array{String,1}:
 "mynum"
 "myval"

julia> f["t1"]["mynum"]
Compressed datastream of 1317 bytes at 6180 (TKey 't1' (TTree))
UnROOT.TBranch
  cursor: UnROOT.Cursor
  fName: String "mynum"
  fTitle: String "mynum/I"
  fFillColor: Int16 0
  fFillStyle: Int16 1001
  fCompress: Int32 101
  fBasketSize: Int32 32000
  fEntryOffsetLen: Int32 0
  fWriteBasket: Int32 1
  fEntryNumber: Int64 25
  fIOFeatures: UnROOT.ROOT_3a3a_TIOFeatures
  fOffset: Int32 0
  fMaxBaskets: UInt32 0x0000000a
  fSplitLevel: Int32 0
  fEntries: Int64 25
  fFirstEntry: Int64 0
  fTotBytes: Int64 170
  fZipBytes: Int64 116
  fBranches: UnROOT.TObjArray
  fLeaves: UnROOT.TObjArray
  fBaskets: UnROOT.TObjArray
  fBasketBytes: Array{Int32}((10,)) Int32[116, 0, 0, 0, 0, 0, 0, 0, 0, 0]
  fBasketEntry: Array{Int64}((10,)) [0, 25, 0, 0, 0, 0, 0, 0, 0, 0]
  fBasketSeek: Array{Int64}((10,)) [238, 0, 0, 0, 0, 0, 0, 0, 0, 0]
  fFileName: String ""


julia> seek(f.fobj, 238)
IOStream(<file test/samples/tree_with_histos.root>)

julia> basketkey = UnROOT.unpack(f.fobj, UnROOT.TKey)
UnROOT.TKey64(116, 1004, 100, 0x6526eafb, 70, 0, 238, 100, "TBasket", "mynum", "t1")

julia> s = UnROOT.datastream(f.fobj, basketkey)
Compressed datastream of 100 bytes at 289 (TKey 'mynum' (TBasket))
IOBuffer(data=UInt8[...], readable=true, writable=false, seekable=true, append=false, size=100, maxsize=Inf, ptr=1, mark=-1)

julia> [UnROOT.readtype(s, Int32) for _ in 1:f["t1"]["mynum"].fEntries]
Compressed datastream of 1317 bytes at 6180 (TKey 't1' (TTree))
25-element Array{Int32,1}:
  0
  1
  2
  3
  4
  5
  6
  7
  8
  9
 10
 10
 10
 10
 10
```
</p>
</details>

## Contributors ✨

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tr>
    <td align="center"><a href="http://www.tamasgal.com"><img src="https://avatars.githubusercontent.com/u/1730350?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Tamas Gal</b></sub></a><br /><a href="https://github.com/tamasgal/UnROOT.jl/commits?author=tamasgal" title="Code">💻</a> <a href="https://github.com/tamasgal/UnROOT.jl/commits?author=tamasgal" title="Documentation">📖</a> <a href="#infra-tamasgal" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="#data-tamasgal" title="Data">🔣</a> <a href="https://github.com/tamasgal/UnROOT.jl/commits?author=tamasgal" title="Tests">⚠️</a></td>
    <td align="center"><a href="https://github.com/Moelf"><img src="https://avatars.githubusercontent.com/u/5306213?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Jerry Ling</b></sub></a><br /><a href="https://github.com/tamasgal/UnROOT.jl/commits?author=Moelf" title="Code">💻</a> <a href="https://github.com/tamasgal/UnROOT.jl/commits?author=Moelf" title="Tests">⚠️</a> <a href="#data-Moelf" title="Data">🔣</a> <a href="https://github.com/tamasgal/UnROOT.jl/commits?author=Moelf" title="Documentation">📖</a></td>
    <td align="center"><a href="https://github.com/8me"><img src="https://avatars.githubusercontent.com/u/17862090?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Johannes Schumann</b></sub></a><br /><a href="https://github.com/tamasgal/UnROOT.jl/commits?author=8me" title="Code">💻</a> <a href="https://github.com/tamasgal/UnROOT.jl/commits?author=8me" title="Tests">⚠️</a> <a href="#data-8me" title="Data">🔣</a></td>
    <td align="center"><a href="https://github.com/aminnj"><img src="https://avatars.githubusercontent.com/u/5760027?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Nick Amin</b></sub></a><br /><a href="https://github.com/tamasgal/UnROOT.jl/commits?author=aminnj" title="Code">💻</a> <a href="https://github.com/tamasgal/UnROOT.jl/commits?author=aminnj" title="Tests">⚠️</a> <a href="#data-aminnj" title="Data">🔣</a></td>
    <td align="center"><a href="https://giordano.github.io"><img src="https://avatars.githubusercontent.com/u/765740?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Mosè Giordano</b></sub></a><br /><a href="#infra-giordano" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a></td>
  </tr>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!
