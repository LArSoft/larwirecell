# wcls wires-generation & consistency scripts

Tooling to generate WCT "wires files" from a LArSoft geometry (so they are
consistent-by-construction with the channel map) and to summarize/validate them.

| script | what it does |
|--------|--------------|
| `wct-wires-file.sh generate` | from a JSON spec → generate FHiCL → run `lar` (WireSchemaFileWriter) → `<era><tech>-<variant>.json.bz2` + a `.sidecar.json` provenance file (input spec + FHiCL; **no** summary) |
| `wct-wires-file.sh summary`  | summarize an existing wires file (wraps `wires-summary.py`); the per-anode and per-plane arrays are excluded unless `--with-anodes` / `--with-planes` |
| `wires-summary.py`  | summarize a wires file: uncompressed-JSON sha256 + totals by default; `--with-planes` adds the per-plane array (wire/channel counts, endpoint bounding boxes, mean/stdev of transverse wire angle, pitch distance, and `dir·X`, mean pitch vector, wire-zero and central-wire centers); `--with-anodes` adds the per-anode array |
| `wires-compare.py`  | compare two wires files (or summaries): exact hash match else structural diff within tolerance; exit 0 = agree |

## Generate a wires file

```bash
cat > fdvd-v6.json <<'EOF'
{ "era": "fd", "technology": "vd", "variant": "3view30deg-v6-wrapped",
  "gdml": "dunevd10kt_3view_30deg_v6_refactored.gdml",
  "gdml_sha1sum": "8edf6b4efeaf8cabcb75e61954e34a4f7f5eea2d",
  "detector_name": "dunevd10kt_3view_30deg_v6_refactored",
  "channel_map": "CRPWireReadoutGeom", "sorter": "GeoObjectSorterCRU60D",
  "pdmap_tool": "PDVDPDMapAlg", "pdmap_file": "PDVD_PDS_Mapping_v09162025.json",
  "grouping": "pair-tpc", "correction": 4 }
EOF
export GDML_PATH=.../dunecore/dunecore/Geometry/gdml
export PDMAP_PATH=.../dunecore/dunecore/ChannelMap
./wct-wires-file.sh generate -i fdvd-v6.json -o out/
```

The side-car `out/<canonical>.sidecar.json` holds the input spec, the generated
FHiCL text, and the channel-map/sorter — but **not** a wire summary (produce that
separately with `summary`).

`grouping`: `tpcset` (default, mapper-faithful one anode per CRU), `per-tpc`
(one single-face anode per CRM), or `pair-tpc` (2 CRM = 1 CRU = one 2-face anode,
the layout of the shipped wrapped v6).

### Runtime requirement — dunecore

The DUNE channel maps (`DUNEWireReadout`, `WireReadoutCRUGeom`,
`CRPWireReadoutGeom`, `ColdBoxWireReadoutGeom`) and sorters (`GeoObjectSorter*`)
live in **dunecore**.  `lar` needs dunecore's plugins on
`LD_LIBRARY_PATH`/`CET_PLUGIN_PATH` and its fcl on `FHICL_FILE_PATH`.  Ensure the
**installed** dunecore is current (an old `DUNEWireReadout` lacking the
`1x8x`/`1x6x` split maps every `dunevd10kt` to the non-wrapping
`WireReadoutCRUGeom`).

## Summarize a wires file

```bash
# default: sha256 + totals only
./wct-wires-file.sh summary out/<canonical>.json.bz2 -o out/<canonical>.summary.json
# add the per-plane and/or per-anode arrays:
./wct-wires-file.sh summary --with-planes out/<canonical>.json.bz2
./wct-wires-file.sh summary --with-anodes --with-planes out/<canonical>.json.bz2
```

## Validate / reproduce an existing file (no dunecore needed)

The reference side of the reproduction test uses only `wcwires`:

```bash
# existing file -> same WCT pitch corrections the writer applies internally
wcwires -c 4 -o ref.json.bz2  existing-wires.json.bz2

# once generated via `wct-wires-file.sh generate`:
wires-compare.py  out/<canonical>.json.bz2  ref.json.bz2
```

`IDENTICAL` = byte-identical uncompressed JSON; `EQUIVALENT` = within numeric
tolerance; `DIFFER` = real differences (prints the first 50).
