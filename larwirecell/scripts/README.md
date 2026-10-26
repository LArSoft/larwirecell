# wcls wires-generation & consistency scripts

Tooling to generate WCT "wires files" from a LArSoft geometry (so they are
consistent-by-construction with the channel map) and to validate them.

| script | what it does |
|--------|--------------|
| `generate-wires.sh` | from a JSON spec → generate FHiCL → run `lar` (WireSchemaFileWriter) → `<era><tech>-<variant>.json.bz2` + `.sidecar.json` provenance |
| `wires-summary.py`  | summarize a wires file: uncompressed-JSON sha256; per-anode & per-plane wire/channel counts, endpoint bounding boxes; per-plane mean/stdev of transverse wire angle, pitch distance, and `dir·X`; mean pitch vector; wire-zero and central-wire centers |
| `wires-compare.py`  | compare two wires files (or summaries): exact hash match else structural diff within tolerance; exit 0 = agree |

## Generate a wires file

```bash
cat > fdvd-v6.json <<'EOF'
{ "era": "fd", "technology": "vd", "variant": "3view30deg-v6",
  "gdml": "/abs/path/dunevd10kt_3view_30deg_v6_refactored_nowires.gdml",
  "detector_name": "dunevd10kt_3view_30deg_v6_refactored",
  "channel_map": "WireReadoutCRUGeom", "sorter": "GeoObjectSorterCRU",
  "grouping": "tpcset", "correction": 4 }
EOF
./generate-wires.sh -i fdvd-v6.json -o out/
```

`grouping`: `tpcset` (default, mapper-faithful), `per-tpc`, or `pair-tpc`
(legacy 2-face wrapped, reproduces the shipped v6 convention).

### Runtime requirement — dunecore

The DUNE channel maps (`DUNEWireReadout`, `WireReadoutCRUGeom`,
`CRPWireReadoutGeom`, `ColdBoxWireReadoutGeom`) and sorters
(`GeoObjectSorter*`) live in **dunecore**.  `lar` needs dunecore's plugins on
`LD_LIBRARY_PATH`/`CET_PLUGIN_PATH` and its fcl on `FHICL_FILE_PATH`.  The
xerosere `wcls` spack env does **not** include dunecore; build/install it (or
point `FHICL_EXTRA` and extend `LD_LIBRARY_PATH` at an existing dunecore) before
running `generate-wires.sh`.

## Validate / reproduce an existing file (no dunecore needed)

The reference side of the reproduction test uses only `wcwires`:

```bash
# existing file -> same WCT pitch corrections the writer applies internally
wcwires -c 4 -o ref.json.bz2  existing-wires.json.bz2
wires-summary.py ref.json.bz2 -o ref.summary.json

# once generated via generate-wires.sh:
wires-compare.py  out/fdvd-3view30deg-v6.json.bz2  ref.json.bz2
```

`IDENTICAL` = byte-identical uncompressed JSON; `EQUIVALENT` = within numeric
tolerance; `DIFFER` = real differences (prints the first 50).
