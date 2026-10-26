#!/bin/bash
# generate-wires.sh -- produce a WCT "wires file" from a LArSoft geometry by
# running the WireSchemaFileWriter art module, plus a provenance side-car.
#
# Given a small JSON spec of the transformation, this:
#   1. derives a canonical name  <era><technology>-<variant>  (e.g. fdvd-v6),
#   2. generates a FHiCL file that configures the LArSoft Geometry + WireReadout
#      (channel map + sorter) and runs the WireSchemaFileWriter analyzer,
#   3. runs `lar` to produce <canonical>.json.bz2 (the wires file),
#   4. writes <canonical>.sidecar.json containing the input spec, the FHiCL text
#      and a full summary of the wires file (hashes, per-anode / per-plane stats).
#
# Input JSON schema (all keys required unless noted):
#   {
#     "era":          "fd" | "pd",          # far detector vs proto-dune
#     "technology":   "vd" | "hd",          # vertical vs horizontal drift
#     "variant":      "3view30deg-v6",      # free-form qualifier
#     "gdml":         "/abs/path/to/detector_nowires.gdml",
#     "detector_name":"dunevd10kt_3view_30deg_v6_refactored",  # drives DUNEWireReadout
#     "channel_map":  "WireReadoutCRUGeom", # provenance (chosen by detector_name)
#     "sorter":       "GeoObjectSorterCRU", # geometry sorter tool_type
#     "grouping":     "tpcset",             # optional: tpcset|per-tpc|pair-tpc (def tpcset)
#     "correction":   4                     # optional: 1..4 (def 4 = pitch)
#   }
#
# Environment (override as needed; defaults assume the xerosere "wcls" env):
#   XERO           repo root (default: autodetected two levels up from this script's repo)
#   WCLS_VIEW      spack view prefix       (default: $XERO/extern/envs/wcls/view)
#   WCT_INSTALL    WCT install w/ correct()(default: $XERO/installs/envs/wcls)
#   LWC_BUILD      larwirecell build dir   (default: $XERO/builds/envs/wcls/larwirecell)
#   FHICL_EXTRA    extra dirs for FHICL_FILE_PATH (e.g. dunecore fcl)  [optional]
#
# NOTE: the DUNE channel maps (DUNEWireReadout / WireReadout*Geom) and sorters
# (GeoObjectSorter*) live in *dunecore*.  dunecore (and its fcl) must be on the
# runtime LD_LIBRARY_PATH / FHICL_FILE_PATH for `lar` to resolve them.  Point
# FHICL_EXTRA / extend LD_LIBRARY_PATH accordingly if dunecore is not in the view.
set -euo pipefail

usage() { sed -n '2,40p' "$0"; exit "${1:-0}"; }

INPUT=""; OUTDIR="."
while [ $# -gt 0 ]; do
  case "$1" in
    -i|--input)  INPUT="$2"; shift 2;;
    -o|--outdir) OUTDIR="$2"; shift 2;;
    -h|--help)   usage 0;;
    *) echo "unknown arg: $1" >&2; usage 1;;
  esac
done
[ -n "$INPUT" ] || { echo "ERROR: --input <spec.json> required" >&2; usage 1; }
[ -f "$INPUT" ] || { echo "ERROR: input not found: $INPUT" >&2; exit 1; }

HERE="$(cd "$(dirname "$0")" && pwd)"
SUMMARY_PY="$HERE/wires-summary.py"

# --- defaults for the xerosere wcls env -----------------------------------
: "${XERO:=$(cd "$HERE/../../../.." && pwd)}"
: "${WCLS_VIEW:=$XERO/extern/envs/wcls/view}"
: "${WCT_INSTALL:=$XERO/installs/envs/wcls}"
: "${LWC_BUILD:=$XERO/builds/envs/wcls/larwirecell}"

# --- parse the input spec --------------------------------------------------
read_json() { python3 -c "import json,sys;print(json.load(open('$INPUT')).get('$1',''))"; }
ERA=$(read_json era); TECH=$(read_json technology); VAR=$(read_json variant)
GDML=$(read_json gdml); DETNAME=$(read_json detector_name)
GDML_SHA1=$(read_json gdml_sha1sum)
SORTER=$(read_json sorter); CMAP=$(read_json channel_map)
PDMAP_TOOL=$(read_json pdmap_tool); PDMAP_FILE=$(read_json pdmap_file)
GROUPING=$(read_json grouping); [ -n "$GROUPING" ] || GROUPING=tpcset
CORR=$(read_json correction);  [ -n "$CORR" ] || CORR=4
for v in ERA TECH VAR GDML DETNAME SORTER; do
  eval "val=\$$v"; [ -n "$val" ] || { echo "ERROR: input missing '$v'" >&2; exit 1; }
done
# Resolve the GDML.  An absolute path is used as-is.  A relative path (including
# a bare basename) is resolved by first checking the current working directory,
# then walking the ':'-separated $GDML_PATH.  The first existing match is taken.
# If the spec provides "gdml_sha1sum", the found file's sha1sum must match it or
# it is a fatal error (if absent, the first filename match is accepted).
resolve_gdml() {
  local g="$1"
  case "$g" in
    /*) [ -f "$g" ] && { printf '%s' "$g"; return 0; } ;;   # absolute
    *)
      [ -f "$PWD/$g" ] && { printf '%s' "$PWD/$g"; return 0; }   # working dir
      local d
      for d in $(printf '%s' "${GDML_PATH:-}" | tr ':' ' '); do
        [ -f "$d/$g" ] && { printf '%s' "$d/$g"; return 0; }     # GDML_PATH entries
      done ;;
  esac
  return 1
}
GDML_RESOLVED="$(resolve_gdml "$GDML")" \
  || { echo "ERROR: gdml '$GDML' not found in CWD or \$GDML_PATH" >&2; exit 1; }
GDML="$GDML_RESOLVED"
if [ -n "$GDML_SHA1" ]; then
  got=$(sha1sum "$GDML" | awk '{print $1}')
  [ "$got" = "$GDML_SHA1" ] || {
    echo "ERROR: gdml sha1 mismatch for $GDML" >&2
    echo "  expected: $GDML_SHA1" >&2
    echo "  found:    $got" >&2
    exit 1; }
  echo ">>> gdml sha1 : verified ($got)"
fi

CANON="${ERA}${TECH}-${VAR}"
mkdir -p "$OUTDIR"
WIRES="$OUTDIR/${CANON}.json.bz2"
FCL="$OUTDIR/${CANON}.fcl"
SIDE="$OUTDIR/${CANON}.sidecar.json"
GDML_DIR="$(cd "$(dirname "$GDML")" && pwd)"
GDML_BASE="$(basename "$GDML")"

# Optional PD (photon-detector) map tool.  DUNEWireReadout requires a PDMapTool
# when it selects a VD/CRP channel map (CRPWireReadoutGeom / ColdBox...), but not
# for the non-wrapping WireReadoutCRUGeom.  Emit the block only when requested.
PDMAP_FCL=""
PDMAP_DIR=""
if [ -n "$PDMAP_TOOL" ]; then
  if [ -n "$PDMAP_FILE" ]; then
    # Resolve the mapping file: absolute as-is, else CWD then $PDMAP_PATH then $GDML_PATH.
    pdres=""
    case "$PDMAP_FILE" in
      /*) [ -f "$PDMAP_FILE" ] && pdres="$PDMAP_FILE" ;;
      *) [ -f "$PWD/$PDMAP_FILE" ] && pdres="$PWD/$PDMAP_FILE"
         for d in $(printf '%s' "${PDMAP_PATH:-}:${GDML_PATH:-}" | tr ':' ' '); do
           [ -z "$pdres" ] && [ -f "$d/$PDMAP_FILE" ] && pdres="$d/$PDMAP_FILE"
         done ;;
    esac
    [ -n "$pdres" ] || { echo "ERROR: pdmap_file '$PDMAP_FILE' not found (set PDMAP_PATH)" >&2; exit 1; }
    PDMAP_DIR="$(cd "$(dirname "$pdres")" && pwd)"
    PDMAP_FCL="    PDMapTool: { tool_type: \"$PDMAP_TOOL\" MappingFile: \"$(basename "$pdres")\" }"
  else
    PDMAP_FCL="    PDMapTool: { tool_type: \"$PDMAP_TOOL\" }"
  fi
fi

# --- generate the FHiCL ----------------------------------------------------
cat > "$FCL" <<EOF
# Auto-generated by generate-wires.sh for canonical name: ${CANON}
# channel map (via detector name): ${CMAP}   sorter: ${SORTER}
process_name: wireswriter

services: {
  # The channel map is chosen by DUNEWireReadout from the geometry Name.
  WireReadout: {
    service_provider: "DUNEWireReadout"
    ChannelsPerOpDet: 1
${PDMAP_FCL}
  }
  GeometryConfigurationWriter: {}
  # The TPC/cryostat sorter (which fixes channel ordering) is on the Geometry
  # service, NOT on WireReadout.
  Geometry: {
    GDML: "${GDML_BASE}"
    Name: "${DETNAME}"
    SurfaceY: 0
    SortingParameters: { tool_type: "${SORTER}" }
  }
}

source: { module_type: EmptyEvent  maxEvents: 1 }

physics: {
  analyzers: {
    wireswriter: {
      module_type: WireSchemaFileWriter
      output:     "${WIRES}"
      grouping:   "${GROUPING}"
      correction: ${CORR}
      validate:   true
    }
  }
  ana:       [ wireswriter ]
  end_paths: [ ana ]
}
EOF

# --- runtime environment ---------------------------------------------------
# Order matters: our WCT install (with correct()) and our larwirecell build
# must precede the view's stock WCT 0.37.0 / larwirecell.
export LD_LIBRARY_PATH="$WCT_INSTALL/lib:$LWC_BUILD/lib:${WCLS_VIEW}/lib:${LD_LIBRARY_PATH:-}"
export CET_PLUGIN_PATH="$LWC_BUILD/lib:$WCT_INSTALL/lib:${WCLS_VIEW}/lib:${CET_PLUGIN_PATH:-}"
export FW_SEARCH_PATH="$GDML_DIR${PDMAP_DIR:+:$PDMAP_DIR}:${FW_SEARCH_PATH:-}"
export FHICL_FILE_PATH="$OUTDIR:${FHICL_EXTRA:-}:${WCLS_VIEW}/fcl:${FHICL_FILE_PATH:-}"

echo ">>> canonical : $CANON"
echo ">>> gdml      : $GDML"
echo ">>> fcl       : $FCL"
echo ">>> running   : lar -c $FCL"
"${WCLS_VIEW}/bin/lar" -c "$FCL" -n 1

[ -f "$WIRES" ] || { echo "ERROR: wires file not produced: $WIRES" >&2; exit 2; }

# --- summary + side-car ----------------------------------------------------
TMP_SUM="$OUTDIR/.${CANON}.summary.json"
python3 "$SUMMARY_PY" "$WIRES" -o "$TMP_SUM"
python3 - "$INPUT" "$FCL" "$TMP_SUM" "$SIDE" "$CANON" "$CMAP" "$SORTER" <<'PY'
import json, sys
inp, fcl, summ, out, canon, cmap, sorter = sys.argv[1:8]
side = {
    "canonical_name": canon,
    "input": json.load(open(inp)),
    "channel_map": cmap,
    "sorter": sorter,
    "fcl_text": open(fcl).read(),
    "summary": json.load(open(summ)),
}
json.dump(side, open(out, "w"), indent=2)
open(out, "a").write("\n")
PY
rm -f "$TMP_SUM"
echo ">>> wires     : $WIRES"
echo ">>> sidecar   : $SIDE"
