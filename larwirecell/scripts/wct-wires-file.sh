#!/bin/bash
# wct-wires-file.sh -- produce and summarize WCT "wires files".
#
# Subcommands:
#
#   generate -i SPEC.json -o OUTDIR
#       Run the WireSchemaFileWriter art module (via `lar`) to build a WCT wires
#       file directly from the LArSoft geometry + channel map described by SPEC,
#       and write a provenance side-car.  The side-car holds the input spec, the
#       generated FHiCL text, and the chosen channel-map/sorter -- but NOT a wire
#       summary (use the `summary` subcommand for that).
#       Produces: OUTDIR/<canonical>.{json.bz2, fcl, sidecar.json}
#       where <canonical> = <era><technology>-<variant>.
#
#   summary [--with-anodes] [--with-planes] WIRES [-o OUT.json]
#       Summarize an existing wires file: a sha256 of the uncompressed JSON and
#       the totals.  The per-anode "anodes" array and the per-wire-plane "planes"
#       array (counts, bounding boxes, wire angle / pitch mean+stdev, dir.X, mean
#       pitch vector, wire-0 & central centers) are BOTH EXCLUDED by default; pass
#       --with-anodes and/or --with-planes to include them.  Writes to stdout
#       unless -o is given.
#
# generate SPEC.json schema (all keys required unless noted):
#   {
#     "era":          "fd" | "pd",          # far detector vs proto-dune
#     "technology":   "vd" | "hd",          # vertical vs horizontal drift
#     "variant":      "3view30deg-v6",      # free-form qualifier
#     "gdml":         "detector.gdml",      # abs path, or relative (see below)
#     "gdml_sha1sum": "....",               # optional: pin/verify the GDML
#     "detector_name":"dunevd10kt_3view_30deg_v6_refactored",  # drives DUNEWireReadout
#     "channel_map":  "CRPWireReadoutGeom", # provenance (chosen by detector_name)
#     "sorter":       "GeoObjectSorterCRU60D", # geometry sorter tool_type
#     "pdmap_tool":   "PDVDPDMapAlg",       # optional: required by CRP/ColdBox maps
#     "pdmap_file":   "PDVD_PDS_Mapping_v09162025.json", # optional PD map file
#     "grouping":     "pair-tpc",           # optional: tpcset|per-tpc|pair-tpc (def tpcset)
#     "correction":   4                     # optional: 1..4 (def 4 = pitch; ignored until WCT has correct())
#   }
# A relative "gdml"/"pdmap_file" is resolved against CWD then $GDML_PATH /
# $PDMAP_PATH.  If "gdml_sha1sum" is given, the resolved GDML's sha1 must match.
#
# Environment (override as needed; defaults assume the xerosere "wcls" env):
#   XERO           repo root (default: autodetected from this script's location)
#   WCLS_VIEW      spack view prefix       (default: $XERO/extern/envs/wcls/view)
#   WCT_INSTALL    WCT install             (default: $XERO/installs/envs/wcls)
#   LWC_BUILD      larwirecell build dir   (default: $XERO/builds/envs/wcls/larwirecell)
#   GDML_PATH      ':'-sep dirs to find a relative gdml
#   PDMAP_PATH     ':'-sep dirs to find a relative pdmap_file
#   FHICL_EXTRA    extra dirs for FHICL_FILE_PATH (e.g. dunecore fcl)  [optional]
#
# NOTE: the DUNE channel maps (DUNEWireReadout / WireReadout*Geom) and sorters
# (GeoObjectSorter*) live in *dunecore*.  dunecore must be on the runtime
# LD_LIBRARY_PATH / CET_PLUGIN_PATH / FHICL_FILE_PATH for `lar` to resolve them.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
SUMMARY_PY="$HERE/wires-summary.py"

usage() { sed -n '2,52p' "$0"; exit "${1:-0}"; }

#============================================================================
cmd_generate() {
  local INPUT="" OUTDIR="."
  while [ $# -gt 0 ]; do
    case "$1" in
      -i|--input)  INPUT="$2"; shift 2;;
      -o|--outdir) OUTDIR="$2"; shift 2;;
      -h|--help)   usage 0;;
      *) echo "unknown arg: $1" >&2; usage 1;;
    esac
  done
  [ -n "$INPUT" ] || { echo "ERROR: generate --input <spec.json> required" >&2; usage 1; }
  [ -f "$INPUT" ] || { echo "ERROR: input not found: $INPUT" >&2; exit 1; }

  # --- defaults for the xerosere wcls env ----------------------------------
  : "${XERO:=$(cd "$HERE/../../../.." && pwd)}"
  : "${WCLS_VIEW:=$XERO/extern/envs/wcls/view}"
  : "${WCT_INSTALL:=$XERO/installs/envs/wcls}"
  : "${LWC_BUILD:=$XERO/builds/envs/wcls/larwirecell}"

  # --- parse the input spec ------------------------------------------------
  read_json() { python3 -c "import json,sys;print(json.load(open('$INPUT')).get('$1',''))"; }
  local ERA TECH VAR GDML DETNAME GDML_SHA1 SORTER CMAP PDMAP_TOOL PDMAP_FILE GROUPING CORR
  ERA=$(read_json era); TECH=$(read_json technology); VAR=$(read_json variant)
  GDML=$(read_json gdml); DETNAME=$(read_json detector_name)
  GDML_SHA1=$(read_json gdml_sha1sum)
  SORTER=$(read_json sorter); CMAP=$(read_json channel_map)
  PDMAP_TOOL=$(read_json pdmap_tool); PDMAP_FILE=$(read_json pdmap_file)
  GROUPING=$(read_json grouping); [ -n "$GROUPING" ] || GROUPING=tpcset
  CORR=$(read_json correction);  [ -n "$CORR" ] || CORR=4
  local v val
  for v in ERA TECH VAR GDML DETNAME SORTER; do
    eval "val=\$$v"; [ -n "$val" ] || { echo "ERROR: input missing '$v'" >&2; exit 1; }
  done

  # Resolve the GDML: absolute as-is, else CWD then $GDML_PATH.
  local g d
  g="$GDML"; GDML=""
  case "$g" in
    /*) [ -f "$g" ] && GDML="$g" ;;
    *)  [ -f "$PWD/$g" ] && GDML="$PWD/$g"
        for d in $(printf '%s' "${GDML_PATH:-}" | tr ':' ' '); do
          [ -z "$GDML" ] && [ -f "$d/$g" ] && GDML="$d/$g"
        done ;;
  esac
  [ -n "$GDML" ] || { echo "ERROR: gdml '$g' not found in CWD or \$GDML_PATH" >&2; exit 1; }
  if [ -n "$GDML_SHA1" ]; then
    local got; got=$(sha1sum "$GDML" | awk '{print $1}')
    [ "$got" = "$GDML_SHA1" ] || {
      echo "ERROR: gdml sha1 mismatch for $GDML" >&2
      echo "  expected: $GDML_SHA1" >&2
      echo "  found:    $got" >&2; exit 1; }
    echo ">>> gdml sha1 : verified ($got)"
  fi

  local CANON WIRES FCL SIDE GDML_DIR GDML_BASE
  CANON="${ERA}${TECH}-${VAR}"
  mkdir -p "$OUTDIR"
  WIRES="$OUTDIR/${CANON}.json.bz2"
  FCL="$OUTDIR/${CANON}.fcl"
  SIDE="$OUTDIR/${CANON}.sidecar.json"
  GDML_DIR="$(cd "$(dirname "$GDML")" && pwd)"
  GDML_BASE="$(basename "$GDML")"

  # Optional PD (photon-detector) map tool: required by CRP/ColdBox channel maps,
  # not by the non-wrapping WireReadoutCRUGeom.
  local PDMAP_FCL="" PDMAP_DIR="" pdres
  if [ -n "$PDMAP_TOOL" ]; then
    if [ -n "$PDMAP_FILE" ]; then
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

  # --- generate the FHiCL --------------------------------------------------
  cat > "$FCL" <<EOF
# Auto-generated by 'wct-wires-file.sh generate' for canonical name: ${CANON}
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

  # --- runtime environment -------------------------------------------------
  # Order matters: our WCT install and our larwirecell build
  # must precede the view's stock WCT / larwirecell.
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

  # --- side-car (input spec + FHiCL; NO summary) ---------------------------
  python3 - "$INPUT" "$FCL" "$SIDE" "$CANON" "$CMAP" "$SORTER" <<'PY'
import json, sys
inp, fcl, out, canon, cmap, sorter = sys.argv[1:7]
side = {
    "canonical_name": canon,
    "input": json.load(open(inp)),
    "channel_map": cmap,
    "sorter": sorter,
    "fcl_text": open(fcl).read(),
}
json.dump(side, open(out, "w"), indent=2)
open(out, "a").write("\n")
PY
  echo ">>> wires     : $WIRES"
  echo ">>> sidecar   : $SIDE"
}

#============================================================================
cmd_summary() {
  local opts="" out="/dev/stdout" wires=""
  while [ $# -gt 0 ]; do
    case "$1" in
      --with-anodes) opts="$opts --with-anodes"; shift;;
      --with-planes) opts="$opts --with-planes"; shift;;
      -o|--output)   out="$2"; shift 2;;
      -h|--help)     usage 0;;
      -*) echo "unknown option: $1" >&2; usage 1;;
      *)  wires="$1"; shift;;
    esac
  done
  [ -n "$wires" ] || { echo "ERROR: summary requires a wires file" >&2; usage 1; }
  [ -f "$wires" ] || { echo "ERROR: wires file not found: $wires" >&2; exit 1; }
  # shellcheck disable=SC2086
  python3 "$SUMMARY_PY" $opts "$wires" -o "$out"
}

#============================================================================
sub="${1:-}"; [ $# -gt 0 ] && shift || true
case "$sub" in
  generate) cmd_generate "$@";;
  summary)  cmd_summary "$@";;
  -h|--help|help|"") usage 0;;
  *) echo "unknown subcommand: '$sub' (want 'generate' or 'summary')" >&2; usage 1;;
esac
