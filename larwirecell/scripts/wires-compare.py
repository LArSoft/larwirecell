#!/usr/bin/env python3
"""Compare two WCT wires files (or two summaries from wires-summary.py).

Intended use: validate that a freshly generated wires file (from
WireSchemaFileWriter) matches an existing "golden" file after both have been put
through the same WCT corrections.  Reports an exact uncompressed-JSON hash match
and, failing that, a structural diff with numeric tolerances.

Exit status 0 if the two agree (hash-identical OR all structural checks within
tolerance), 1 otherwise.
"""
import argparse
import json
import math
import sys

import importlib.util
import os

# import sibling wires-summary.py (hyphenated name -> load by path)
_spec = importlib.util.spec_from_file_location(
    "wires_summary", os.path.join(os.path.dirname(__file__), "wires-summary.py"))
_wsum = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_wsum)


def get_summary(path):
    """Accept either a wires file or a summary JSON; return a summary dict."""
    if path.endswith((".json.bz2", ".json.gz")):
        store, raw = _wsum.load_store(path)
        s = _wsum.summarize(store)
        s["uncompressed_sha256"] = __import__("hashlib").sha256(raw).hexdigest()
        return s
    d = json.load(open(path))
    return d.get("summary", d)  # allow a side-car file


def _close(a, b, atol, rtol):
    if a is None or b is None:
        return a == b
    return abs(a - b) <= atol + rtol * max(abs(a), abs(b))


def _vclose(a, b, atol, rtol):
    if a is None or b is None:
        return a == b
    return all(_close(x, y, atol, rtol) for x, y in zip(a, b))


def compare(sa, sb, atol, rtol):
    diffs = []
    # exact hash
    ha, hb = sa.get("uncompressed_sha256"), sb.get("uncompressed_sha256")
    hash_match = (ha is not None and ha == hb)

    if sa.get("totals") != sb.get("totals"):
        diffs.append(("totals", sa.get("totals"), sb.get("totals")))

    pa = {(p["anode"], p.get("face", 0), p["plane"]): p for p in sa.get("planes", [])}
    pb = {(p["anode"], p.get("face", 0), p["plane"]): p for p in sb.get("planes", [])}
    if set(pa) != set(pb):
        diffs.append(("plane_set", sorted(set(pa) - set(pb)), sorted(set(pb) - set(pa))))

    for k in sorted(set(pa) & set(pb)):
        A, B = pa[k], pb[k]
        for fld in ("n_wires", "n_channels"):
            if A[fld] != B[fld]:
                diffs.append((f"plane{k}.{fld}", A[fld], B[fld]))
        for fld in ("wire_angle_deg", "pitch_dist", "wiredir_dot_x"):
            if not _close(A[fld]["mean"], B[fld]["mean"], atol, rtol):
                diffs.append((f"plane{k}.{fld}.mean", A[fld]["mean"], B[fld]["mean"]))
        for fld in ("wire_zero_center", "central_wire_center", "mean_pitch_vector"):
            if not _vclose(A.get(fld), B.get(fld), atol, rtol):
                diffs.append((f"plane{k}.{fld}", A.get(fld), B.get(fld)))
        if A.get("bbox") and B.get("bbox"):
            if not (_vclose(A["bbox"]["min"], B["bbox"]["min"], atol, rtol) and
                    _vclose(A["bbox"]["max"], B["bbox"]["max"], atol, rtol)):
                diffs.append((f"plane{k}.bbox", A["bbox"], B["bbox"]))
    return hash_match, diffs


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a", help="wires file or summary/sidecar JSON (A)")
    ap.add_argument("b", help="wires file or summary/sidecar JSON (B)")
    ap.add_argument("--atol", type=float, default=1e-6, help="absolute tolerance (file units / deg)")
    ap.add_argument("--rtol", type=float, default=1e-6, help="relative tolerance")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args(argv)

    sa, sb = get_summary(args.a), get_summary(args.b)
    hash_match, diffs = compare(sa, sb, args.atol, args.rtol)

    if hash_match:
        if not args.quiet:
            print("IDENTICAL: uncompressed-JSON sha256 match")
        return 0
    if not diffs:
        if not args.quiet:
            print("EQUIVALENT: hashes differ but all structural checks within tolerance "
                  f"(atol={args.atol}, rtol={args.rtol})")
        return 0
    def _short(v, n=160):
        s = repr(v)
        return s if len(s) <= n else s[:n] + f"... (+{len(s)-n} chars)"
    print(f"DIFFER: {len(diffs)} difference(s):")
    for name, va, vb in diffs[:50]:
        print(f"  {name}:\n    A={_short(va)}\n    B={_short(vb)}")
    if len(diffs) > 50:
        print(f"  ... and {len(diffs)-50} more")
    return 1


if __name__ == "__main__":
    sys.exit(main())
