#!/usr/bin/env python3
"""Summarize a Wire-Cell Toolkit "wires file" (.json / .json.bz2 / .json.gz).

Produces a JSON summary used both for provenance side-cars (see
wct-wires-file.sh) and for comparing two wires files (see wires-compare.py).

The summary contains, as specified for the wcls wires-consistency tooling:

  * a sha256 hash of the *uncompressed* wires JSON bytes,
  * per-anode: ident, per-plane wire/channel counts and totals, endpoint
    bounding box,
  * per-wire-plane: anode/plane ident, wire & channel counts, endpoint
    bounding box, mean/stdev of the transverse (Y-Z) wire angle and of the
    pitch distance, mean/stdev of dot(wire_direction, +X), the mean pitch
    vector, and the center points of "wire zero" and of the central wire.

All lengths are in the native units of the file (WCT base unit = mm); angles
in degrees.  Run after WCT pitch-correction (wcwires -c 4 or the
WireSchemaFileWriter default) so the per-plane angle/pitch variances are small.
"""
import argparse
import bz2
import gzip
import hashlib
import json
import math
import sys


def _open_bytes(path):
    if path.endswith(".bz2"):
        return bz2.open(path, "rb").read()
    if path.endswith(".gz"):
        return gzip.open(path, "rb").read()
    return open(path, "rb").read()


def load_store(path):
    """Return (store_dict, uncompressed_json_bytes)."""
    raw = _open_bytes(path)
    d = json.loads(raw)
    store = d["Store"] if "Store" in d else d
    return store, raw


def _unwrap(lst, key):
    """WCT JSON wraps each element as {"<Key>": {...}}; unwrap to the inner dict."""
    out = []
    for e in lst:
        if isinstance(e, dict) and len(e) == 1 and key in e:
            out.append(e[key])
        else:
            out.append(e)
    return out


def _mean_std(vals):
    n = len(vals)
    if n == 0:
        return (0.0, 0.0)
    m = sum(vals) / n
    if n == 1:
        return (m, 0.0)
    var = sum((v - m) ** 2 for v in vals) / (n - 1)
    return (m, math.sqrt(max(0.0, var)))


def _bbox():
    return [ [math.inf]*3, [-math.inf]*3 ]


def _bbox_add(bb, p):
    for i in range(3):
        if p[i] < bb[0][i]:
            bb[0][i] = p[i]
        if p[i] > bb[1][i]:
            bb[1][i] = p[i]


def _bbox_out(bb):
    if bb[0][0] == math.inf:
        return None
    return {"min": bb[0], "max": bb[1]}


def summarize(store):
    detectors = _unwrap(store["detectors"], "Detector")
    anodes = _unwrap(store["anodes"], "Anode")
    faces = _unwrap(store["faces"], "Face")
    planes = _unwrap(store["planes"], "Plane")
    wires = _unwrap(store["wires"], "Wire")
    points = _unwrap(store["points"], "Point") if "points" in store else None

    def wire_ends(w):
        """Return (tail_xyz, head_xyz) for a wire, handling point-index or inline forms."""
        t, h = w["tail"], w["head"]
        if points is not None and isinstance(t, int):
            pt, ph = points[t], points[h]
            return ([pt["x"], pt["y"], pt["z"]], [ph["x"], ph["y"], ph["z"]])
        # inline Point dicts
        if isinstance(t, dict):
            return ([t["x"], t["y"], t["z"]], [h["x"], h["y"], h["z"]])
        # inline list
        return (list(t), list(h))

    def plane_summary(pl, aident):
        widx = pl["wires"]
        n = len(widx)
        chans = set()
        angs = []        # transverse Y-Z wire angle, degrees
        dotx = []        # dot(wire_dir_hat, +X)
        pitchvecs = []   # 3D pitch vectors between consecutive wire centers
        centers = []     # wire center points, in list (pitch) order
        bb = _bbox()
        for wi in widx:
            w = wires[wi]
            chans.add(w["channel"])
            tail, head = wire_ends(w)
            _bbox_add(bb, tail)
            _bbox_add(bb, head)
            d = [head[i] - tail[i] for i in range(3)]
            L = math.sqrt(sum(c*c for c in d)) or 1.0
            dhat = [c / L for c in d]
            dotx.append(dhat[0])
            angs.append(math.degrees(math.atan2(d[2], d[1])))  # atan2(dz,dy)
            centers.append([0.5*(tail[i]+head[i]) for i in range(3)])
        # pitch vectors: consecutive center differences (list order == WIP order)
        for i in range(1, len(centers)):
            pitchvecs.append([centers[i][j] - centers[i-1][j] for j in range(3)])
        pitch_dists = [math.sqrt(sum(c*c for c in pv)) for pv in pitchvecs]
        ang_m, ang_s = _mean_std(angs)
        pd_m, pd_s = _mean_std(pitch_dists)
        dx_m, dx_s = _mean_std(dotx)
        mean_pitch_vec = [0.0, 0.0, 0.0]
        if pitchvecs:
            for pv in pitchvecs:
                for j in range(3):
                    mean_pitch_vec[j] += pv[j]
            mean_pitch_vec = [c / len(pitchvecs) for c in mean_pitch_vec]
        return {
            "anode": aident,
            "plane": pl["ident"],
            "n_wires": n,
            "n_channels": len(chans),
            "bbox": _bbox_out(bb),
            "wire_angle_deg": {"mean": ang_m, "stdev": ang_s},
            "pitch_dist": {"mean": pd_m, "stdev": pd_s},
            "wiredir_dot_x": {"mean": dx_m, "stdev": dx_s},
            "mean_pitch_vector": mean_pitch_vec,
            "wire_zero_center": centers[0] if centers else None,
            "central_wire_center": centers[n // 2] if centers else None,
        }

    plane_summaries = []
    anode_summaries = []
    tot_wires = 0
    all_chans = set()
    for a in anodes:
        aident = a["ident"]
        abb = _bbox()
        a_nwires = 0
        a_chans = set()
        a_planes = []
        for fidx in a["faces"]:
            f = faces[fidx]
            for pidx in f["planes"]:
                ps = plane_summary(planes[pidx], aident)
                ps["face"] = f["ident"]
                plane_summaries.append(ps)
                a_planes.append({"face": f["ident"], "plane": ps["plane"],
                                 "n_wires": ps["n_wires"], "n_channels": ps["n_channels"]})
                a_nwires += ps["n_wires"]
                if ps["bbox"]:
                    _bbox_add(abb, ps["bbox"]["min"])
                    _bbox_add(abb, ps["bbox"]["max"])
                for pidx2 in planes[pidx]["wires"]:
                    c = wires[pidx2]["channel"]
                    a_chans.add(c)
                    all_chans.add(c)
        tot_wires += a_nwires
        anode_summaries.append({
            "ident": aident,
            "n_faces": len(a["faces"]),
            "planes": a_planes,
            "total_wires": a_nwires,
            "total_channels": len(a_chans),
            "bbox": _bbox_out(abb),
        })

    return {
        "totals": {
            "n_detectors": len(detectors),
            "n_anodes": len(anodes),
            "n_faces": len(faces),
            "n_planes": len(planes),
            "n_wires": len(wires),
            "n_unique_channels": len(all_chans),
        },
        "anodes": anode_summaries,
        "planes": plane_summaries,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("wires", help="wires file (.json/.json.bz2/.json.gz)")
    ap.add_argument("-o", "--output", default="/dev/stdout",
                    help="output summary JSON (default: stdout)")
    ap.add_argument("--with-anodes", action="store_true",
                    help="include the per-anode 'anodes' array (excluded by default)")
    ap.add_argument("--with-planes", action="store_true",
                    help="include the per-plane 'planes' array (excluded by default)")
    args = ap.parse_args(argv)

    store, raw = load_store(args.wires)
    summary = summarize(store)
    if not args.with_anodes:
        summary.pop("anodes", None)
    if not args.with_planes:
        summary.pop("planes", None)
    summary["source_file"] = args.wires
    summary["uncompressed_sha256"] = hashlib.sha256(raw).hexdigest()
    summary["uncompressed_bytes"] = len(raw)

    with open(args.output, "w") as fp:
        json.dump(summary, fp, indent=2)
        fp.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
