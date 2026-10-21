#!/usr/bin/env python3
"""Summary statistics of a ContAngle run, for regression references (compare with `git diff`).

Usage: post_process.py [caseDir=subvolume] [logFile=AllRunContAngle.log] > ref/tutorial_subvolume_stats.txt

Prints `key value` lines: the contact angles reported in the log, and the number, mean, standard deviation, min and max
of the non-zero values (points on the interface) of the contact angles (contactAngles.txt) and the curvature (Kc.txt).
"""
import math
import re
import sys
from pathlib import Path


def read_vtk_scalars(fil, name):
    """values of the VTK point-data array `name` (header line: `name 1 N float`)"""
    toks = Path(fil).read_text().split()
    for i, t in enumerate(toks):
        if t == name and toks[i + 1] == "1" and toks[i + 3] == "float":
            n = int(toks[i + 2])
            return [float(v) for v in toks[i + 4 : i + 4 + n]]
    raise ValueError(f"{name} not found in {fil}")


def stats(vals):
    n = len(vals)
    mean = sum(vals) / n
    std = math.sqrt(sum((v - mean) ** 2 for v in vals) / n)
    return n, mean, std, min(vals), max(vals)


def main():
    case = Path(sys.argv[1] if len(sys.argv) > 1 else "subvolume")
    log = Path(sys.argv[2] if len(sys.argv) > 2 else "AllRunContAngle.log")
    out = []
    text = log.read_text(errors="ignore")
    for key in ("ContactAngMax", "ContactAngMin", "ContactAngAvg"):
        m = re.findall(rf"{key}:\s*(\S+)", text)
        out.append((key, m[-1]))
    ca = read_vtk_scalars(case / "contactAngles.txt", "CA")
    kc = read_vtk_scalars(case / "Kc.txt", "Kc")
    out.append(("nPoints", len(ca)))
    for name, vals in (("CA", [v for v in ca if v != 0]), ("Kc", [v for v in kc if v != 0])):
        for key, v in zip(("count", "mean", "std", "min", "max"), stats(vals)):
            out.append((f"{name}_{key}", v if isinstance(v, int) else f"{v:.6g}"))
    print("\n".join(f"{k} {v}" for k, v in out))


if __name__ == "__main__":
    main()
