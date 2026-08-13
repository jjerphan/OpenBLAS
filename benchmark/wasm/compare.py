#!/usr/bin/env python3
"""Compare two WASM bench CSVs. Speedup > 1 means the new file is faster."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path


def load(path: Path) -> dict[tuple[str, str], dict[str, float]]:
    out: dict[tuple[str, str], dict[str, float]] = {}
    with path.open() as f:
        for row in csv.DictReader(f):
            out[(row["op"], row["n"])] = {
                "median_s": float(row["median_s"]),
                "mflops": float(row["mflops"]),
            }
    return out


def geomean(xs: list[float]) -> float:
    if not xs:
        return float("nan")
    return math.exp(sum(math.log(x) for x in xs) / len(xs))


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("baseline", type=Path)
    p.add_argument("new", type=Path)
    args = p.parse_args()

    base = load(args.baseline)
    new = load(args.new)
    keys = sorted(set(base) & set(new), key=lambda k: (k[0], int(k[1])))
    missing = set(base) ^ set(new)
    if missing:
        print(f"warning: unmatched keys {sorted(missing)}", file=sys.stderr)

    print(f"{'op':<10} {'n':>8} {'base_mflops':>12} {'new_mflops':>12} {'speedup':>8}")
    speedups: list[float] = []
    by_op: dict[str, list[float]] = {}
    for op, n in keys:
        b = base[(op, n)]["mflops"]
        nv = new[(op, n)]["mflops"]
        sp = nv / b if b > 0 else float("nan")
        speedups.append(sp)
        by_op.setdefault(op, []).append(sp)
        print(f"{op:<10} {n:>8} {b:12.4f} {nv:12.4f} {sp:8.3f}")

    print()
    print(f"{'group':<10} {'geomean_speedup':>16}")
    for op, xs in by_op.items():
        print(f"{op:<10} {geomean(xs):16.3f}")
    print(f"{'ALL':<10} {geomean(speedups):16.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
