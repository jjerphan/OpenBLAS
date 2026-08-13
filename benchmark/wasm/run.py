#!/usr/bin/env python3
"""Run the frozen WASM OpenBLAS suite under Node and write a CSV."""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path

SUITE: list[tuple[str, int]] = [
    # Level 3
    *[(op, n) for op in ("sgemm", "dgemm") for n in (32, 64, 128, 256, 512)],
    *[(op, n) for op in ("strsm", "dtrsm") for n in (32, 64, 128, 256)],
    *[(op, n) for op in ("ssyrk", "dsyrk") for n in (32, 64, 128, 256)],
    # Level 2
    *[(op, n) for op in ("sgemv", "dgemv", "sgemv_t", "dgemv_t") for n in (64, 256, 1024)],
    # Level 1
    *[(op, n) for op in ("saxpy", "daxpy", "sdot", "ddot", "sscal", "dscal") for n in (10_000, 100_000, 1_000_000)],
]

FIELDS = ["op", "n", "inner", "median_s", "mflops", "min_s", "max_s"]


def main() -> int:
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser()
    p.add_argument("-o", "--output", type=Path, default=here / "results" / "latest.csv")
    p.add_argument("--bench-js", type=Path, default=here / "out" / "bench.js")
    p.add_argument("--node", default=os.environ.get("NODE", "node"))
    p.add_argument("--label", default="")
    args = p.parse_args()

    if not args.bench_js.exists():
        print(f"missing {args.bench_js}; run build.sh first", file=sys.stderr)
        return 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    env = os.environ.copy()
    print("op,n,inner,median_s,mflops,min_s,max_s", flush=True)
    for op, n in SUITE:
        cmd = [args.node, str(args.bench_js), op, str(n)]
        proc = subprocess.run(cmd, check=False, capture_output=True, text=True, env=env)
        line = proc.stdout.strip().splitlines()[-1] if proc.stdout.strip() else ""
        if proc.returncode != 0 or not line:
            print(proc.stderr, file=sys.stderr)
            print(f"FAILED {op} {n} rc={proc.returncode}", file=sys.stderr)
            return proc.returncode or 1
        print(line, flush=True)
        parts = line.split(",")
        rows.append(dict(zip(FIELDS, parts, strict=True)))

    with args.output.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {args.output} ({len(rows)} rows)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
