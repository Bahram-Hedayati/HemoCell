#!/usr/bin/env python3
"""Check the supplied short vessel run using HemoCell's per-type CSV output."""
import argparse
import csv
import math
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("output", type=Path)
args = parser.parse_args()
for role, expected in (("RBC", None), ("PLT", None), ("TX", 1), ("RX", 1)):
    files = sorted((args.output / "csv").glob(role + ".*.csv"))
    if len(files) < 2:
        raise SystemExit(f"{role}: need initial and later CSV snapshots")
    samples = []
    ids = None
    for path in files:
        with path.open() as stream:
            rows = list(csv.DictReader(stream))
        if expected is None:
            expected = len(rows)
            if expected == 0:
                raise SystemExit(f"{path}: no {role} cells loaded")
        if len(rows) != expected:
            raise SystemExit(f"{path}: expected {expected} cells, found {len(rows)}")
        for row in rows:
            if not all(math.isfinite(float(value)) for value in row.values()):
                raise SystemExit(f"{path}: non-finite value")
        current_ids = {row["cellId"] for row in rows}
        if ids is not None and ids != current_ids:
            raise SystemExit(f"{path}: cell IDs changed")
        ids = current_ids
        samples.append(rows)
    if role in ("TX", "RX"):
        displacement = [float(samples[-1][0][axis]) - float(samples[0][0][axis])
                        for axis in ("X", "Y", "Z")]
        distance = math.sqrt(sum(value * value for value in displacement))
        if distance <= 0:
            raise SystemExit(f"{role}: no first-to-last displacement detected")
        print(f"{role}: stable ID {next(iter(ids))}; displacement [um] = "
              f"{[value * 1e6 for value in displacement]}")
    print(f"{role}: {expected} cells, {len(files)} finite snapshots")
