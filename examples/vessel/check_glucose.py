#!/usr/bin/env python3
"""Check exported glucose dose, positivity and plasma exclusion (NumPy/h5py)."""
import argparse
import csv
from pathlib import Path

import h5py
import numpy as np


def check(output, rate):
    directory = output / "glucose"
    with (directory / "mass_balance.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError("No glucose diagnostics found")
    previous = -1
    for row in rows:
        iteration = int(row["iteration"])
        if iteration <= previous:
            raise ValueError("Duplicate or out-of-order glucose output")
        previous = iteration
        with h5py.File(directory / f"Glucose.{iteration:012d}.h5", "r") as data:
            concentration = data["Concentration_molecules_per_m3"][:]
            molar = data["Concentration_mol_per_m3"][:]
            plasma = data["PlasmaMask"][:]
            spacing = float(np.asarray(data.attrs["dx_m"]).item())
            time = float(np.asarray(data.attrs["time_s"]).item())
            released = float(np.asarray(data.attrs["released_molecules"]).item())
        expected = rate * time
        total = concentration.sum() * spacing**3
        tolerance = 1e-10 * max(expected, 1e-8)
        if not np.isfinite(concentration).all() or np.any(concentration < 0):
            raise ValueError(f"Invalid concentration at {iteration}")
        if not np.isin(plasma, [0, 1]).all() or np.any(concentration[plasma == 0] != 0):
            raise ValueError(f"Glucose outside plasma at {iteration}")
        if not np.allclose(molar, concentration / 6.02214076e23, rtol=1e-14, atol=0):
            raise ValueError(f"Molar conversion failed at {iteration}")
        values = [released, total, float(row["total_molecules"]),
                  float(row["released_molecules"]), float(row["expected_molecules"])]
        if any(not np.isfinite(value) or abs(value - expected) > tolerance for value in values):
            raise ValueError(f"Dose mismatch at {iteration}: expected {expected}, measured {values}")
        if float(row["excluded_molecules"]) != 0 or abs(float(row["time_s"]) - time) > 1e-15:
            raise ValueError(f"Inconsistent diagnostics at {iteration}")
    print(f"Passed {len(rows)} glucose snapshots; final time={time:g} s, "
          f"dose={total:.12g} expected molecules; excluded concentration is zero.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="Simulation output directory")
    parser.add_argument("--rate", type=float, default=100, help="Configured release rate (molecules/s)")
    args = parser.parse_args()
    if not np.isfinite(args.rate) or args.rate < 0:
        parser.error("rate must be finite and nonnegative")
    check(args.output, args.rate)
