"""
Add a ``std_dev`` column to every summary / PMF data file under
``figures/data/distribution_data`` and ``figures/data/pmf_data``.

The standard deviation is computed from the standard error (already present
in the files) as::

    std_dev = sqrt(9) * se = 3 * se

The script targets the per-residue summary files only:

* ``figures/data/distribution_data/{total,monomer_4.5A,multimer_4.5A}/<acid>/summary_<acid>.dat``
* ``figures/data/pmf_data/{total,monomer_4.5A,multimer_4.5A}/<acid>/pmf_<acid>.dat``

Raw per-trajectory files (``trajectory1.dat`` etc.) are left untouched.
Files that already contain a ``std_dev`` column are recomputed in place so
the script is idempotent.

Run from the repository root::

    python add_std_dev.py
"""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent
DATA_ROOT = REPO_ROOT / "figures" / "data"

# (sub-folder under figures/data, glob pattern for the summary files to edit)
TARGETS = [
    ("distribution_data", "summary_*.dat"),
    ("pmf_data", "pmf_*.dat"),
]
SUBSETS = ["total", "monomer_4.5A", "multimer_4.5A"]

# std_dev = sqrt(N) * se, with N=9 (3 trajectories x 3 batches of 200 ns)
N_SAMPLES = 9
SCALE = math.sqrt(N_SAMPLES)

# Output formatting — keep the look of the original files (right-aligned
# fixed-width columns, 6 decimals).
COL_WIDTH = 13
DECIMALS = 6


def format_row(values: list[float]) -> str:
    return "".join(f"{v:>{COL_WIDTH}.{DECIMALS}f}" for v in values)


def format_header(names: list[str]) -> str:
    return "".join(f"{n:>{COL_WIDTH}}" for n in names)


def process_file(path: Path) -> None:
    df = pd.read_csv(path, sep=r"\s+", engine="python")

    if "se" not in df.columns:
        print(f"  skip (no 'se' column): {path.relative_to(REPO_ROOT)}")
        return

    df["std_dev"] = SCALE * df["se"]

    # Preferred column order; keep any extra columns at the end just in case.
    preferred = ["z", "mean", "se", "std_dev"]
    extra = [c for c in df.columns if c not in preferred]
    cols = [c for c in preferred if c in df.columns] + extra
    df = df[cols]

    lines = [format_header(cols)]
    lines.extend(format_row([row[c] for c in cols]) for _, row in df.iterrows())
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    if not DATA_ROOT.is_dir():
        raise SystemExit(f"data directory not found: {DATA_ROOT}")

    total = 0
    for sub, pattern in TARGETS:
        base = DATA_ROOT / sub
        if not base.is_dir():
            print(f"missing: {base.relative_to(REPO_ROOT)} (skipped)")
            continue
        for subset in SUBSETS:
            root = base / subset
            if not root.is_dir():
                print(f"missing: {root.relative_to(REPO_ROOT)} (skipped)")
                continue
            files = sorted(root.glob(f"*/{pattern}"))
            print(f"{root.relative_to(REPO_ROOT)}: {len(files)} file(s)")
            for f in files:
                process_file(f)
                total += 1

    print(f"\nDone. Updated {total} file(s) with std_dev = sqrt({N_SAMPLES}) * se.")


if __name__ == "__main__":
    main()
