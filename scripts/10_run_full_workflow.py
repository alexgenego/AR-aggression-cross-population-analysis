#!/usr/bin/env python3
"""
Run the complete reconstructed analysis workflow in the intended order.

This wrapper assumes that the restricted input workbooks have been placed in
the local gitignored data/ directory using the filenames documented in
scripts/README.md.

Each component script performs its own internal reconciliation checks against
selected frozen targets and exits non-zero if validation fails.

By default all scripts 01-09 are run, including the stochastic LASSO and
bootstrap reconstructions with their fixed seeds.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


SCRIPTS = [
    "01_three_population_shared_loci.py",
    "02_armenian_single_locus_fdr.py",
    "03_armenian_joint_seven_locus.py",
    "04_armenian_pairwise_interactions.py",
    "05_armenian_lasso_reconstruction.py",
    "06_hadza_datoga_nine_variant_block.py",
    "07_variance_dispersion_skewness.py",
    "08_age_models.py",
    "09_power_analysis.py",
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs"),
        help="Shared output directory passed to each analysis script.",
    )
    parser.add_argument(
        "--skip-lasso",
        action="store_true",
        help="Skip the stochastic LASSO reconstruction (script 05).",
    )
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    args.output_dir.mkdir(parents=True, exist_ok=True)

    selected = [
        name for name in SCRIPTS
        if not (args.skip_lasso and name.startswith("05_"))
    ]

    print("Running reconstructed aggression-genetics workflow")
    print(f"Python executable: {sys.executable}")
    print(f"Output directory: {args.output_dir.resolve()}")

    for index, name in enumerate(selected, start=1):
        path = script_dir / name
        if not path.exists():
            raise FileNotFoundError(f"Missing workflow script: {path}")

        command = [
            sys.executable,
            str(path),
            "--output-dir",
            str(args.output_dir),
        ]
        print(f"\n[{index}/{len(selected)}] {name}")
        completed = subprocess.run(command, check=False)
        if completed.returncode != 0:
            raise SystemExit(
                f"Workflow stopped: {name} exited with code "
                f"{completed.returncode}."
            )

    print(
        "\nWorkflow completed successfully. "
        "All component-level reconciliation checks passed."
    )


if __name__ == "__main__":
    main()
