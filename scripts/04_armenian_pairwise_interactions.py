#!/usr/bin/env python3
"""
Reconstructed reproducibility script for the 252 exploratory pairwise
variant-variant interaction tests in the Armenian sample.

Analysis family:
  21 variant pairs × 6 aggression traits × 2 sexes = 252 tests.

For each model:
  outcome ~ age + centered(g1) + centered(g2) + centered(g1)*centered(g2)

HC3 heteroscedasticity-robust covariance estimates are used. The interaction
coefficient, 95% CI, p value, partial R2, genotype-cell counts, and strict
sensitivity result are reported. BH-FDR is applied over the full family of
252 tests.

The script reproduces Supporting Information S4.

Restricted participant-level data must remain outside the public repository.

Required local file:
  data/Armenian_full_analysis.xlsx
  sheet "первичные данные"
"""

from __future__ import annotations

import argparse
import itertools
import re
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


TRAITS = ("PhA", "VerbA", "Ang", "Host", "Proact_ag", "React_ag")
GENES = (
    "AR_CAG",
    "HTR2C",
    "OXTRrs53576",
    "HTR2Ars6311",
    "BDNF rs6265",
    "ERCC5",
    "DRD4eIII",
)
SEX_LABEL = {1: "men", 2: "women"}
STRICT_OUTLIER_EXCEL_ROWS = {69, 72, 125, 181}

EXPECTED = {
    ("men", "HTR2C", "ERCC5", "React_ag"):
        (66, 4.004748385621753, 0.00039787198491945296),
    ("men", "AR_CAG", "ERCC5", "React_ag"):
        (65, 4.688523916078002, 0.0008879671422767653),
    ("women", "BDNF rs6265", "ERCC5", "Proact_ag"):
        (113, 1.7337856743757787, 0.012924774092917165),
}


def bh(values):
    return multipletests(np.asarray(values, dtype=float), method="fdr_bh")[1]


def parse_ar(value):
    if pd.isna(value):
        return None
    nums = [float(x) for x in re.findall(r"\d+(?:\.\d+)?", str(value))]
    return nums or None


def parse_drd4(value):
    if pd.isna(value):
        return None
    nums = [int(x) for x in re.findall(r"\d+", str(value))]
    if len(nums) < 2:
        return None
    return float(max(nums[:2]) >= 7)


def allele_dosage(value, allele):
    if pd.isna(value):
        return None
    s = str(value).upper().replace("/", "").replace("-", "").replace(" ", "")
    if len(s) != 2 or any(ch not in "ACGT" for ch in s):
        return None
    return float(s.count(allele))


def genotype_value(row, gene, sex):
    value = row.get(gene)

    if gene == "AR_CAG":
        alleles = parse_ar(value)
        if not alleles:
            return None
        if sex == 1:
            return float(alleles[0] >= 23)
        if len(alleles) < 2:
            return None
        return float(((alleles[0] + alleles[1]) / 2.0) >= 22)

    if gene == "DRD4eIII":
        return parse_drd4(value)

    if gene == "HTR2C":
        if pd.isna(value):
            return None
        s = str(value).strip().upper()
        if sex == 1:
            if s.startswith("C"):
                return 1.0
            if s.startswith("G"):
                return 0.0
            return None
        return allele_dosage(value, "C")

    if gene == "OXTRrs53576":
        return allele_dosage(value, "A")
    if gene == "HTR2Ars6311":
        return allele_dosage(value, "A")
    if gene == "BDNF rs6265":
        return allele_dosage(value, "A")
    if gene == "ERCC5":
        return allele_dosage(value, "T")

    raise KeyError(gene)


def clean_outcome(value, trait):
    if pd.isna(value):
        return None
    x = float(value)
    if trait == "Ang" and x == 50:
        return None
    return x


def collect(df, sex, gene1, gene2, trait, exclude_strict=False):
    y: List[float] = []
    age: List[float] = []
    g1: List[float] = []
    g2: List[float] = []

    for idx, row in df.iterrows():
        excel_row = int(idx) + 2
        if int(row["SEX"]) != sex:
            continue
        if exclude_strict and excel_row in STRICT_OUTLIER_EXCEL_ROWS:
            continue
        if pd.isna(row["AGE"]):
            continue

        outcome = clean_outcome(row[trait], trait)
        a = genotype_value(row, gene1, sex)
        b = genotype_value(row, gene2, sex)

        if outcome is None or a is None or b is None:
            continue

        y.append(outcome)
        age.append(float(row["AGE"]))
        g1.append(float(a))
        g2.append(float(b))

    return (
        np.asarray(y),
        np.asarray(age),
        np.asarray(g1),
        np.asarray(g2),
    )


def cell_counts(g1, g2):
    counts: Dict[str, int] = {}
    for a, b in zip(g1, g2):
        key = f"{a:g}x{b:g}"
        counts[key] = counts.get(key, 0) + 1
    return counts


def fit_one(df, sex, gene1, gene2, trait, exclude_strict=False):
    y, age, g1, g2 = collect(
        df, sex, gene1, gene2, trait, exclude_strict
    )

    c1 = g1 - g1.mean()
    c2 = g2 - g2.mean()
    interaction = c1 * c2

    X = np.column_stack([np.ones(len(y)), age, c1, c2, interaction])
    robust = sm.OLS(y, X).fit(cov_type="HC3")

    full = sm.OLS(y, X).fit()
    reduced = sm.OLS(
        y, np.column_stack([np.ones(len(y)), age, c1, c2])
    ).fit()
    partial_r2 = (
        (full.rsquared - reduced.rsquared) / (1.0 - reduced.rsquared)
    )

    beta = float(robust.params[4])
    se = float(robust.bse[4])
    counts = cell_counts(g1, g2)
    min_cell = min(counts.values()) if counts else 0

    return {
        "N": int(len(y)),
        "beta": beta,
        "SE": se,
        "CI_low": beta - 1.96 * se,
        "CI_high": beta + 1.96 * se,
        "p": float(robust.pvalues[4]),
        "partial_R2": float(partial_r2),
        "cell_count": len(counts),
        "min_cell_n": int(min_cell),
        "stable_cells": "yes" if min_cell >= 4 else "no",
        "cell_sizes": "; ".join(
            f"{key}: {value}" for key, value in sorted(counts.items())
        ),
    }


def run(df):
    rows: List[Dict] = []

    for sex in (1, 2):
        for gene1, gene2 in itertools.combinations(GENES, 2):
            for trait in TRAITS:
                full = fit_one(
                    df, sex, gene1, gene2, trait, exclude_strict=False
                )
                strict = fit_one(
                    df, sex, gene1, gene2, trait, exclude_strict=True
                )

                rows.append({
                    "Sex": SEX_LABEL[sex],
                    "Gene 1": gene1,
                    "Gene 2": gene2,
                    "Trait": trait,
                    "N": full["N"],
                    "interaction beta": full["beta"],
                    "SE": full["SE"],
                    "95% CI lower": full["CI_low"],
                    "95% CI upper": full["CI_high"],
                    "p": full["p"],
                    "partial R2": full["partial_R2"],
                    "Number of genotype cells": full["cell_count"],
                    "Minimum cell n": full["min_cell_n"],
                    "Stable cells": full["stable_cells"],
                    "p without strict phenodeviants": strict["p"],
                    "beta without strict phenodeviants": strict["beta"],
                    "Cell sizes": full["cell_sizes"],
                })

    out = pd.DataFrame(rows)
    out["q_BH_252"] = bh(out["p"])
    out["q_BH_252_without_strict"] = bh(
        out["p without strict phenodeviants"]
    )
    out = out.sort_values(
        ["p", "Sex", "Gene 1", "Gene 2", "Trait"]
    ).reset_index(drop=True)
    return out


def validate(out):
    if len(out) != 252:
        raise RuntimeError(f"Expected 252 interaction tests, found {len(out)}")

    for key, expected in EXPECTED.items():
        sex, g1, g2, trait = key
        row = out[
            (out["Sex"] == sex)
            & (out["Gene 1"] == g1)
            & (out["Gene 2"] == g2)
            & (out["Trait"] == trait)
        ]
        if len(row) != 1:
            raise RuntimeError(f"Expected one row for {key}")
        row = row.iloc[0]
        exp_n, exp_beta, exp_p = expected
        if int(row["N"]) != exp_n:
            raise RuntimeError(f"{key}: N mismatch")
        if abs(float(row["interaction beta"]) - exp_beta) > 5e-10:
            raise RuntimeError(f"{key}: beta mismatch")
        if abs(float(row["p"]) - exp_p) > 5e-10:
            raise RuntimeError(f"{key}: p mismatch")

    if int((out["q_BH_252"] < 0.05).sum()) != 0:
        raise RuntimeError("Unexpected FDR-confirmed interaction.")

    print("Validation PASSED: frozen S4 interaction targets reproduced.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--armenian",
        type=Path,
        default=Path("data/Armenian_full_analysis.xlsx"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs"),
    )
    args = parser.parse_args()

    df = pd.read_excel(args.armenian, sheet_name="первичные данные")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    out = run(df)
    out.to_csv(args.output_dir / "S4_Interactions_252.csv", index=False)

    validate(out)
    print(f"Outputs written to: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
