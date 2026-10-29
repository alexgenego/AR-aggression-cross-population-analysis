#!/usr/bin/env python3
"""
Reconstructed reproducibility script for the Armenian joint seven-variant models.

This script reproduces the joint multilocus analysis summarized in Supporting
Information S3a for the manuscript:

"Androgen Receptor CAG Repeat Length and Aggressive Behavior Across Armenian,
Hadza, and Datoga Men: A Cross-Population Genetic Association Analysis"

IMPORTANT
---------
This is a reconstructed implementation based on the frozen analytical workbook
and final statistical specification. It is not claimed to be the original
historical source code.

Restricted participant-level data remain outside the public repository.

Required local file
-------------------
data/Armenian_full_analysis.xlsx
  sheet: "первичные данные"

Model
-----
Within each sex and each of six aggression traits:
  outcome ~ age + AR_CAG + HTR2C + OXTR + HTR2A + BDNF + ERCC5 + DRD4

All seven primary genetic codings are entered simultaneously. HC3-robust
covariance estimates are used. Primary and strict sensitivity models are
reported; the strict sensitivity set excludes Excel rows 69, 72, 125, and 181.
Ang=50 in Excel row 79 is treated as missing.

Outputs
-------
outputs/S3a_Seven_Locus.csv
outputs/S3a_Model_Diagnostics.csv
"""

from __future__ import annotations

import argparse
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

CODING_LABEL = {
    "AR_CAG": "AR: high category vs low category",
    "HTR2C": "HTR2C: C-allele dosage",
    "OXTRrs53576": "OXTR: A-allele dosage",
    "HTR2Ars6311": "HTR2A: A-allele dosage",
    "BDNF rs6265": "BDNF: A-allele dosage",
    "ERCC5": "ERCC5: T-allele dosage",
    "DRD4eIII": "DRD4: L+ vs SS",
}

EXPECTED = {
    ("men", "PhA", "AR_CAG"): (-2.5131682126128836, 0.004718679635161784, 64),
    ("women", "VerbA", "BDNF rs6265"): (-2.0118344537516264, 0.012407306205232102, 107),
    ("men", "Host", "HTR2Ars6311"): (2.3575513585288244, 0.025478883817329125, 64),
    ("women", "Proact_ag", "AR_CAG"): (-1.4069253734132459, 0.0324892679093122, 106),
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


def collect_complete_cases(df, sex, trait, exclude_strict):
    y: List[float] = []
    age: List[float] = []
    genotypes: List[List[float]] = []

    for idx, row in df.iterrows():
        excel_row = int(idx) + 2
        if int(row["SEX"]) != sex:
            continue
        if exclude_strict and excel_row in STRICT_OUTLIER_EXCEL_ROWS:
            continue
        if pd.isna(row["AGE"]):
            continue

        outcome = clean_outcome(row[trait], trait)
        if outcome is None:
            continue

        g = [genotype_value(row, gene, sex) for gene in GENES]
        if any(v is None for v in g):
            continue

        y.append(outcome)
        age.append(float(row["AGE"]))
        genotypes.append([float(v) for v in g])

    return np.asarray(y), np.asarray(age), np.asarray(genotypes)


def fit_joint(df, sex, trait, exclude_strict):
    y, age, G = collect_complete_cases(df, sex, trait, exclude_strict)
    X = np.column_stack([np.ones(len(y)), G, age])

    robust = sm.OLS(y, X).fit(cov_type="HC3")
    ols = sm.OLS(y, X).fit()

    standardized_beta = (
        robust.params[1:1 + len(GENES)]
        * G.std(axis=0, ddof=1)
        / y.std(ddof=1)
    )

    return {
        "N": len(y),
        "robust": robust,
        "ols": ols,
        "standardized_beta": standardized_beta,
    }


def run(df):
    coefficient_rows: List[Dict] = []
    diagnostic_rows: List[Dict] = []

    for sex in (1, 2):
        for trait in TRAITS:
            full = fit_joint(df, sex, trait, exclude_strict=False)
            strict = fit_joint(df, sex, trait, exclude_strict=True)

            m = full["robust"]
            ms = strict["robust"]

            for j, gene in enumerate(GENES):
                beta = float(m.params[1 + j])
                se = float(m.bse[1 + j])
                coefficient_rows.append({
                    "Sex": SEX_LABEL[sex],
                    "Trait": trait,
                    "Gene": gene,
                    "Coding": CODING_LABEL[gene],
                    "N": int(full["N"]),
                    "beta controlling other variants": beta,
                    "SE": se,
                    "95% CI lower": beta - 1.96 * se,
                    "95% CI upper": beta + 1.96 * se,
                    "standardized beta": float(full["standardized_beta"][j]),
                    "p": float(m.pvalues[1 + j]),
                    "N without strict phenodeviants": int(strict["N"]),
                    "beta without strict phenodeviants": float(ms.params[1 + j]),
                    "p without strict phenodeviants": float(ms.pvalues[1 + j]),
                })

            diagnostic_rows.append({
                "Sex": SEX_LABEL[sex],
                "Trait": trait,
                "N": int(full["N"]),
                "R2": float(full["ols"].rsquared),
                "Adjusted R2": float(full["ols"].rsquared_adj),
                "N without strict phenodeviants": int(strict["N"]),
                "Adjusted R2 without strict phenodeviants": float(strict["ols"].rsquared_adj),
            })

    coef = pd.DataFrame(coefficient_rows)
    coef["q_BH"] = bh(coef["p"])
    coef["q_BH without strict"] = bh(coef["p without strict phenodeviants"])
    coef = coef.sort_values(["p", "Sex", "Trait", "Gene"]).reset_index(drop=True)

    diag = pd.DataFrame(diagnostic_rows)

    return coef, diag


def validate(coef):
    for key, (exp_beta, exp_p, exp_n) in EXPECTED.items():
        sex, trait, gene = key
        row = coef[
            (coef["Sex"] == sex)
            & (coef["Trait"] == trait)
            & (coef["Gene"] == gene)
        ]
        if len(row) != 1:
            raise RuntimeError(f"Expected exactly one row for {key}")
        row = row.iloc[0]
        if int(row["N"]) != exp_n:
            raise RuntimeError(f"{key}: N mismatch")
        if abs(float(row["beta controlling other variants"]) - exp_beta) > 5e-10:
            raise RuntimeError(f"{key}: beta mismatch")
        if abs(float(row["p"]) - exp_p) > 5e-10:
            raise RuntimeError(f"{key}: p mismatch")

    if int((coef["q_BH"] < 0.05).sum()) != 0:
        raise RuntimeError("Unexpected FDR-confirmed coefficient in joint model.")

    print("Validation PASSED: frozen S3a coefficient targets reproduced.")


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

    coef, diag = run(df)
    coef.to_csv(args.output_dir / "S3a_Seven_Locus.csv", index=False)
    diag.to_csv(args.output_dir / "S3a_Model_Diagnostics.csv", index=False)

    validate(coef)
    print(f"Outputs written to: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
