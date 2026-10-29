#!/usr/bin/env python3
"""
Reproduce the pooled Hadza-Datoga nine-variant analyses summarized in S5a-S5c.

Required local input (restricted; do not commit):
  data/Hadza_Datoga_2012_cleaned_males.xlsx
  sheet: "Объединенные_мужчины"

Outputs:
  outputs/S5a_HD_Genetic_Block.csv
  outputs/S5b_HD_Coefficients.csv
  outputs/S5c_HD_Interaction_Block.csv
  outputs/S5c_HD_Individual_Interactions.csv

The source workbook is the cleaned male archive. Its documented cleaning rules
already convert isolated sentinel value 21 to missing while retaining true
AR-21 calls, and convert failed genotypes to missing.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

TRAITS = {
    "physag": "Physical aggression",
    "verrbag": "Verbal aggression",
    "anger": "Anger",
    "hostil": "Hostility",
}

MARKERS = (
    "DRD4P",
    "SERTPr",
    "AR",
    "SERT_In2",
    "DRD4E3",
    "MAOA",
    "DAT1_VNTR",
    "DRD2Taq1",
    "COMT",
)

LABELS_CONT = (
    "DRD4 promoter: S-allele dosage",
    "SERT promoter: low-expression allele dosage",
    "AR: CAG repeat count",
    "SERT intron 2: 12R dosage",
    "DRD4 exon III: L+ versus SS",
    "MAOA repeat count",
    "DAT1: 10R dosage",
    "DRD2: A1 dosage",
    "COMT: Met dosage",
)
LABELS_CAT = list(LABELS_CONT)
LABELS_CAT[2] = "AR: >=23 versus <=22 repeats"

CONTINUOUS_COLS_CONT = (0, 1, 2, 3, 5, 6, 7, 8)
CONTINUOUS_COLS_CAT = (0, 1, 3, 5, 6, 7, 8)

EXPECTED = {
    ("continuous", "physag"): {
        "N": 160,
        "R2": 0.14477606979157498,
        "partial": 0.06044459803252651,
        "block_p": 0.3484721425017169,
        "population_beta": 2.758584543247182,
        "population_p": 0.006428287479827859,
    },
    ("category", "anger"): {
        "N": 154,
        "R2": 0.2954929840163031,
        "partial": 0.1677736187048707,
        "block_p": 0.0001459187155080442,
        "population_beta": 3.7240616563234763,
        "population_p": 0.000013507792783475899,
    },
}


def bh(values):
    return multipletests(np.asarray(values, dtype=float), method="fdr_bh")[1]


def parse_drd4p(value) -> Optional[float]:
    if pd.isna(value):
        return None
    return {"LL": 0.0, "LS": 1.0, "SL": 1.0, "SS": 2.0}.get(
        str(value).strip().upper()
    )


def parse_sertpr(value) -> Optional[float]:
    if pd.isna(value):
        return None
    text = str(value).strip().replace(" ", "/").replace("//", "/")
    parts = [x for x in text.split("/") if x]
    if len(parts) != 2:
        return None
    normalized = []
    mapping = {"s": "S", "sa": "Sa", "lg": "Lg", "la": "La", "l": "L", "xl": "xL"}
    for part in parts:
        token = mapping.get(part.strip().lower())
        if token is None:
            return None
        normalized.append(token)
    if any(token in {"L", "xL"} for token in normalized):
        return None
    low = {"S", "Sa", "Lg"}
    return float(sum(token in low for token in normalized))


def parse_ar(value) -> Optional[float]:
    if pd.isna(value):
        return None
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    return x if 10 <= x <= 40 else None


def parse_repeat_dosage(value, target: int) -> Optional[float]:
    if pd.isna(value):
        return None
    numbers = [int(x) for x in re.findall(r"\d+", str(value))]
    if len(numbers) < 2:
        return None
    return float(sum(x == target for x in numbers[:2]))


def parse_drd4e3(value) -> Optional[float]:
    if pd.isna(value):
        return None
    numbers = [int(x) for x in re.findall(r"\d+", str(value))]
    if len(numbers) < 2:
        return None
    return float(max(numbers[:2]) >= 7)


def parse_maoa(value) -> Optional[float]:
    if pd.isna(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def parse_drd2(value) -> Optional[float]:
    if pd.isna(value):
        return None
    return {
        "A2/A2": 0.0,
        "A1/A2": 1.0,
        "A2/A1": 1.0,
        "A1/A1": 2.0,
    }.get(str(value).strip().upper())


def parse_comt(value) -> Optional[float]:
    if pd.isna(value):
        return None
    return {
        "V/V": 0.0,
        "V/M": 1.0,
        "M/V": 1.0,
        "M/M": 2.0,
    }.get(str(value).strip().upper())


def encode_row(row: pd.Series, ar_mode: str) -> Optional[np.ndarray]:
    values = [
        parse_drd4p(row["DRD4P"]),
        parse_sertpr(row["SERTPr"]),
        parse_ar(row["AR"]),
        parse_repeat_dosage(row["SERT_In2"], 12),
        parse_drd4e3(row["DRD4E3"]),
        parse_maoa(row["MAOA"]),
        parse_repeat_dosage(row["DAT1_VNTR"], 10),
        parse_drd2(row["DRD2Taq1"]),
        parse_comt(row["COMT"]),
    ]
    if any(value is None for value in values):
        return None
    if ar_mode == "category":
        values[2] = float(values[2] >= 23)
    return np.asarray(values, dtype=float)


def collect(df: pd.DataFrame, trait: str, ar_mode: str):
    rows: List[Tuple[str, float, float, np.ndarray]] = []
    for _, row in df.iterrows():
        if int(row["sex"]) != 1:
            continue
        if pd.isna(row["age"]) or pd.isna(row[trait]):
            continue
        g = encode_row(row, ar_mode)
        if g is None:
            continue
        population = "Hadza" if str(row["population"]).strip() == "Хадза" else "Datoga"
        rows.append((population, float(row["age"]), float(row[trait]), g))
    return rows


def global_predictors(rows, ar_mode: str):
    G = np.vstack([row[3] for row in rows])
    cols = CONTINUOUS_COLS_CONT if ar_mode == "continuous" else CONTINUOUS_COLS_CAT
    out = G.copy()
    for j in cols:
        out[:, j] = (G[:, j] - G[:, j].mean()) / G[:, j].std(ddof=1)
    return out


def within_population_standardization(rows, ar_mode: str):
    pops = np.asarray([row[0] for row in rows])
    age = np.asarray([row[1] for row in rows], dtype=float)
    y = np.asarray([row[2] for row in rows], dtype=float)
    G = np.vstack([row[3] for row in rows])
    cols = CONTINUOUS_COLS_CONT if ar_mode == "continuous" else CONTINUOUS_COLS_CAT

    age_z = age.copy()
    y_z = y.copy()
    G_z = G.copy()
    for pop in ("Hadza", "Datoga"):
        mask = pops == pop
        age_z[mask] = (age[mask] - age[mask].mean()) / age[mask].std(ddof=1)
        y_z[mask] = (y[mask] - y[mask].mean()) / y[mask].std(ddof=1)
        for j in cols:
            G_z[mask, j] = (G[mask, j] - G[mask, j].mean()) / G[mask, j].std(ddof=1)
    return pops, age_z, y_z, G_z


def fit_primary(rows, ar_mode: str):
    pops = np.asarray([row[0] for row in rows])
    y = np.asarray([row[2] for row in rows], dtype=float)
    age = np.asarray([row[1] for row in rows], dtype=float)
    G = global_predictors(rows, ar_mode)
    pop_d = (pops == "Datoga").astype(float)

    X = np.column_stack([np.ones(len(y)), age, pop_d, G])
    robust = sm.OLS(y, X).fit(cov_type="HC3")
    ols = sm.OLS(y, X).fit()

    X_reduced = np.column_stack([np.ones(len(y)), age, pop_d])
    reduced = sm.OLS(y, X_reduced).fit()
    partial_r2 = (ols.rsquared - reduced.rsquared) / (1.0 - reduced.rsquared)

    R = np.zeros((9, X.shape[1]))
    for j in range(9):
        R[j, 3 + j] = 1.0
    block_test = robust.wald_test(R, scalar=True)

    return robust, ols, reduced, float(partial_r2), float(block_test.statistic), float(block_test.pvalue), G, pops


def fit_interactions(rows, ar_mode: str):
    pops, age_z, y_z, G_z = within_population_standardization(rows, ar_mode)
    pop_d = (pops == "Datoga").astype(float)

    X_base = np.column_stack([np.ones(len(y_z)), age_z, pop_d, G_z])
    X_full = np.column_stack([X_base, G_z * pop_d[:, None]])

    base = sm.OLS(y_z, X_base).fit(cov_type="HC3")
    full = sm.OLS(y_z, X_full).fit(cov_type="HC3")
    partial_r2 = (full.rsquared - base.rsquared) / (1.0 - base.rsquared)

    R = np.zeros((9, X_full.shape[1]))
    for j in range(9):
        R[j, 12 + j] = 1.0
    block_test = full.wald_test(R, scalar=True)

    return base, full, float(partial_r2), float(block_test.statistic), float(block_test.pvalue), G_z, pops


def run(df: pd.DataFrame):
    block_rows: List[Dict] = []
    coefficient_rows: List[Dict] = []
    interaction_block_rows: List[Dict] = []
    individual_interaction_rows: List[Dict] = []

    for ar_mode in ("continuous", "category"):
        labels = LABELS_CONT if ar_mode == "continuous" else LABELS_CAT
        for trait, trait_label in TRAITS.items():
            rows = collect(df, trait, ar_mode)
            robust, ols, reduced, partial, wald, block_p, G, pops = fit_primary(rows, ar_mode)

            pop_ci = robust.conf_int()[2]
            block_rows.append({
                "AR coding": ar_mode,
                "Trait": trait_label,
                "N": len(rows),
                "Hadza N": int(np.sum(pops == "Hadza")),
                "Datoga N": int(np.sum(pops == "Datoga")),
                "R2 full model": float(ols.rsquared),
                "Adjusted R2": float(ols.rsquared_adj),
                "R2 age+population": float(reduced.rsquared),
                "Partial R2 genetic block": partial,
                "Wald chi2 genetic block": wald,
                "df": 9,
                "Block p": block_p,
                "Datoga-Hadza effect, points": float(robust.params[2]),
                "Population CI low": float(pop_ci[0]),
                "Population CI high": float(pop_ci[1]),
                "Population p": float(robust.pvalues[2]),
            })

            for j, label in enumerate(labels):
                beta = float(robust.params[3 + j])
                se = float(robust.bse[3 + j])
                coefficient_rows.append({
                    "AR coding": ar_mode,
                    "Trait": trait_label,
                    "Marker": label,
                    "N": len(rows),
                    "beta, points": beta,
                    "SE HC3": se,
                    "95% CI lower": beta - 1.96 * se,
                    "95% CI upper": beta + 1.96 * se,
                    "p": float(robust.pvalues[3 + j]),
                })

            base, full, ipartial, iwald, ip, Gz, ipops = fit_interactions(rows, ar_mode)
            interaction_block_rows.append({
                "AR coding": ar_mode,
                "Trait": trait_label,
                "N": len(rows),
                "Parameters": 21,
                "R2 without interactions": float(base.rsquared),
                "R2 with interactions": float(full.rsquared),
                "Partial R2 interactions": ipartial,
                "Wald chi2": iwald,
                "df": 9,
                "Interaction-block p": ip,
                "Condition number": float(np.linalg.cond(full.model.exog)),
            })

            for j, label in enumerate(labels):
                hadza_beta = float(full.params[3 + j])
                interaction_beta = float(full.params[12 + j])
                datoga_beta = hadza_beta + interaction_beta

                contrast = np.zeros(len(full.params))
                contrast[3 + j] = 1.0
                contrast[12 + j] = 1.0
                datoga_test = full.t_test(contrast)
                datoga_ci = datoga_test.conf_int()[0]
                int_ci = full.conf_int()[12 + j]

                individual_interaction_rows.append({
                    "AR coding": ar_mode,
                    "Trait": trait_label,
                    "Marker": label,
                    "N": len(rows),
                    "Effect in Hadza, SD": hadza_beta,
                    "Hadza p": float(full.pvalues[3 + j]),
                    "Effect in Datoga, SD": datoga_beta,
                    "Datoga CI lower": float(datoga_ci[0]),
                    "Datoga CI upper": float(datoga_ci[1]),
                    "Datoga p": float(datoga_test.pvalue),
                    "Datoga-Hadza interaction, SD": interaction_beta,
                    "Interaction CI lower": float(int_ci[0]),
                    "Interaction CI upper": float(int_ci[1]),
                    "Interaction p": float(full.pvalues[12 + j]),
                })

    blocks = pd.DataFrame(block_rows)
    coefs = pd.DataFrame(coefficient_rows)
    iblocks = pd.DataFrame(interaction_block_rows)
    icoefs = pd.DataFrame(individual_interaction_rows)

    blocks["Block q BH (4)"] = np.nan
    blocks["Population q BH (4)"] = np.nan
    iblocks["Interaction-block q BH (4)"] = np.nan
    coefs["q BH (36)"] = np.nan
    icoefs["Interaction q BH (36)"] = np.nan

    for ar_mode in ("continuous", "category"):
        mask = blocks["AR coding"].eq(ar_mode)
        blocks.loc[mask, "Block q BH (4)"] = bh(blocks.loc[mask, "Block p"])
        blocks.loc[mask, "Population q BH (4)"] = bh(blocks.loc[mask, "Population p"])

        mask_i = iblocks["AR coding"].eq(ar_mode)
        iblocks.loc[mask_i, "Interaction-block q BH (4)"] = bh(
            iblocks.loc[mask_i, "Interaction-block p"]
        )

        mask_c = coefs["AR coding"].eq(ar_mode)
        coefs.loc[mask_c, "q BH (36)"] = bh(coefs.loc[mask_c, "p"])

        mask_ic = icoefs["AR coding"].eq(ar_mode)
        icoefs.loc[mask_ic, "Interaction q BH (36)"] = bh(
            icoefs.loc[mask_ic, "Interaction p"]
        )

    return blocks, coefs, iblocks, icoefs


def validate(blocks, coefs, iblocks, icoefs):
    for (ar_mode, trait), exp in EXPECTED.items():
        label = TRAITS[trait]
        row = blocks[(blocks["AR coding"] == ar_mode) & (blocks["Trait"] == label)].iloc[0]
        for key, column in (
            ("N", "N"),
            ("R2", "R2 full model"),
            ("partial", "Partial R2 genetic block"),
            ("block_p", "Block p"),
            ("population_beta", "Datoga-Hadza effect, points"),
            ("population_p", "Population p"),
        ):
            observed = float(row[column])
            expected = float(exp[key])
            if abs(observed - expected) > 5e-10:
                raise RuntimeError(
                    f"Validation failed for {ar_mode} {trait} {column}: "
                    f"observed={observed}, expected={expected}"
                )

    checks = {
        "DRD4 promoter: S-allele dosage": 0.3342070999160691,
        "AR: CAG repeat count": -0.796197670729479,
        "DRD4 exon III: L+ versus SS": 0.4851003746423457,
    }
    for marker, expected in checks.items():
        row = coefs[
            (coefs["AR coding"] == "continuous")
            & (coefs["Trait"] == "Physical aggression")
            & (coefs["Marker"] == marker)
        ].iloc[0]
        if abs(float(row["beta, points"]) - expected) > 5e-10:
            raise RuntimeError(f"Coefficient validation failed for {marker}")

    row = iblocks[
        (iblocks["AR coding"] == "continuous")
        & (iblocks["Trait"] == "Physical aggression")
    ].iloc[0]
    if abs(float(row["R2 without interactions"]) - 0.07622522939622756) > 5e-10:
        raise RuntimeError("Interaction base-model R2 validation failed")
    if abs(float(row["Interaction-block p"]) - 0.5676196146862817) > 5e-10:
        raise RuntimeError("Interaction block p validation failed")

    row = icoefs[
        (icoefs["AR coding"] == "continuous")
        & (icoefs["Trait"] == "Physical aggression")
        & (icoefs["Marker"] == "DRD4 promoter: S-allele dosage")
    ].iloc[0]
    if abs(float(row["Effect in Hadza, SD"]) - 0.09411129267220981) > 5e-10:
        raise RuntimeError("Hadza interaction coefficient validation failed")
    if abs(float(row["Datoga-Hadza interaction, SD"]) - (-0.07741218790784031)) > 5e-10:
        raise RuntimeError("Interaction coefficient validation failed")

    print("Validation PASSED: frozen S5a-S5c targets reproduced.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/Hadza_Datoga_2012_cleaned_males.xlsx"),
    )
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()

    df = pd.read_excel(args.input, sheet_name="Объединенные_мужчины")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    blocks, coefs, iblocks, icoefs = run(df)
    blocks.to_csv(args.output_dir / "S5a_HD_Genetic_Block.csv", index=False)
    coefs.to_csv(args.output_dir / "S5b_HD_Coefficients.csv", index=False)
    iblocks.to_csv(args.output_dir / "S5c_HD_Interaction_Block.csv", index=False)
    icoefs.to_csv(args.output_dir / "S5c_HD_Individual_Interactions.csv", index=False)

    validate(blocks, coefs, iblocks, icoefs)
    print(f"Outputs written to: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
