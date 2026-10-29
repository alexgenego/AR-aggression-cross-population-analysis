#!/usr/bin/env python3
"""
Reproduce the complete age models summarized in Supporting Information S8.

Restricted local inputs
-----------------------
data/Armenian_full_analysis.xlsx
  sheet: "первичные данные"

data/Hadza_Datoga_2012_cleaned_males.xlsx
  sheet: "Объединенные_мужчины"

All reported coefficients are standardized-effect coefficients. HC3 robust
covariance estimates are used.

Analysis families
-----------------
Armenians, both sexes:
- standardized outcome ~ standardized age + sex
- standardized outcome ~ standardized age + sex + age*sex
- standardized outcome ~ standardized age + sex + age^2

Hadza and Datoga men:
- simple within-population standardized age slope
- standardized age slope adjusted for the complete nine-variant block
- within-population quadratic age term

Hadza + Datoga:
- age-by-population interaction adjusted for the nine-variant block, using
  within-population standardized outcome, age, and continuous/dosage genetic
  predictors.

Output
------
outputs/S8_Age.csv
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


ARM_TRAITS = (
    ("Physical aggression", "PhA"),
    ("Verbal aggression", "VerbA"),
    ("Anger", "Ang"),
    ("Hostility", "Host"),
    ("Proactive aggression", "Proact_ag"),
    ("Reactive aggression", "React_ag"),
)

HD_TRAITS = (
    ("Physical aggression", "physag"),
    ("Verbal aggression", "verrbag"),
    ("Anger", "anger"),
    ("Hostility", "hostil"),
)

CONTINUOUS_GENETIC_COLUMNS = (0, 1, 2, 3, 5, 6, 7, 8)


def bh(values):
    return multipletests(np.asarray(values, dtype=float), method="fdr_bh")[1]


def zscore(values):
    values = np.asarray(values, dtype=float)
    return (values - values.mean()) / values.std(ddof=1)


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
    mapping = {"s": "S", "sa": "Sa", "lg": "Lg", "la": "La", "l": "L", "xl": "xL"}
    normalized = [mapping.get(part.strip().lower()) for part in parts]
    if any(x is None for x in normalized):
        return None
    if any(x in {"L", "xL"} for x in normalized):
        return None
    return float(sum(x in {"S", "Sa", "Lg"} for x in normalized))


def parse_ar(value) -> Optional[float]:
    if pd.isna(value):
        return None
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    return x if 10 <= x <= 40 else None


def repeat_dosage(value, target: int) -> Optional[float]:
    if pd.isna(value):
        return None
    numbers = [int(x) for x in re.findall(r"\d+", str(value))]
    if len(numbers) < 2:
        return None
    return float(sum(x == target for x in numbers[:2]))


def drd4_long(value) -> Optional[float]:
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


def encode_hd(row) -> Optional[np.ndarray]:
    values = [
        parse_drd4p(row["DRD4P"]),
        parse_sertpr(row["SERTPr"]),
        parse_ar(row["AR"]),
        repeat_dosage(row["SERT_In2"], 12),
        drd4_long(row["DRD4E3"]),
        parse_maoa(row["MAOA"]),
        repeat_dosage(row["DAT1_VNTR"], 10),
        parse_drd2(row["DRD2Taq1"]),
        parse_comt(row["COMT"]),
    ]
    if any(value is None for value in values):
        return None
    return np.asarray(values, dtype=float)


def add_row(rows, sample, sex, trait, term, n, ages, beta, ci, p, interpretation, family):
    rows.append({
        "Sample": sample,
        "Sex": sex,
        "Trait": trait,
        "Model/tested term": term,
        "N": int(n),
        "Min age": float(np.min(ages)),
        "Max age": float(np.max(ages)),
        "beta, SD": float(beta),
        "95% CI lower": float(ci[0]),
        "95% CI upper": float(ci[1]),
        "p": float(p),
        "Interpretation": interpretation,
        "_family": family,
    })


def armenian_models(df, rows):
    for trait_label, col in ARM_TRAITS:
        sub = df.copy()
        sub["age_value"] = pd.to_numeric(sub["AGE"], errors="coerce")
        sub["outcome"] = pd.to_numeric(sub[col], errors="coerce")
        if col == "Ang":
            sub.loc[sub["outcome"].eq(50), "outcome"] = np.nan
        sub = sub[
            sub["age_value"].notna()
            & sub["outcome"].notna()
            & sub["SEX"].notna()
        ].copy()

        age_z = zscore(sub["age_value"])
        y_z = zscore(sub["outcome"])
        sex = sub["SEX"].eq(2).astype(float).to_numpy()

        # Linear age + sex.
        X = np.column_stack([np.ones(len(sub)), age_z, sex])
        model = sm.OLS(y_z, X).fit(cov_type="HC3")
        add_row(
            rows, "Armenians", "Both sexes", trait_label,
            "Linear age + sex", len(sub), sub["age_value"],
            model.params[1], model.conf_int()[1], model.pvalues[1],
            "Age and outcome standardized", "arm_linear",
        )

        # Age x sex.
        X = np.column_stack([
            np.ones(len(sub)), age_z, sex, age_z * sex
        ])
        model = sm.OLS(y_z, X).fit(cov_type="HC3")
        add_row(
            rows, "Armenians", "Both sexes", trait_label,
            "Age x sex", len(sub), sub["age_value"],
            model.params[3], model.conf_int()[3], model.pvalues[3],
            "Difference in age slope for women relative to men", "arm_interaction",
        )

        # Quadratic term.
        X = np.column_stack([
            np.ones(len(sub)), age_z, sex, age_z ** 2
        ])
        model = sm.OLS(y_z, X).fit(cov_type="HC3")
        add_row(
            rows, "Armenians", "Both sexes", trait_label,
            "Nonlinearity: age^2", len(sub), sub["age_value"],
            model.params[3], model.conf_int()[3], model.pvalues[3],
            "Test of a U-shaped or inverted U-shaped association", "arm_quadratic",
        )


def hd_simple_and_adjusted(df, rows):
    for population, source_name in (("Hadza", "Хадза"), ("Datoga", "Датога")):
        for trait_label, col in HD_TRAITS:
            sub = df[df["population"].eq(source_name)].copy()
            sub["age_value"] = pd.to_numeric(sub["age"], errors="coerce")
            sub["outcome"] = pd.to_numeric(sub[col], errors="coerce")
            valid = sub["age_value"].notna() & sub["outcome"].notna()
            simple = sub.loc[valid].copy()

            age_z = zscore(simple["age_value"])
            y_z = zscore(simple["outcome"])

            model = sm.OLS(
                y_z, np.column_stack([np.ones(len(simple)), age_z])
            ).fit(cov_type="HC3")
            add_row(
                rows, population, "Men", trait_label,
                "Linear age", len(simple), simple["age_value"],
                model.params[1], model.conf_int()[1], model.pvalues[1],
                "Within-population model", f"{population}_linear",
            )

            Xq = np.column_stack([
                np.ones(len(simple)), age_z, age_z ** 2
            ])
            qmodel = sm.OLS(y_z, Xq).fit(cov_type="HC3")
            add_row(
                rows, population, "Men", trait_label,
                "Nonlinearity: age^2", len(simple), simple["age_value"],
                qmodel.params[2], qmodel.conf_int()[2], qmodel.pvalues[2],
                "Test of a curvilinear association", f"{population}_quadratic",
            )

            complete_records = []
            for _, row in simple.iterrows():
                g = encode_hd(row)
                if g is not None:
                    complete_records.append(
                        (float(row["age_value"]), float(row["outcome"]), g)
                    )

            ages = np.asarray([r[0] for r in complete_records])
            outcomes = np.asarray([r[1] for r in complete_records])
            G = np.vstack([r[2] for r in complete_records])
            age_z_c = zscore(ages)
            y_z_c = zscore(outcomes)

            X = np.column_stack([np.ones(len(ages)), age_z_c, G])
            adjusted = sm.OLS(y_z_c, X).fit(cov_type="HC3")
            add_row(
                rows, population, "Men", trait_label,
                "Age + 9 genetic variables", len(ages), ages,
                adjusted.params[1], adjusted.conf_int()[1], adjusted.pvalues[1],
                "Complete-case; AR as continuous CAG repeat number",
                f"{population}_genetic",
            )


def hd_interactions(df, rows):
    for trait_label, col in HD_TRAITS:
        records = []
        for _, row in df.iterrows():
            if pd.isna(row["age"]) or pd.isna(row[col]):
                continue
            g = encode_hd(row)
            if g is None:
                continue
            population = "Hadza" if row["population"] == "Хадза" else "Datoga"
            records.append(
                (population, float(row["age"]), float(row[col]), g)
            )

        pops = np.asarray([r[0] for r in records])
        ages = np.asarray([r[1] for r in records])
        outcomes = np.asarray([r[2] for r in records])
        G = np.vstack([r[3] for r in records])

        age_z = ages.copy()
        y_z = outcomes.copy()
        G_z = G.copy()

        for population in ("Hadza", "Datoga"):
            mask = pops == population
            age_z[mask] = zscore(ages[mask])
            y_z[mask] = zscore(outcomes[mask])
            for j in CONTINUOUS_GENETIC_COLUMNS:
                G_z[mask, j] = zscore(G[mask, j])

        pop_dummy = (pops == "Datoga").astype(float)
        interaction = age_z * pop_dummy
        X = np.column_stack([
            np.ones(len(ages)), age_z, pop_dummy, G_z, interaction
        ])
        model = sm.OLS(y_z, X).fit(cov_type="HC3")

        add_row(
            rows, "Hadza + Datoga", "Men", trait_label,
            "Age x population + 9 markers", len(ages), ages,
            model.params[-1], model.conf_int()[-1], model.pvalues[-1],
            "Difference in age slopes, Datoga-Hadza", "HD_interaction",
        )


def run(armenian, hd):
    rows: List[Dict] = []
    armenian_models(armenian, rows)
    hd_simple_and_adjusted(hd, rows)
    hd_interactions(hd, rows)

    out = pd.DataFrame(rows)
    out["within-family q"] = np.nan
    for family in out["_family"].unique():
        mask = out["_family"].eq(family)
        out.loc[mask, "within-family q"] = bh(out.loc[mask, "p"])

    return out.drop(columns=["_family"])


def validate(out):
    checks = [
        (
            ("Armenians", "Physical aggression", "Linear age + sex"),
            -0.04583417433615363, 0.339327723013471, 185,
        ),
        (
            ("Hadza", "Anger", "Linear age"),
            -0.21186854918160658, 0.021876931553997305, 116,
        ),
        (
            ("Datoga", "Physical aggression", "Linear age"),
            0.23793621679793786, 0.005340177869266138, 85,
        ),
        (
            ("Hadza + Datoga", "Anger", "Age x population + 9 markers"),
            0.2779121153435288, 0.08526626800685032, 154,
        ),
    ]
    for (sample, trait, term), beta, p, n in checks:
        row = out[
            (out["Sample"] == sample)
            & (out["Trait"] == trait)
            & (out["Model/tested term"] == term)
        ].iloc[0]
        if int(row["N"]) != n:
            raise RuntimeError(f"S8 N validation failed for {sample} {trait}")
        if abs(float(row["beta, SD"]) - beta) > 5e-10:
            raise RuntimeError(f"S8 beta validation failed for {sample} {trait}")
        if abs(float(row["p"]) - p) > 5e-10:
            raise RuntimeError(f"S8 p validation failed for {sample} {trait}")

    print("Validation PASSED: frozen S8 targets reproduced.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--armenian",
        type=Path,
        default=Path("data/Armenian_full_analysis.xlsx"),
    )
    parser.add_argument(
        "--hadza-datoga",
        type=Path,
        default=Path("data/Hadza_Datoga_2012_cleaned_males.xlsx"),
    )
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()

    armenian = pd.read_excel(args.armenian, sheet_name="первичные данные")
    hd = pd.read_excel(
        args.hadza_datoga, sheet_name="Объединенные_мужчины"
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)

    out = run(armenian, hd)
    out.to_csv(args.output_dir / "S8_Age.csv", index=False)
    validate(out)
    print(f"Outputs written to: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
