#!/usr/bin/env python3
"""
Reproduce the power calculations summarized in Supporting Information S9.

The historical calculation can be recovered exactly from the frozen results:

1. Fit the age-adjusted standardized single-variant model with HC3 robust SE.
2. Treat |beta / robust_SE| as the estimated noncentrality parameter.
3. Current post-hoc design power is evaluated from a noncentral t distribution
   using df = N - 3 (intercept, age, genotype).
4. For a future total sample size N_new, noncentrality scales as
      delta_new = delta_current * sqrt(N_new / N_current)
5. Conservative planning reduces the estimated signal by 25%.
6. Familywise planning scenarios use alpha=0.0125 (four traits in one locus)
   and alpha=0.00625 (eight tests for two loci).

Current power is descriptive only and is not independent evidence for an
association.

Restricted local inputs
-----------------------
data/Armenian_full_analysis.xlsx
  sheet: "первичные данные"

data/Hadza_Datoga_2012_cleaned_males.xlsx
  sheet: "Объединенные_мужчины"

Outputs
-------
outputs/S9_Power_within.csv
outputs/S9_Power_between_population.csv
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import nct, norm, t
from statsmodels.stats.multitest import multipletests


TRAITS = (
    ("Physical aggression", "physag", "PhA"),
    ("Verbal aggression", "verrbag", "VerbA"),
    ("Anger", "anger", "Ang"),
    ("Hostility", "hostil", "Host"),
)
POPULATIONS = ("Hadza", "Datoga", "Armenians")


def bh(values):
    return multipletests(np.asarray(values, dtype=float), method="fdr_bh")[1]


def parse_drd4(value) -> Optional[float]:
    if pd.isna(value):
        return None
    numbers = [int(x) for x in re.findall(r"\d+", str(value))]
    if len(numbers) < 2:
        return None
    return float(max(numbers[:2]) >= 7)


def ar_binary(value) -> Optional[float]:
    if pd.isna(value):
        return None
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    if not 10 <= x <= 40:
        return None
    return float(x >= 23)


def zscore(values):
    values = np.asarray(values, dtype=float)
    return (values - values.mean()) / values.std(ddof=1)


def population_dataset(
    population: str,
    trait_label: str,
    marker: str,
    armenian: pd.DataFrame,
    hd: pd.DataFrame,
):
    trait_info = {label: (hd_col, arm_col) for label, hd_col, arm_col in TRAITS}
    hd_col, arm_col = trait_info[trait_label]

    records = []

    if population in ("Hadza", "Datoga"):
        source_name = "Хадза" if population == "Hadza" else "Датога"
        sub = hd[hd["population"].eq(source_name)].copy()
        for _, row in sub.iterrows():
            if pd.isna(row["age"]) or pd.isna(row[hd_col]):
                continue
            genotype = (
                ar_binary(row["AR"])
                if marker == "AR"
                else parse_drd4(row["DRD4E3"])
            )
            if genotype is None:
                continue
            records.append(
                (float(row["age"]), float(row[hd_col]), genotype)
            )
    else:
        sub = armenian[armenian["SEX"].eq(1)].copy()
        for _, row in sub.iterrows():
            outcome = row[arm_col]
            if arm_col == "Ang" and outcome == 50:
                outcome = np.nan
            if pd.isna(row["AGE"]) or pd.isna(outcome):
                continue
            genotype = (
                ar_binary(row["AR_CAG"])
                if marker == "AR"
                else parse_drd4(row["DRD4eIII"])
            )
            if genotype is None:
                continue
            records.append(
                (float(row["AGE"]), float(outcome), genotype)
            )

    age = np.asarray([r[0] for r in records])
    outcome = np.asarray([r[1] for r in records])
    genotype = np.asarray([r[2] for r in records])
    return age, outcome, genotype


def fit_simple_effect(
    population: str,
    trait_label: str,
    marker: str,
    armenian: pd.DataFrame,
    hd: pd.DataFrame,
):
    age, outcome, genotype = population_dataset(
        population, trait_label, marker, armenian, hd
    )
    age_z = zscore(age)
    outcome_z = zscore(outcome)

    X = np.column_stack([np.ones(len(outcome)), genotype, age_z])
    model = sm.OLS(outcome_z, X).fit(cov_type="HC3")

    return {
        "N": len(outcome),
        "group0": int(np.sum(genotype == 0)),
        "group1": int(np.sum(genotype == 1)),
        "beta": float(model.params[1]),
        "SE": float(model.bse[1]),
        "p": float(model.pvalues[1]),
        "CI": np.asarray(model.conf_int())[1],
    }


def current_power(delta: float, n: int, alpha: float = 0.05) -> float:
    df = n - 3
    critical = t.ppf(1 - alpha / 2, df)
    return float(
        nct.cdf(-critical, df, delta)
        + 1
        - nct.cdf(critical, df, delta)
    )


def power_at_future_n(
    delta_current: float,
    n_current: int,
    n_new: int,
    alpha: float,
    signal_scale: float = 1.0,
) -> float:
    delta_new = (
        abs(delta_current)
        * signal_scale
        * np.sqrt(n_new / n_current)
    )
    return current_power(delta_new, n_new, alpha)


def required_n(
    delta_current: float,
    n_current: int,
    target_power: float,
    alpha: float,
    signal_scale: float = 1.0,
    max_n: int = 100000,
) -> int:
    for n in range(4, max_n + 1):
        if power_at_future_n(
            delta_current,
            n_current,
            n,
            alpha,
            signal_scale,
        ) >= target_power:
            return n
    raise RuntimeError("Required N exceeds search bound")


def all_marker_effects(marker, armenian, hd):
    rows = []
    for population in POPULATIONS:
        for trait_label, _, _ in TRAITS:
            effect = fit_simple_effect(
                population, trait_label, marker, armenian, hd
            )
            rows.append({
                "Marker": marker,
                "Population": population,
                "Trait": trait_label,
                **effect,
            })
    out = pd.DataFrame(rows)
    out["Gene-specific BH q (12)"] = bh(out["p"])
    return out


def planning_row(effect, alpha):
    delta = abs(effect["beta"] / effect["SE"])
    n = int(effect["N"])
    return {
        "N80": required_n(delta, n, 0.80, alpha, 1.0),
        "N90": required_n(delta, n, 0.90, alpha, 1.0),
    }


def run_within(armenian, hd):
    ar = all_marker_effects("AR", armenian, hd)
    dr = all_marker_effects("DRD4", armenian, hd)
    all_effects = pd.concat([ar, dr], ignore_index=True)

    selections = [
        ("DRD4", "Armenians", "Physical aggression"),
        ("AR", "Datoga", "Anger"),
        ("AR", "Datoga", "Physical aggression"),
        ("DRD4", "Armenians", "Verbal aggression"),
    ]

    rows: List[Dict] = []
    for marker, population, trait in selections:
        row = all_effects[
            (all_effects["Marker"] == marker)
            & (all_effects["Population"] == population)
            & (all_effects["Trait"] == trait)
        ].iloc[0]

        delta = abs(float(row["beta"]) / float(row["SE"]))
        n = int(row["N"])
        p05 = {
            target: required_n(delta, n, target, 0.05)
            for target in (0.80, 0.90)
        }
        p0125 = {
            target: required_n(delta, n, target, 0.0125)
            for target in (0.80, 0.90)
        }
        p00625 = {
            target: required_n(delta, n, target, 0.00625)
            for target in (0.80, 0.90)
        }
        conservative_0125 = {
            target: required_n(delta, n, target, 0.0125, 0.75)
            for target in (0.80, 0.90)
        }
        conservative_00625 = {
            target: required_n(delta, n, target, 0.00625, 0.75)
            for target in (0.80, 0.90)
        }

        n90_conservative = conservative_0125[0.90]
        group0_future = int(round(
            n90_conservative * int(row["group0"]) / n
        ))
        group1_future = n90_conservative - group0_future

        rows.append({
            "Marker": "AR >=23 vs <=22" if marker == "AR" else "DRD4E3 L+ vs SS",
            "Population": population,
            "Trait": trait,
            "Current N": n,
            "Group 0": int(row["group0"]),
            "Group 1": int(row["group1"]),
            "Effect, SD": float(row["beta"]),
            "SE, SD": float(row["SE"]),
            "p": float(row["p"]),
            "Gene-specific BH q": float(row["Gene-specific BH q (12)"]),
            "Current power, alpha=0.05": current_power(delta, n, 0.05),
            "N80, alpha=0.05": p05[0.80],
            "N90, alpha=0.05": p05[0.90],
            "N80, alpha=0.0125": p0125[0.80],
            "N90, alpha=0.0125": p0125[0.90],
            "Conservative N80, alpha=0.0125": conservative_0125[0.80],
            "Conservative N90, alpha=0.0125": n90_conservative,
            "Group 0 at conservative N90": group0_future,
            "Group 1 at conservative N90": group1_future,
            "N80, alpha=0.00625": p00625[0.80],
            "N90, alpha=0.00625": p00625[0.90],
            "Conservative N80, alpha=0.00625": conservative_00625[0.80],
            "Conservative N90, alpha=0.00625": conservative_00625[0.90],
        })

    return pd.DataFrame(rows)


def pair_interaction(
    population_a,
    population_b,
    trait,
    marker,
    armenian,
    hd,
):
    records = []
    separate_effects = {}

    for population in (population_a, population_b):
        age, outcome, genotype = population_dataset(
            population, trait, marker, armenian, hd
        )
        age_z = zscore(age)
        outcome_z = zscore(outcome)

        # Separate population effect for table display.
        sep = sm.OLS(
            outcome_z,
            np.column_stack([np.ones(len(outcome)), genotype, age_z]),
        ).fit(cov_type="HC3")
        separate_effects[population] = float(sep.params[1])

        for y, a, g in zip(outcome_z, age_z, genotype):
            records.append((population, y, a, g))

    frame = pd.DataFrame(
        records, columns=["population", "outcome", "age", "genotype"]
    )
    pop = frame["population"].eq(population_b).astype(float).to_numpy()
    genotype = frame["genotype"].to_numpy()
    X = np.column_stack([
        np.ones(len(frame)),
        frame["age"].to_numpy(),
        pop,
        genotype,
        genotype * pop,
    ])
    model = sm.OLS(frame["outcome"].to_numpy(), X).fit(cov_type="HC3")

    return {
        "N": len(frame),
        "effect_A": separate_effects[population_a],
        "effect_B": separate_effects[population_b],
        "difference": float(model.params[4]),
        "SE": float(model.bse[4]),
        "CI": np.asarray(model.conf_int())[4],
        "p": float(model.pvalues[4]),
    }


def run_between(armenian, hd):
    tests = [
        ("AR: hostility, Armenians - Datoga", "Datoga", "Armenians", "Hostility", "AR"),
        ("DRD4: physical aggression, Armenians - Datoga", "Datoga", "Armenians", "Physical aggression", "DRD4"),
        ("DRD4: verbal aggression, Armenians - Datoga", "Datoga", "Armenians", "Verbal aggression", "DRD4"),
        ("DRD4: verbal aggression, Armenians - Hadza", "Hadza", "Armenians", "Verbal aggression", "DRD4"),
    ]

    rows: List[Dict] = []
    for label, pop_a, pop_b, trait, marker in tests:
        effect = pair_interaction(
            pop_a, pop_b, trait, marker, armenian, hd
        )
        delta = abs(effect["difference"] / effect["SE"])
        n = effect["N"]

        rows.append({
            "Test": label,
            "N": n,
            "Effect A, SD": effect["effect_A"],
            "Effect B, SD": effect["effect_B"],
            "Difference, SD": effect["difference"],
            "95% CI lower": effect["CI"][0],
            "95% CI upper": effect["CI"][1],
            "p": effect["p"],
            "Current power": current_power(delta, n, 0.05),
            "N80": required_n(delta, n, 0.80, 0.05),
            "N90": required_n(delta, n, 0.90, 0.05),
            "N90 with 25% effect-size reduction": required_n(
                delta, n, 0.90, 0.05, 0.75
            ),
        })

    return pd.DataFrame(rows)


def validate(within, between):
    first = within.iloc[0]
    checks = {
        "Current power, alpha=0.05": 0.6686928822594669,
        "N80, alpha=0.05": 93,
        "N90, alpha=0.05": 123,
        "N80, alpha=0.0125": 132,
        "N90, alpha=0.0125": 168,
        "Conservative N80, alpha=0.0125": 232,
        "Conservative N90, alpha=0.0125": 295,
        "N80, alpha=0.00625": 151,
        "N90, alpha=0.00625": 190,
        "Conservative N80, alpha=0.00625": 266,
        "Conservative N90, alpha=0.00625": 334,
    }
    for column, expected in checks.items():
        observed = float(first[column])
        if abs(observed - float(expected)) > 5e-10:
            raise RuntimeError(f"S9 validation failed for {column}")

    row = between[
        between["Test"].eq(
            "DRD4: physical aggression, Armenians - Datoga"
        )
    ].iloc[0]
    for column, expected in (
        ("Current power", 0.61038258194462),
        ("N80", 232),
        ("N90", 310),
        ("N90 with 25% effect-size reduction", 550),
    ):
        if abs(float(row[column]) - float(expected)) > 5e-10:
            raise RuntimeError(f"S9 between validation failed for {column}")

    print("Validation PASSED: frozen S9 power targets reproduced.")


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

    within = run_within(armenian, hd)
    between = run_between(armenian, hd)

    within.to_csv(args.output_dir / "S9_Power_within.csv", index=False)
    between.to_csv(
        args.output_dir / "S9_Power_between_population.csv",
        index=False,
    )

    validate(within, between)
    print(f"Outputs written to: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
