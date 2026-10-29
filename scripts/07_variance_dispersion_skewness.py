#!/usr/bin/env python3
"""
Reconstruct the cross-population variance, dispersion, and skewness analyses
summarized in Supporting Information S7a-S7d.

Required restricted local inputs
--------------------------------
data/Armenian_full_analysis.xlsx
  sheet: "первичные данные"

data/Hadza_Datoga_2012_cleaned_males.xlsx
  sheet: "Объединенные_мужчины"

Design
------
For each population and aggression trait, a simple linear age trajectory is
fitted using all men with valid age and that trait, irrespective of genotype.
Residuals remain in the original score units. AR and DRD4 categories are then
applied to those age-adjusted residuals.

Within-population variance:
  Brown-Forsythe (median-centered Levene), Fligner-Killeen, SD ratio, and
  normal-consistency-scaled MAD ratio.

Between-population variance:
  omnibus and pairwise Brown-Forsythe/Fligner-Killeen tests within the same
  genotype category.

Skewness:
  Fisher-Pearson unbiased skewness, scipy skewtest, Bowley quartile skewness,
  and category location summaries.

BOOTSTRAP NOTE
--------------
The historical bootstrap random state was not preserved. Bootstrap intervals
and bootstrap skewness-difference p values are therefore transparent seeded
reconstructions (default seed 2026; 2000 resamples). Deterministic quantities
such as N, age residuals, SDs, Brown-Forsythe tests, Fisher skewness, and
skewtest p values reproduce the frozen analysis. Small Fligner-Killeen
differences can occur across SciPy versions because of implementation/tie
handling.

Outputs
-------
outputs/S7a_WithinPop_Variance_reconstructed.csv
outputs/S7b_BetweenPop_Variance_reconstructed.csv
outputs/S7c_Skewness_reconstructed.csv
outputs/S7d_Skewness_Differences_reconstructed.csv
"""

from __future__ import annotations

import argparse
import itertools
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests


TRAITS = (
    ("Physical aggression", "physag", "PhA"),
    ("Verbal aggression", "verrbag", "VerbA"),
    ("Anger", "anger", "Ang"),
    ("Hostility", "hostil", "Host"),
)

POPULATIONS = ("Hadza", "Datoga", "Armenians")
MARKERS = ("AR", "DRD4")

EXPECTED = {
    ("within", "Datoga", "DRD4", "Verbal aggression"): {
        "N0": 51,
        "N1": 31,
        "SD0": 3.3361352313891723,
        "SD1": 4.758264948261124,
        "ratio": 1.4262805966291046,
        "bf_p": 0.01963083274014686,
    },
    ("within", "Datoga", "AR", "Physical aggression"): {
        "N0": 32,
        "N1": 40,
        "ratio": 0.74485388284043,
        "bf_p": 0.07714452448773353,
    },
    ("within", "Datoga", "AR", "Anger"): {
        "N0": 29,
        "N1": 38,
        "ratio": 0.7226537927477451,
        "bf_p": 0.06174906786991511,
    },
    ("between", "Armenians-Datoga", "DRD4", "Physical aggression", "L+"): {
        "ratio": 0.5856045133393328,
        "bf_p": 0.009748187761084473,
    },
    ("between", "Armenians-Hadza", "DRD4", "Physical aggression", "L+"): {
        "ratio": 0.5404880259550117,
        "bf_p": 0.00828228776098191,
    },
    ("skew", "Armenians", "AR", "Physical aggression", "<=22"): {
        "N": 36,
        "skew": 1.1900401808330576,
        "p": 0.005286179662302692,
    },
    ("skew", "Armenians", "AR", "Physical aggression", ">=23"): {
        "N": 30,
        "skew": -1.3518225452899286,
        "p": 0.004029717824370728,
    },
}


def bh(values: Iterable[float]) -> np.ndarray:
    return multipletests(np.asarray(list(values), dtype=float), method="fdr_bh")[1]


def parse_drd4(value) -> Optional[str]:
    if pd.isna(value):
        return None
    numbers = [int(x) for x in re.findall(r"\d+", str(value))]
    if len(numbers) < 2:
        return None
    return "L+" if max(numbers[:2]) >= 7 else "SS"


def ar_category(value) -> Optional[str]:
    if pd.isna(value):
        return None
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    if not 10 <= x <= 40:
        return None
    return ">=23" if x >= 23 else "<=22"


def build_population_tables(armenian: pd.DataFrame, hd: pd.DataFrame):
    result: Dict[Tuple[str, str], pd.DataFrame] = {}

    for trait_label, hd_col, arm_col in TRAITS:
        for population in ("Hadza", "Datoga"):
            source_name = "Хадза" if population == "Hadza" else "Датога"
            sub = hd[hd["population"].eq(source_name)].copy()
            sub["age_value"] = pd.to_numeric(sub["age"], errors="coerce")
            sub["outcome"] = pd.to_numeric(sub[hd_col], errors="coerce")

            valid = sub["age_value"].notna() & sub["outcome"].notna()
            model = sm.OLS(
                sub.loc[valid, "outcome"],
                sm.add_constant(sub.loc[valid, "age_value"]),
            ).fit()

            sub["residual"] = np.nan
            sub.loc[valid, "residual"] = (
                sub.loc[valid, "outcome"]
                - model.predict(sm.add_constant(sub.loc[valid, "age_value"]))
            )
            sub["AR_category"] = sub["AR"].apply(ar_category)
            sub["DRD4_category"] = sub["DRD4E3"].apply(parse_drd4)
            result[(population, trait_label)] = sub

        sub = armenian[armenian["SEX"].eq(1)].copy()
        sub["age_value"] = pd.to_numeric(sub["AGE"], errors="coerce")
        sub["outcome"] = pd.to_numeric(sub[arm_col], errors="coerce")
        if arm_col == "Ang":
            sub.loc[sub["outcome"].eq(50), "outcome"] = np.nan

        valid = sub["age_value"].notna() & sub["outcome"].notna()
        model = sm.OLS(
            sub.loc[valid, "outcome"],
            sm.add_constant(sub.loc[valid, "age_value"]),
        ).fit()

        sub["residual"] = np.nan
        sub.loc[valid, "residual"] = (
            sub.loc[valid, "outcome"]
            - model.predict(sm.add_constant(sub.loc[valid, "age_value"]))
        )
        sub["AR_category"] = sub["AR_CAG"].apply(ar_category)
        sub["DRD4_category"] = sub["DRD4eIII"].apply(parse_drd4)
        result[("Armenians", trait_label)] = sub

    return result


def groups_for(table: pd.DataFrame, marker: str):
    column = "AR_category" if marker == "AR" else "DRD4_category"
    labels = ("<=22", ">=23") if marker == "AR" else ("SS", "L+")
    groups = {}
    for label in labels:
        groups[label] = table.loc[
            table[column].eq(label), "residual"
        ].dropna().to_numpy(dtype=float)
    return labels, groups


def mad_normal(x: np.ndarray) -> float:
    return float(stats.median_abs_deviation(x, scale="normal"))


def bootstrap_ratio(
    group0: np.ndarray,
    group1: np.ndarray,
    statistic: str,
    rng: np.random.Generator,
    n_boot: int,
):
    values = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        a = rng.choice(group0, size=len(group0), replace=True)
        b = rng.choice(group1, size=len(group1), replace=True)
        if statistic == "sd":
            values[i] = np.std(b, ddof=1) / np.std(a, ddof=1)
        elif statistic == "mad":
            denom = mad_normal(a)
            values[i] = np.nan if denom == 0 else mad_normal(b) / denom
        else:
            raise KeyError(statistic)
    values = values[np.isfinite(values)]
    return np.quantile(values, [0.025, 0.975])


def run_within(tables, seed: int, n_boot: int):
    rows: List[Dict] = []
    seed_counter = 0

    for marker in MARKERS:
        for population in POPULATIONS:
            for trait_label, _, _ in TRAITS:
                labels, groups = groups_for(
                    tables[(population, trait_label)], marker
                )
                label0, label1 = labels
                x0, x1 = groups[label0], groups[label1]

                sd0, sd1 = np.std(x0, ddof=1), np.std(x1, ddof=1)
                mad0, mad1 = mad_normal(x0), mad_normal(x1)
                bf = stats.levene(x0, x1, center="median")
                fl = stats.fligner(x0, x1, center="median")

                rng_sd = np.random.default_rng(seed + 10000 + seed_counter)
                rng_mad = np.random.default_rng(seed + 20000 + seed_counter)
                seed_counter += 1
                sd_ci = bootstrap_ratio(x0, x1, "sd", rng_sd, n_boot)
                mad_ci = bootstrap_ratio(x0, x1, "mad", rng_mad, n_boot)

                rows.append({
                    "Marker": marker,
                    "Population": population,
                    "Trait": trait_label,
                    "Category 0": label0,
                    "Category 1": label1,
                    "N0": len(x0),
                    "N1": len(x1),
                    "SD0": sd0,
                    "SD1": sd1,
                    "SD1/SD0": sd1 / sd0,
                    "Bootstrap SD ratio CI low": sd_ci[0],
                    "Bootstrap SD ratio CI high": sd_ci[1],
                    "MAD0": mad0,
                    "MAD1": mad1,
                    "MAD1/MAD0": mad1 / mad0 if mad0 != 0 else np.nan,
                    "Bootstrap MAD ratio CI low": mad_ci[0],
                    "Bootstrap MAD ratio CI high": mad_ci[1],
                    "Brown-Forsythe p": bf.pvalue,
                    "Fligner-Killeen p": fl.pvalue,
                    "Bootstrap seed": seed,
                    "Bootstrap resamples": n_boot,
                })

    out = pd.DataFrame(rows)
    out["Brown-Forsythe q (24)"] = bh(out["Brown-Forsythe p"])
    out["Fligner-Killeen q (24)"] = bh(out["Fligner-Killeen p"])
    return out


def run_between(tables):
    rows: List[Dict] = []

    for marker in MARKERS:
        categories = ("<=22", ">=23") if marker == "AR" else ("SS", "L+")
        for trait_label, _, _ in TRAITS:
            for category in categories:
                samples = {}
                for population in POPULATIONS:
                    labels, groups = groups_for(
                        tables[(population, trait_label)], marker
                    )
                    samples[population] = groups[category]

                bf = stats.levene(
                    *(samples[p] for p in POPULATIONS),
                    center="median",
                )
                fl = stats.fligner(
                    *(samples[p] for p in POPULATIONS),
                    center="median",
                )
                rows.append({
                    "Marker": marker,
                    "Trait": trait_label,
                    "Category": category,
                    "Test type": "omnibus",
                    "Comparison": "3 populations",
                    "SD B / SD A": np.nan,
                    "Brown-Forsythe p": bf.pvalue,
                    "Fligner-Killeen p": fl.pvalue,
                })

                for pop_a, pop_b in (
                    ("Hadza", "Datoga"),
                    ("Hadza", "Armenians"),
                    ("Datoga", "Armenians"),
                ):
                    a, b = samples[pop_a], samples[pop_b]
                    bf2 = stats.levene(a, b, center="median")
                    fl2 = stats.fligner(a, b, center="median")
                    rows.append({
                        "Marker": marker,
                        "Trait": trait_label,
                        "Category": category,
                        "Test type": "pair",
                        "Comparison": f"{pop_b}-{pop_a}",
                        "SD B / SD A": np.std(b, ddof=1) / np.std(a, ddof=1),
                        "Brown-Forsythe p": bf2.pvalue,
                        "Fligner-Killeen p": fl2.pvalue,
                    })

    out = pd.DataFrame(rows)
    out["Brown-Forsythe q overall"] = bh(out["Brown-Forsythe p"])
    out["Fligner-Killeen q overall"] = bh(out["Fligner-Killeen p"])
    return out


def bowley_skew(x: np.ndarray) -> float:
    q1, med, q3 = np.quantile(x, [0.25, 0.50, 0.75])
    denom = q3 - q1
    return np.nan if denom == 0 else float((q3 + q1 - 2 * med) / denom)


def bootstrap_skew_ci(
    x: np.ndarray,
    rng: np.random.Generator,
    n_boot: int,
):
    values = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        sample = rng.choice(x, size=len(x), replace=True)
        values[i] = stats.skew(sample, bias=False)
    values = values[np.isfinite(values)]
    return np.quantile(values, [0.025, 0.975])


def run_skewness(tables, seed: int, n_boot: int):
    rows: List[Dict] = []
    counter = 0

    for marker in MARKERS:
        for population in POPULATIONS:
            for trait_label, _, _ in TRAITS:
                table = tables[(population, trait_label)]
                population_mean = float(table["outcome"].dropna().mean())
                labels, groups = groups_for(table, marker)
                category_column = (
                    "AR_category" if marker == "AR" else "DRD4_category"
                )

                for category in labels:
                    mask = table[category_column].eq(category) & table["residual"].notna()
                    residual = table.loc[mask, "residual"].to_numpy(dtype=float)
                    raw = table.loc[mask, "outcome"].to_numpy(dtype=float)

                    rng = np.random.default_rng(seed + 30000 + counter)
                    counter += 1
                    skew_ci = bootstrap_skew_ci(residual, rng, n_boot)

                    fisher = float(stats.skew(residual, bias=False))
                    skew_p = float(stats.skewtest(residual).pvalue)

                    rows.append({
                        "Marker": marker,
                        "Population": population,
                        "Trait": trait_label,
                        "Category": category,
                        "N": len(residual),
                        "Population raw mean": population_mean,
                        "Category raw mean": np.mean(raw),
                        "Category raw median": np.median(raw),
                        "Raw mean deviation": np.mean(raw) - population_mean,
                        "Age-adjusted mean deviation": np.mean(residual),
                        "Age-adjusted median deviation": np.median(residual),
                        "Fisher skewness": fisher,
                        "Bootstrap skew CI low": skew_ci[0],
                        "Bootstrap skew CI high": skew_ci[1],
                        "skewtest p": skew_p,
                        "Bowley skewness": bowley_skew(residual),
                        "Q1": np.quantile(raw, 0.25),
                        "Q3": np.quantile(raw, 0.75),
                        "P10": np.quantile(raw, 0.10),
                        "P90": np.quantile(raw, 0.90),
                        "Minimum": np.min(raw),
                        "Maximum": np.max(raw),
                        "Bootstrap seed": seed,
                        "Bootstrap resamples": n_boot,
                    })

    out = pd.DataFrame(rows)
    out["skewtest q (48)"] = bh(out["skewtest p"])
    return out


def run_skew_differences(tables, seed: int, n_boot: int):
    rows: List[Dict] = []
    counter = 0

    for marker in MARKERS:
        for population in POPULATIONS:
            for trait_label, _, _ in TRAITS:
                labels, groups = groups_for(
                    tables[(population, trait_label)], marker
                )
                label0, label1 = labels
                x0, x1 = groups[label0], groups[label1]
                observed = (
                    stats.skew(x1, bias=False)
                    - stats.skew(x0, bias=False)
                )

                rng = np.random.default_rng(seed + 40000 + counter)
                counter += 1
                bootstrap = np.empty(n_boot, dtype=float)
                for i in range(n_boot):
                    a = rng.choice(x0, size=len(x0), replace=True)
                    b = rng.choice(x1, size=len(x1), replace=True)
                    bootstrap[i] = (
                        stats.skew(b, bias=False)
                        - stats.skew(a, bias=False)
                    )
                bootstrap = bootstrap[np.isfinite(bootstrap)]
                ci = np.quantile(bootstrap, [0.025, 0.975])

                p_sign = 2 * min(
                    np.mean(bootstrap <= 0),
                    np.mean(bootstrap >= 0),
                )
                p_sign = min(float(p_sign), 1.0)

                rows.append({
                    "Marker": marker,
                    "Population": population,
                    "Trait": trait_label,
                    "Category 0": label0,
                    "Category 1": label1,
                    "Observed skewness difference 1-0": observed,
                    "Bootstrap mean skewness difference": np.mean(bootstrap),
                    "Bootstrap CI low": ci[0],
                    "Bootstrap CI high": ci[1],
                    "Reconstructed two-sided sign p": p_sign,
                    "Bootstrap seed": seed,
                    "Bootstrap resamples": n_boot,
                })

    out = pd.DataFrame(rows)
    out["Reconstructed q (24)"] = bh(
        out["Reconstructed two-sided sign p"]
    )
    return out


def validate(within, between, skewness):
    def single(df, **conditions):
        mask = np.ones(len(df), dtype=bool)
        for col, val in conditions.items():
            mask &= df[col].eq(val).to_numpy()
        rows = df.loc[mask]
        if len(rows) != 1:
            raise RuntimeError(
                f"Expected one row for {conditions}; found {len(rows)}"
            )
        return rows.iloc[0]

    row = single(
        within,
        Population="Datoga",
        Marker="DRD4",
        Trait="Verbal aggression",
    )
    exp = EXPECTED[("within", "Datoga", "DRD4", "Verbal aggression")]
    for col, key in (
        ("N0", "N0"),
        ("N1", "N1"),
        ("SD0", "SD0"),
        ("SD1", "SD1"),
        ("SD1/SD0", "ratio"),
        ("Brown-Forsythe p", "bf_p"),
    ):
        if abs(float(row[col]) - float(exp[key])) > 5e-10:
            raise RuntimeError(f"S7 within validation failed: {col}")

    for pop_pair, expected_key in (
        ("Armenians-Datoga", ("between", "Armenians-Datoga", "DRD4", "Physical aggression", "L+")),
        ("Armenians-Hadza", ("between", "Armenians-Hadza", "DRD4", "Physical aggression", "L+")),
    ):
        row = single(
            between,
            Marker="DRD4",
            Trait="Physical aggression",
            Category="L+",
            **{"Test type": "pair", "Comparison": pop_pair},
        )
        exp = EXPECTED[expected_key]
        if abs(float(row["SD B / SD A"]) - exp["ratio"]) > 5e-10:
            raise RuntimeError(f"S7 between ratio validation failed: {pop_pair}")
        if abs(float(row["Brown-Forsythe p"]) - exp["bf_p"]) > 5e-10:
            raise RuntimeError(f"S7 between p validation failed: {pop_pair}")

    for category in ("<=22", ">=23"):
        row = single(
            skewness,
            Population="Armenians",
            Marker="AR",
            Trait="Physical aggression",
            Category=category,
        )
        exp = EXPECTED[
            ("skew", "Armenians", "AR", "Physical aggression", category)
        ]
        if int(row["N"]) != exp["N"]:
            raise RuntimeError("S7 skew N validation failed")
        if abs(float(row["Fisher skewness"]) - exp["skew"]) > 5e-10:
            raise RuntimeError("S7 skew validation failed")
        if abs(float(row["skewtest p"]) - exp["p"]) > 5e-10:
            raise RuntimeError("S7 skewtest validation failed")

    print(
        "Validation PASSED: deterministic S7 targets reproduce the frozen analysis."
    )


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
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--bootstrap-resamples", type=int, default=2000)
    args = parser.parse_args()

    armenian = pd.read_excel(args.armenian, sheet_name="первичные данные")
    hd = pd.read_excel(
        args.hadza_datoga,
        sheet_name="Объединенные_мужчины",
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)

    tables = build_population_tables(armenian, hd)
    within = run_within(tables, args.seed, args.bootstrap_resamples)
    between = run_between(tables)
    skewness = run_skewness(tables, args.seed, args.bootstrap_resamples)
    skew_diff = run_skew_differences(
        tables, args.seed, args.bootstrap_resamples
    )

    within.to_csv(
        args.output_dir / "S7a_WithinPop_Variance_reconstructed.csv",
        index=False,
    )
    between.to_csv(
        args.output_dir / "S7b_BetweenPop_Variance_reconstructed.csv",
        index=False,
    )
    skewness.to_csv(
        args.output_dir / "S7c_Skewness_reconstructed.csv",
        index=False,
    )
    skew_diff.to_csv(
        args.output_dir / "S7d_Skewness_Differences_reconstructed.csv",
        index=False,
    )

    validate(within, between, skewness)
    print(
        "Bootstrap outputs are seeded reconstructions because the historical "
        "bootstrap state was not preserved."
    )


if __name__ == "__main__":
    main()
