#!/usr/bin/env python3
"""
Reconstructed reproducibility script for the cross-population shared-locus analysis.

This script reproduces the analyses summarized in Supporting Information S6a-S6d
for the manuscript:

"Androgen Receptor CAG Repeat Length and Aggressive Behavior Across Armenian,
Hadza, and Datoga Men: A Cross-Population Genetic Association Analysis"

IMPORTANT
---------
This is a reconstructed reproducibility implementation based on the frozen
analytical workbooks and the final statistical specification. It is not claimed
to be the original historical analysis source code.

Restricted participant-level source files must remain outside the public
repository. By default the script expects them in ./data/, which is gitignored.

Required local files
--------------------
1. data/Armenian_full_analysis.xlsx
   - raw participant-level data on sheet: "первичные данные"
2. data/Hadza_Datoga_adults_enriched_by_name_v2.xlsx
   - participant-level archive on sheet: "Adults_enriched"

Outputs
-------
outputs/S6a_ThreePop_Primary.csv
outputs/S6b_ThreeLocus_Block.csv
outputs/S6c_AR_Continuous.csv
outputs/S6d_AgeSensitivity_All12.csv
outputs/S6d_AgeSensitivity_Selected.csv
outputs/S6f_Figure1_AR_SourceData.csv

The script uses HC3 heteroscedasticity-robust covariance estimates and
Benjamini-Hochberg FDR correction.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


POPULATIONS = ("Hadza", "Datoga", "Armenians")
TRAITS = (
    ("physical_aggression", "Physical aggression"),
    ("verbal_aggression", "Verbal aggression"),
    ("anger", "Anger"),
    ("hostility", "Hostility"),
)
ARM_TRAIT = {
    "physical_aggression": "PhA",
    "verbal_aggression": "VerbA",
    "anger": "Ang",
    "hostility": "Host",
}

LOCUS_LABELS = {
    "AR": "AR ≥23 vs ≤22",
    "DRD4": "DRD4 L+ vs SS",
    "HTR2A": "HTR2A A-allele dosage",
}

EXPECTED = {
    "S6a_AR_PhA_N": 380,
    "S6a_AR_PhA_beta": -0.3924666663,
    "S6a_AR_PhA_p": 0.0002691628,
    "S6a_AR_PhA_q": 0.0032299539,
    "S6a_AR_PhA_heterogeneity_p": 0.3862528293,
    "S6b_PhA_block_p": 0.0052147507,
    "S6b_PhA_block_q": 0.0208590028,
    "S6c_PhA_beta": -0.0580165486,
    "S6c_PhA_p": 0.0021793401,
    "S6d_AR_PhA_beta": -0.3879738610,
    "S6d_AR_PhA_p": 0.0003754794,
    "S6d_Block_PhA_p": 0.0097250733,
}


def bh(pvalues: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(pvalues), dtype=float)
    return multipletests(p, alpha=0.05, method="fdr_bh")[1]


def parse_ar(value) -> Optional[float]:
    """Parse male AR CAG repeat number from the enriched archive."""
    if pd.isna(value):
        return None
    if isinstance(value, (int, float, np.integer, np.floating)):
        x = float(value)
        return x if 10 <= x <= 40 else None
    s = str(value).strip()
    if not s or "/" in s:
        # Slash values can denote paired alleles / non-male-style entries.
        return None
    m = re.search(r"(\d+)", s)
    if not m:
        return None
    x = float(m.group(1))
    return x if 10 <= x <= 40 else None


def parse_arm_ar(value) -> Optional[float]:
    if pd.isna(value):
        return None
    try:
        x = float(str(value).strip())
    except ValueError:
        return None
    return x if 10 <= x <= 40 else None


def parse_drd4_lplus(value) -> Optional[float]:
    """
    DRD4 exon III VNTR coding:
    SS = 0 if both alleles <7 repeats
    L+ = 1 if at least one allele has >=7 repeats.

    Historical Excel files occasionally render genotypes such as 7/7 or 5/7
    as calendar dates. Any July date is therefore interpreted as containing a
    7-repeat allele.
    """
    if pd.isna(value):
        return None

    if isinstance(value, (pd.Timestamp,)):
        return 1.0 if value.month == 7 else None

    # Python datetime/date objects can also appear depending on the Excel reader.
    if hasattr(value, "month") and hasattr(value, "day"):
        try:
            return 1.0 if int(value.month) == 7 else None
        except Exception:
            pass

    s = str(value).strip().lower()
    if not s:
        return None

    if "июл" in s or "jul" in s:
        return 1.0

    nums = [int(x) for x in re.findall(r"\d+", s)]
    plausible = [x for x in nums if 2 <= x <= 15]
    if len(plausible) < 2:
        return None
    a, b = plausible[:2]
    return float(max(a, b) >= 7)


def parse_htr2a_dosage(value) -> Optional[float]:
    """A-allele dosage from archived A/G genotype labels."""
    if pd.isna(value):
        return None
    s = str(value).upper().replace("/", "").replace(" ", "")
    if len(s) != 2 or any(ch not in "AG" for ch in s):
        return None
    return float(s.count("A"))


def load_data(armenian_file: Path, hd_file: Path):
    arm = pd.read_excel(armenian_file, sheet_name="первичные данные")
    hd = pd.read_excel(hd_file, sheet_name="Adults_enriched", header=3)
    return arm, hd


def hd_is_eligible_adult_male(row: pd.Series) -> bool:
    return (
        row.get("sex") == 1
        and row.get("age_status") in ("18+", "Age missing")
        and row.get("population") in ("Hadza", "Datoga")
    )


def phenotype_reference_stats(arm: pd.DataFrame, hd: pd.DataFrame, trait: str):
    """
    Standardization used in the frozen three-population analysis:
    phenotype is z-standardized within population using ALL eligible male
    phenotype observations for that population, independent of genotype
    availability. Sample SD (ddof=1) is used.
    """
    out = {}

    arm_col = ARM_TRAIT[trait]
    x = pd.to_numeric(
        arm.loc[(arm["SEX"] == 1) & arm[arm_col].notna(), arm_col],
        errors="coerce",
    ).dropna()
    out["Armenians"] = (float(x.mean()), float(x.std(ddof=1)), int(x.shape[0]))

    for pop in ("Hadza", "Datoga"):
        mask = hd.apply(hd_is_eligible_adult_male, axis=1)
        x = pd.to_numeric(
            hd.loc[mask & (hd["population"] == pop) & hd[trait].notna(), trait],
            errors="coerce",
        ).dropna()
        out[pop] = (float(x.mean()), float(x.std(ddof=1)), int(x.shape[0]))

    return out


def locus_value_arm(row: pd.Series, locus: str) -> Optional[float]:
    if locus == "AR":
        x = parse_arm_ar(row.get("AR_CAG"))
        return None if x is None else float(x >= 23)
    if locus == "DRD4":
        return parse_drd4_lplus(row.get("DRD4eIII"))
    if locus == "HTR2A":
        return parse_htr2a_dosage(row.get("HTR2Ars6311"))
    raise KeyError(locus)


def locus_value_hd(row: pd.Series, locus: str) -> Optional[float]:
    if locus == "AR":
        x = parse_ar(row.get("AR2_final"))
        return None if x is None else float(x >= 23)
    if locus == "DRD4":
        return parse_drd4_lplus(row.get("DRD4E3_final"))
    if locus == "HTR2A":
        x = row.get("HTR2A_A_dosage_adult_strand")
        if pd.isna(x):
            return parse_htr2a_dosage(row.get("HTR2A_final"))
        return float(x)
    raise KeyError(locus)


def ar_repeat_arm(row: pd.Series) -> Optional[float]:
    return parse_arm_ar(row.get("AR_CAG"))


def ar_repeat_hd(row: pd.Series) -> Optional[float]:
    return parse_ar(row.get("AR2_final"))


def collect_single_locus(
    arm: pd.DataFrame,
    hd: pd.DataFrame,
    trait: str,
    locus: str,
    continuous_ar: bool = False,
) -> List[Dict]:
    stats = phenotype_reference_stats(arm, hd, trait)
    records: List[Dict] = []

    arm_col = ARM_TRAIT[trait]
    for _, row in arm.iterrows():
        if row.get("SEX") != 1 or pd.isna(row.get(arm_col)):
            continue
        g = ar_repeat_arm(row) if continuous_ar else locus_value_arm(row, locus)
        if g is None:
            continue
        mu, sd, _ = stats["Armenians"]
        records.append(
            {
                "population": "Armenians",
                "g": float(g),
                "y_z": (float(row[arm_col]) - mu) / sd,
                "age": None if pd.isna(row.get("AGE")) else float(row.get("AGE")),
            }
        )

    for _, row in hd.iterrows():
        if not hd_is_eligible_adult_male(row) or pd.isna(row.get(trait)):
            continue
        pop = row.get("population")
        g = ar_repeat_hd(row) if continuous_ar else locus_value_hd(row, locus)
        if g is None:
            continue
        mu, sd, _ = stats[pop]
        records.append(
            {
                "population": pop,
                "g": float(g),
                "y_z": (float(row[trait]) - mu) / sd,
                "age": None if pd.isna(row.get("age")) else float(row.get("age")),
            }
        )

    return records


def collect_three_locus(arm: pd.DataFrame, hd: pd.DataFrame, trait: str) -> List[Dict]:
    stats = phenotype_reference_stats(arm, hd, trait)
    records: List[Dict] = []
    arm_col = ARM_TRAIT[trait]

    for _, row in arm.iterrows():
        if row.get("SEX") != 1 or pd.isna(row.get(arm_col)):
            continue
        ar_count = ar_repeat_arm(row)
        dr = locus_value_arm(row, "DRD4")
        ht = locus_value_arm(row, "HTR2A")
        if ar_count is None or dr is None or ht is None:
            continue
        mu, sd, _ = stats["Armenians"]
        records.append(
            {
                "population": "Armenians",
                "ar": float(ar_count >= 23),
                "drd4": float(dr),
                "htr2a": float(ht),
                "y_z": (float(row[arm_col]) - mu) / sd,
                "age": None if pd.isna(row.get("AGE")) else float(row.get("AGE")),
            }
        )

    for _, row in hd.iterrows():
        if not hd_is_eligible_adult_male(row) or pd.isna(row.get(trait)):
            continue
        pop = row.get("population")
        ar_count = ar_repeat_hd(row)
        dr = locus_value_hd(row, "DRD4")
        ht = locus_value_hd(row, "HTR2A")
        if ar_count is None or dr is None or ht is None:
            continue
        mu, sd, _ = stats[pop]
        records.append(
            {
                "population": pop,
                "ar": float(ar_count >= 23),
                "drd4": float(dr),
                "htr2a": float(ht),
                "y_z": (float(row[trait]) - mu) / sd,
                "age": None if pd.isna(row.get("age")) else float(row.get("age")),
            }
        )

    return records


def population_counts(records: List[Dict]) -> Dict[str, int]:
    return {p: sum(r["population"] == p for r in records) for p in POPULATIONS}


def fit_population_specific(records: List[Dict]) -> Dict[str, float]:
    effects = {}
    for pop in POPULATIONS:
        sub = [r for r in records if r["population"] == pop]
        X = sm.add_constant(np.array([r["g"] for r in sub], dtype=float))
        y = np.array([r["y_z"] for r in sub], dtype=float)
        model = sm.OLS(y, X).fit(cov_type="HC3")
        effects[pop] = float(model.params[1])
    return effects


def fit_common_and_heterogeneity(records: List[Dict]):
    y = np.array([r["y_z"] for r in records], dtype=float)

    X = []
    for r in records:
        d = float(r["population"] == "Datoga")
        a = float(r["population"] == "Armenians")
        X.append([1.0, r["g"], d, a])
    common = sm.OLS(y, np.asarray(X, dtype=float)).fit(cov_type="HC3")

    Xh = []
    for r in records:
        d = float(r["population"] == "Datoga")
        a = float(r["population"] == "Armenians")
        Xh.append([1.0, r["g"], d, a, r["g"] * d, r["g"] * a])
    hetero = sm.OLS(y, np.asarray(Xh, dtype=float)).fit(cov_type="HC3")
    R = np.zeros((2, 6))
    R[0, 4] = 1
    R[1, 5] = 1
    w = hetero.wald_test(R, scalar=True)

    return common, float(w.statistic), float(w.pvalue)


def fit_three_locus(records: List[Dict], age_adjusted: bool = False):
    counts = population_counts(records)
    y = np.array([r["y_z"] for r in records], dtype=float)

    age_stats = None
    if age_adjusted:
        age_stats = {}
        for pop in POPULATIONS:
            ages = np.array(
                [r["age"] for r in records if r["population"] == pop and r["age"] is not None],
                dtype=float,
            )
            age_stats[pop] = (float(ages.mean()), float(ages.std(ddof=1)))

    X = []
    for r in records:
        pop = r["population"]
        d = float(pop == "Datoga")
        a = float(pop == "Armenians")
        base = [1.0, r["ar"], r["drd4"], r["htr2a"], d, a]
        if age_adjusted:
            missing = float(r["age"] is None)
            az = 0.0 if missing else (r["age"] - age_stats[pop][0]) / age_stats[pop][1]
            base += [
                az if pop == "Hadza" else 0.0,
                az if pop == "Datoga" else 0.0,
                az if pop == "Armenians" else 0.0,
                missing,
            ]
        X.append(base)

    model = sm.OLS(y, np.asarray(X, dtype=float)).fit(cov_type="HC3")
    p = len(X[0])
    R = np.zeros((3, p))
    R[0, 1] = 1
    R[1, 2] = 1
    R[2, 3] = 1
    block = model.wald_test(R, scalar=True)

    Xh = []
    for r in records:
        pop = r["population"]
        d = float(pop == "Datoga")
        a = float(pop == "Armenians")
        base = [1.0, r["ar"], r["drd4"], r["htr2a"], d, a]
        if age_adjusted:
            missing = float(r["age"] is None)
            az = 0.0 if missing else (r["age"] - age_stats[pop][0]) / age_stats[pop][1]
            base += [
                az if pop == "Hadza" else 0.0,
                az if pop == "Datoga" else 0.0,
                az if pop == "Armenians" else 0.0,
                missing,
            ]
        base += [
            r["ar"] * d,
            r["ar"] * a,
            r["drd4"] * d,
            r["drd4"] * a,
            r["htr2a"] * d,
            r["htr2a"] * a,
        ]
        Xh.append(base)

    hetero = sm.OLS(y, np.asarray(Xh, dtype=float)).fit(cov_type="HC3")
    start = len(X[0])
    R2 = np.zeros((6, len(Xh[0])))
    for i in range(6):
        R2[i, start + i] = 1
    iw = hetero.wald_test(R2, scalar=True)

    return model, float(block.statistic), float(block.pvalue), float(iw.statistic), float(iw.pvalue), counts


def fit_age_adjusted_single(records: List[Dict]):
    """
    Age sensitivity specification:
    - age z-standardized within each population in the CURRENT analytic set,
      using non-missing ages and sample SD (ddof=1);
    - missing age assigned z=0;
    - one missing-age indicator;
    - population-specific age slopes;
    - population indicators retained.
    """
    age_stats = {}
    for pop in POPULATIONS:
        ages = np.array(
            [r["age"] for r in records if r["population"] == pop and r["age"] is not None],
            dtype=float,
        )
        age_stats[pop] = (float(ages.mean()), float(ages.std(ddof=1)))

    y = np.array([r["y_z"] for r in records], dtype=float)
    X = []
    for r in records:
        pop = r["population"]
        d = float(pop == "Datoga")
        a = float(pop == "Armenians")
        missing = float(r["age"] is None)
        az = 0.0 if missing else (r["age"] - age_stats[pop][0]) / age_stats[pop][1]
        X.append(
            [
                1.0,
                r["g"],
                d,
                a,
                az if pop == "Hadza" else 0.0,
                az if pop == "Datoga" else 0.0,
                az if pop == "Armenians" else 0.0,
                missing,
            ]
        )
    model = sm.OLS(y, np.asarray(X, dtype=float)).fit(cov_type="HC3")

    Xh = []
    for r in records:
        pop = r["population"]
        d = float(pop == "Datoga")
        a = float(pop == "Armenians")
        missing = float(r["age"] is None)
        az = 0.0 if missing else (r["age"] - age_stats[pop][0]) / age_stats[pop][1]
        Xh.append(
            [
                1.0,
                r["g"],
                d,
                a,
                az if pop == "Hadza" else 0.0,
                az if pop == "Datoga" else 0.0,
                az if pop == "Armenians" else 0.0,
                missing,
                r["g"] * d,
                r["g"] * a,
            ]
        )
    hetero = sm.OLS(y, np.asarray(Xh, dtype=float)).fit(cov_type="HC3")
    R = np.zeros((2, 10))
    R[0, 8] = 1
    R[1, 9] = 1
    w = hetero.wald_test(R, scalar=True)
    return model, float(w.statistic), float(w.pvalue)


def run_s6a(arm, hd) -> pd.DataFrame:
    rows = []
    for locus in ("AR", "DRD4", "HTR2A"):
        for trait, trait_label in TRAITS:
            records = collect_single_locus(arm, hd, trait, locus)
            common, hchi, hp = fit_common_and_heterogeneity(records)
            eff = fit_population_specific(records)
            counts = population_counts(records)
            ci = common.conf_int()[1]
            rows.append(
                {
                    "Locus / coding": LOCUS_LABELS[locus],
                    "Trait": trait_label,
                    "Total N": len(records),
                    "Hadza N": counts["Hadza"],
                    "Datoga N": counts["Datoga"],
                    "Armenian N": counts["Armenians"],
                    "Common beta, SD": common.params[1],
                    "95% CI lower": ci[0],
                    "95% CI upper": ci[1],
                    "p": common.pvalues[1],
                    "Heterogeneity chi2 (2 df)": hchi,
                    "Heterogeneity p": hp,
                    "Hadza beta": eff["Hadza"],
                    "Datoga beta": eff["Datoga"],
                    "Armenian beta": eff["Armenians"],
                }
            )
    out = pd.DataFrame(rows)
    out["BH-FDR q (12)"] = bh(out["p"])
    out["Heterogeneity q (12)"] = bh(out["Heterogeneity p"])
    return out


def run_s6b(arm, hd) -> pd.DataFrame:
    rows = []
    for trait, trait_label in TRAITS:
        records = collect_three_locus(arm, hd, trait)
        model, chi, p, ichi, ip, counts = fit_three_locus(records, age_adjusted=False)
        rows.append(
            {
                "Trait": trait_label,
                "Total N": len(records),
                "Hadza N": counts["Hadza"],
                "Datoga N": counts["Datoga"],
                "Armenian N": counts["Armenians"],
                "Block chi2 (3 df)": chi,
                "Block p": p,
                "AR beta": model.params[1],
                "AR p": model.pvalues[1],
                "DRD4 beta": model.params[2],
                "DRD4 p": model.pvalues[2],
                "HTR2A beta": model.params[3],
                "HTR2A p": model.pvalues[3],
                "Interaction chi2 (6 df)": ichi,
                "Interaction p": ip,
            }
        )
    out = pd.DataFrame(rows)
    out["Block q (4)"] = bh(out["Block p"])
    out["Interaction q (4)"] = bh(out["Interaction p"])
    return out


def run_s6c(arm, hd) -> pd.DataFrame:
    rows = []
    for trait, trait_label in TRAITS:
        records = collect_single_locus(arm, hd, trait, "AR", continuous_ar=True)
        common, hchi, hp = fit_common_and_heterogeneity(records)
        eff = fit_population_specific(records)
        counts = population_counts(records)
        ci = common.conf_int()[1]
        rows.append(
            {
                "Trait": trait_label,
                "Total N": len(records),
                "Hadza N": counts["Hadza"],
                "Datoga N": counts["Datoga"],
                "Armenian N": counts["Armenians"],
                "Common beta per CAG repeat, SD": common.params[1],
                "95% CI lower": ci[0],
                "95% CI upper": ci[1],
                "p": common.pvalues[1],
                "Heterogeneity chi2": hchi,
                "Heterogeneity p": hp,
                "Hadza beta": eff["Hadza"],
                "Datoga beta": eff["Datoga"],
                "Armenian beta": eff["Armenians"],
            }
        )
    out = pd.DataFrame(rows)
    out["BH-FDR q (4)"] = bh(out["p"])
    out["Heterogeneity q (4)"] = bh(out["Heterogeneity p"])
    return out


def run_s6d(arm, hd):
    all_rows = []
    for locus in ("AR", "DRD4", "HTR2A"):
        for trait, trait_label in TRAITS:
            records = collect_single_locus(arm, hd, trait, locus)
            model, _, hp = fit_age_adjusted_single(records)
            counts = population_counts(records)
            ci = model.conf_int()[1]
            all_rows.append(
                {
                    "Locus / coding": LOCUS_LABELS[locus],
                    "Trait": trait_label,
                    "Total N": len(records),
                    "Hadza N": counts["Hadza"],
                    "Datoga N": counts["Datoga"],
                    "Armenian N": counts["Armenians"],
                    "beta": model.params[1],
                    "95% CI lower": ci[0],
                    "95% CI upper": ci[1],
                    "p": model.pvalues[1],
                    "Heterogeneity p": hp,
                }
            )

    all_out = pd.DataFrame(all_rows)
    all_out["q12"] = bh(all_out["p"])
    all_out["Heterogeneity q12"] = bh(all_out["Heterogeneity p"])

    # Age-adjusted three-locus block across four traits; BH q across the four block tests.
    block_rows = []
    for trait, trait_label in TRAITS:
        records = collect_three_locus(arm, hd, trait)
        model, chi, p, _, ip, counts = fit_three_locus(records, age_adjusted=True)
        block_rows.append(
            {
                "Trait": trait_label,
                "Total N": len(records),
                "Hadza N": counts["Hadza"],
                "Datoga N": counts["Datoga"],
                "Armenian N": counts["Armenians"],
                "Block chi2": chi,
                "Block p": p,
                "AR p": model.pvalues[1],
                "Interaction p": ip,
            }
        )
    block = pd.DataFrame(block_rows)
    block["Block q (4)"] = bh(block["Block p"])
    block["Interaction q (4)"] = bh(block["Interaction p"])

    # Compact table matching the selected age-sensitivity rows emphasized in S6d.
    keep = (
        ((all_out["Locus / coding"] == LOCUS_LABELS["AR"]))
        | (
            (all_out["Locus / coding"].isin([LOCUS_LABELS["DRD4"], LOCUS_LABELS["HTR2A"]]))
            & (all_out["Trait"] == "Physical aggression")
        )
    )
    selected = all_out.loc[keep].copy()
    ph_block = block.loc[block["Trait"] == "Physical aggression"].iloc[0]
    selected = pd.concat(
        [
            selected,
            pd.DataFrame(
                [
                    {
                        "Locus / coding": "Three-locus block",
                        "Trait": "Physical aggression",
                        "Total N": ph_block["Total N"],
                        "Hadza N": ph_block["Hadza N"],
                        "Datoga N": ph_block["Datoga N"],
                        "Armenian N": ph_block["Armenian N"],
                        "beta": ph_block["Block chi2"],
                        "95% CI lower": np.nan,
                        "95% CI upper": np.nan,
                        "p": ph_block["Block p"],
                        "Heterogeneity p": ph_block["Interaction p"],
                        "q12": ph_block["Block q (4)"],
                        "Heterogeneity q12": ph_block["Interaction q (4)"],
                    }
                ]
            ),
        ],
        ignore_index=True,
    )
    return all_out, selected, block


def validate(s6a, s6b, s6c, s6d_selected):
    checks = [
        ("S6a_AR_PhA_N", int(s6a.iloc[0]["Total N"]), 0),
        ("S6a_AR_PhA_beta", float(s6a.iloc[0]["Common beta, SD"]), 5e-10),
        ("S6a_AR_PhA_p", float(s6a.iloc[0]["p"]), 5e-10),
        ("S6a_AR_PhA_q", float(s6a.iloc[0]["BH-FDR q (12)"]), 5e-10),
        ("S6a_AR_PhA_heterogeneity_p", float(s6a.iloc[0]["Heterogeneity p"]), 5e-10),
        ("S6b_PhA_block_p", float(s6b.iloc[0]["Block p"]), 5e-10),
        ("S6b_PhA_block_q", float(s6b.iloc[0]["Block q (4)"]), 5e-10),
        ("S6c_PhA_beta", float(s6c.iloc[0]["Common beta per CAG repeat, SD"]), 5e-10),
        ("S6c_PhA_p", float(s6c.iloc[0]["p"]), 5e-10),
        ("S6d_AR_PhA_beta", float(s6d_selected.iloc[0]["beta"]), 5e-10),
        ("S6d_AR_PhA_p", float(s6d_selected.iloc[0]["p"]), 5e-10),
        ("S6d_Block_PhA_p", float(s6d_selected.iloc[-1]["p"]), 5e-10),
    ]
    failed = []
    for name, observed, tol in checks:
        expected = EXPECTED[name]
        ok = observed == expected if tol == 0 else abs(observed - expected) <= tol
        if not ok:
            failed.append((name, observed, expected))
    if failed:
        msg = "\n".join(f"{n}: observed={o!r}, expected={e!r}" for n, o, e in failed)
        raise RuntimeError("Validation against the frozen analytical targets failed:\n" + msg)
    print("Validation PASSED: frozen S6 targets reproduced.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--armenian",
        type=Path,
        default=Path("data/Armenian_full_analysis.xlsx"),
        help="Restricted Armenian analysis workbook containing sheet 'первичные данные'.",
    )
    ap.add_argument(
        "--hadza-datoga",
        type=Path,
        default=Path("data/Hadza_Datoga_adults_enriched_by_name_v2.xlsx"),
        help="Restricted Hadza/Datoga enriched adult archive.",
    )
    ap.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = ap.parse_args()

    arm, hd = load_data(args.armenian, args.hadza_datoga)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    s6a = run_s6a(arm, hd)
    s6b = run_s6b(arm, hd)
    s6c = run_s6c(arm, hd)
    s6d_all, s6d_selected, s6d_block = run_s6d(arm, hd)

    s6a.to_csv(args.output_dir / "S6a_ThreePop_Primary.csv", index=False)
    s6b.to_csv(args.output_dir / "S6b_ThreeLocus_Block.csv", index=False)
    s6c.to_csv(args.output_dir / "S6c_AR_Continuous.csv", index=False)
    s6d_all.to_csv(args.output_dir / "S6d_AgeSensitivity_All12.csv", index=False)
    s6d_selected.to_csv(args.output_dir / "S6d_AgeSensitivity_Selected.csv", index=False)

    figure1 = s6a.loc[
        (s6a["Locus / coding"] == LOCUS_LABELS["AR"])
        & (s6a["Trait"] == "Physical aggression"),
        ["Hadza beta", "Datoga beta", "Armenian beta"],
    ]
    # Recreate population-specific point estimates and their HC3 intervals/p values.
    fig_rows = []
    records = collect_single_locus(arm, hd, "physical_aggression", "AR")
    for pop in POPULATIONS:
        sub = [r for r in records if r["population"] == pop]
        X = sm.add_constant(np.array([r["g"] for r in sub], dtype=float))
        y = np.array([r["y_z"] for r in sub], dtype=float)
        model = sm.OLS(y, X).fit(cov_type="HC3")
        ci = model.conf_int()[1]
        fig_rows.append(
            {
                "Trait": "Physical aggression",
                "Population": pop,
                "Expanded N": len(sub),
                "Expanded beta (SD)": model.params[1],
                "CI low": ci[0],
                "CI high": ci[1],
                "p": model.pvalues[1],
            }
        )
    pd.DataFrame(fig_rows).to_csv(
        args.output_dir / "S6f_Figure1_AR_SourceData.csv", index=False
    )

    validate(s6a, s6b, s6c, s6d_selected)

    print(f"Outputs written to: {args.output_dir.resolve()}")
    print("Primary AR physical-aggression result:")
    print(s6a.iloc[0].to_string())


if __name__ == "__main__":
    main()
