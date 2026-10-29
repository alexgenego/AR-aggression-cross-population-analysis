#!/usr/bin/env python3
"""
Reconstructed reproducibility script for the Armenian single-locus and FDR analyses.

Manuscript:
"Androgen Receptor CAG Repeat Length and Aggressive Behavior Across Armenian,
Hadza, and Datoga Men: A Cross-Population Genetic Association Analysis"

This script reproduces the analysis families summarized in Supporting
Information S2a-S2d:

- S2a: 84 age-adjusted single-locus models
- S2b: omnibus tests of genotype/category models
- S2c: hierarchical FDR at the gene level
- S2d: hierarchical FDR at the individual-trait level

IMPORTANT
---------
This is a reconstructed reproducibility implementation based on the frozen
analytical workbook and final statistical specification. It is not claimed to
be the original historical source code.

Restricted participant-level data must remain outside the public repository.
By default the script expects the input workbook in ./data/, which is excluded
from Git version control.

Required local file
-------------------
data/Armenian_full_analysis.xlsx
  raw participant-level data on sheet: "первичные данные"

Frozen QC rules reproduced here
--------------------------------
- SEX: 1 = men, 2 = women
- age is included as a covariate in every single-locus model
- Ang = 50 is outside the permitted questionnaire range and is treated as missing
- strict multivariate phenodeviants are retained in the primary analysis and
  excluded only in the sensitivity analysis
- strict phenodeviant Excel rows in the frozen workbook: 69, 72, 125, 181
- HC3 heteroscedasticity-robust covariance estimates
- primary family: 7 variants × 6 traits × 2 sexes = 84 tests
- BH-FDR across all 84 primary tests
- hierarchical gene-level FDR uses the six strict/sensitivity p-values for
  each gene, combined by Simes; BH is then applied to seven genes separately
  in men and women

Primary single-locus coding
---------------------------
AR_CAG
  men: >=23 vs <=22 CAG repeats
  women: mean repeat length >=22 vs <=21.5
DRD4 exon III VNTR
  L+ = at least one allele with >=7 repeats; SS otherwise
HTR2C
  men: C- vs G-; women: C-allele dosage 0/1/2
OXTR rs53576
  A-allele dosage 0/1/2
HTR2A rs6311
  A-allele dosage 0/1/2 in the archived assay orientation
BDNF rs6265
  A-allele dosage 0/1/2
ERCC5 rs1047768
  T-allele dosage 0/1/2

Outputs
-------
outputs/S2a_Single_Locus.csv
outputs/S2b_Genotype_Tests.csv
outputs/S2c_FDR_Genes.csv
outputs/S2d_FDR_Traits.csv
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


TRAITS: Tuple[str, ...] = (
    "PhA",
    "VerbA",
    "Ang",
    "Host",
    "Proact_ag",
    "React_ag",
)

GENES: Tuple[str, ...] = (
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

PRIMARY_CODING_LABEL = {
    (1, "AR_CAG"): "AR: high category vs low category",
    (2, "AR_CAG"): "AR: mean >=22 vs <=21.5",
    (1, "DRD4eIII"): "DRD4: L+ vs SS",
    (2, "DRD4eIII"): "DRD4: L+ vs SS",
    (1, "HTR2C"): "HTR2C: C-allele dosage",
    (2, "HTR2C"): "HTR2C: C-allele dosage",
    (1, "OXTRrs53576"): "OXTR: A-allele dosage",
    (2, "OXTRrs53576"): "OXTR: A-allele dosage",
    (1, "HTR2Ars6311"): "HTR2A: A-allele dosage",
    (2, "HTR2Ars6311"): "HTR2A: A-allele dosage",
    (1, "BDNF rs6265"): "BDNF: A-allele dosage",
    (2, "BDNF rs6265"): "BDNF: A-allele dosage",
    (1, "ERCC5"): "ERCC5: T-allele dosage",
    (2, "ERCC5"): "ERCC5: T-allele dosage",
}

EXPECTED = {
    "AR_men_PhA_p": 0.004542930475003729,
    "AR_men_PhA_q84": 0.3816061599003132,
    "AR_men_PhA_strict_p": 0.016279575556002585,
    "DRD4_men_PhA_p": 0.015002651605563416,
    "BDNF_women_VerbA_p": 0.019573923052110954,
    "Omnibus_OXTR_women_PhA_p": 0.0010482265962477022,
    "Omnibus_OXTR_women_PhA_q": 0.10691911281726563,
    "Simes_HTR2A_men": 0.060204123923619804,
    "Simes_AR_men": 0.09767745333601552,
    "Gene_q_AR_men": 0.3418710866760543,
}


def bh(pvalues: Iterable[float]) -> np.ndarray:
    values = np.asarray(list(pvalues), dtype=float)
    return multipletests(values, alpha=0.05, method="fdr_bh")[1]


def simes(pvalues: Sequence[float]) -> float:
    ps = sorted(float(p) for p in pvalues)
    m = len(ps)
    return min(m * p / (i + 1) for i, p in enumerate(ps))


def parse_ar_alleles(value) -> Optional[List[float]]:
    if pd.isna(value):
        return None
    nums = [float(x) for x in re.findall(r"\d+(?:\.\d+)?", str(value).strip())]
    return nums or None


def parse_drd4_long_count(value) -> Optional[int]:
    if pd.isna(value):
        return None
    nums = [int(x) for x in re.findall(r"\d+", str(value).strip())]
    if len(nums) < 2:
        return None
    return int(nums[0] >= 7) + int(nums[1] >= 7)


def allele_dosage(value, allele: str) -> Optional[float]:
    if pd.isna(value):
        return None
    s = str(value).upper().replace("/", "").replace("-", "").replace(" ", "")
    if len(s) != 2 or any(ch not in "ACGT" for ch in s):
        return None
    return float(s.count(allele.upper()))


def clean_outcome(value, trait: str) -> Optional[float]:
    if pd.isna(value):
        return None
    x = float(value)
    if trait == "Ang" and x == 50:
        return None
    return x


def primary_genotype_value(row: pd.Series, gene: str, sex: int) -> Optional[float]:
    value = row.get(gene)

    if gene == "AR_CAG":
        alleles = parse_ar_alleles(value)
        if not alleles:
            return None
        if sex == 1:
            return float(alleles[0] >= 23)
        if len(alleles) < 2:
            return None
        mean_repeat = (alleles[0] + alleles[1]) / 2.0
        return float(mean_repeat >= 22)

    if gene == "DRD4eIII":
        n_long = parse_drd4_long_count(value)
        return None if n_long is None else float(n_long >= 1)

    if gene == "HTR2C":
        if pd.isna(value):
            return None
        s = str(value).upper().strip()
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


def collect_primary(
    df: pd.DataFrame,
    sex: int,
    gene: str,
    trait: str,
    exclude_strict: bool,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    y: List[float] = []
    g: List[float] = []
    age: List[float] = []

    for idx, row in df.iterrows():
        excel_row = int(idx) + 2
        if int(row["SEX"]) != sex:
            continue
        if exclude_strict and excel_row in STRICT_OUTLIER_EXCEL_ROWS:
            continue
        if pd.isna(row["AGE"]):
            continue

        outcome = clean_outcome(row[trait], trait)
        genotype = primary_genotype_value(row, gene, sex)
        if outcome is None or genotype is None:
            continue

        y.append(outcome)
        g.append(genotype)
        age.append(float(row["AGE"]))

    return np.asarray(y), np.asarray(g), np.asarray(age)


def partial_r2(y: np.ndarray, g: np.ndarray, age: np.ndarray) -> float:
    full = sm.OLS(y, np.column_stack([np.ones(len(y)), g, age])).fit()
    reduced = sm.OLS(y, np.column_stack([np.ones(len(y)), age])).fit()
    return float((full.rsquared - reduced.rsquared) / (1.0 - reduced.rsquared))


def fit_primary(
    df: pd.DataFrame,
    sex: int,
    gene: str,
    trait: str,
    exclude_strict: bool,
) -> Dict[str, float]:
    y, g, age = collect_primary(df, sex, gene, trait, exclude_strict)
    X = np.column_stack([np.ones(len(y)), g, age])
    model = sm.OLS(y, X).fit(cov_type="HC3")
    beta = float(model.params[1])
    se = float(model.bse[1])

    return {
        "N": int(len(y)),
        "beta": beta,
        "SE": se,
        "CI_low": beta - 1.96 * se,
        "CI_high": beta + 1.96 * se,
        "p": float(model.pvalues[1]),
        "partial_R2": partial_r2(y, g, age),
    }


def run_s2a(df: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict] = []
    for sex in (1, 2):
        for gene in GENES:
            for trait in TRAITS:
                full = fit_primary(df, sex, gene, trait, exclude_strict=False)
                strict = fit_primary(df, sex, gene, trait, exclude_strict=True)
                rows.append(
                    {
                        "Sex": SEX_LABEL[sex],
                        "Gene": gene,
                        "Coding": PRIMARY_CODING_LABEL[(sex, gene)],
                        "Trait": trait,
                        "N": full["N"],
                        "beta": full["beta"],
                        "SE": full["SE"],
                        "95% CI lower": full["CI_low"],
                        "95% CI upper": full["CI_high"],
                        "p": full["p"],
                        "partial R2": full["partial_R2"],
                        "N without strict phenodeviants": strict["N"],
                        "beta without strict phenodeviants": strict["beta"],
                        "p without strict phenodeviants": strict["p"],
                    }
                )

    out = pd.DataFrame(rows)
    out["q_BH_84"] = bh(out["p"])
    out["q_BH_84_without_strict"] = bh(out["p without strict phenodeviants"])
    return out


def ar_category(value, model: str, sex: int) -> Optional[str]:
    alleles = parse_ar_alleles(value)
    if not alleles:
        return None

    if sex == 1:
        x = alleles[0]
        if model == "AR 2 categories":
            return ">=23" if x >= 23 else "<=22"
        if model == "AR 3 categories":
            return "<=20" if x <= 20 else ("21-23" if x <= 23 else ">=24")
        raise KeyError(model)

    if len(alleles) < 2:
        return None
    a1, a2 = alleles[:2]
    amin = min(a1, a2)
    amean = (a1 + a2) / 2.0
    amax = max(a1, a2)

    if model == "AR min 2":
        return "min>=21" if amin >= 21 else "min<=20"
    if model == "AR min 3":
        return "min<=20" if amin <= 20 else ("min=21" if amin == 21 else "min>=22")
    if model == "AR mean 2":
        return "mean>=22" if amean >= 22 else "mean<=21.5"
    if model == "AR mean 3":
        return "mean<=21" if amean <= 21 else ("mean=21.5-22.5" if amean <= 22.5 else "mean>=23")
    if model == "AR max 2":
        return "max>=24" if amax >= 24 else "max<=23"
    if model == "AR max 3":
        return "max<=22" if amax <= 22 else ("max=23-24" if amax <= 24 else "max>=25")
    raise KeyError(model)


def category_value(row: pd.Series, gene: str, model: str, sex: int) -> Optional[str]:
    value = row.get(gene)
    if gene == "AR_CAG":
        return ar_category(value, model, sex)

    if gene == "DRD4eIII":
        n_long = parse_drd4_long_count(value)
        if n_long is None:
            return None
        if model == "DRD4 2 categories":
            return "L+" if n_long >= 1 else "SS"
        if model == "DRD4 3 categories":
            return "SS" if n_long == 0 else ("SL" if n_long == 1 else "LL")
        raise KeyError(model)

    if pd.isna(value):
        return None
    return str(value).strip()


OMNIBUS_MODELS = (
    (1, "AR_CAG", "AR 2 categories"),
    (1, "AR_CAG", "AR 3 categories"),
    (1, "BDNF rs6265", "original genotypes"),
    (1, "DRD4eIII", "DRD4 2 categories"),
    (1, "DRD4eIII", "DRD4 3 categories"),
    (1, "ERCC5", "original genotypes"),
    (1, "HTR2Ars6311", "original genotypes"),
    (1, "HTR2C", "original genotypes"),
    (1, "OXTRrs53576", "original genotypes"),
    (2, "AR_CAG", "AR max 2"),
    (2, "AR_CAG", "AR max 3"),
    (2, "AR_CAG", "AR mean 2"),
    (2, "AR_CAG", "AR mean 3"),
    (2, "AR_CAG", "AR min 2"),
    (2, "AR_CAG", "AR min 3"),
    (2, "BDNF rs6265", "original genotypes"),
    (2, "DRD4eIII", "DRD4 2 categories"),
    (2, "DRD4eIII", "DRD4 3 categories"),
    (2, "ERCC5", "original genotypes"),
    (2, "HTR2Ars6311", "original genotypes"),
    (2, "HTR2C", "original genotypes"),
    (2, "OXTRrs53576", "original genotypes"),
)


def collect_omnibus(
    df: pd.DataFrame,
    sex: int,
    gene: str,
    model_name: str,
    trait: str,
    exclude_strict: bool,
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    y: List[float] = []
    age: List[float] = []
    categories: List[str] = []

    for idx, row in df.iterrows():
        excel_row = int(idx) + 2
        if int(row["SEX"]) != sex:
            continue
        if exclude_strict and excel_row in STRICT_OUTLIER_EXCEL_ROWS:
            continue
        if pd.isna(row["AGE"]):
            continue

        outcome = clean_outcome(row[trait], trait)
        category = category_value(row, gene, model_name, sex)
        if outcome is None or category is None:
            continue

        y.append(outcome)
        age.append(float(row["AGE"]))
        categories.append(category)

    return np.asarray(y), np.asarray(age), categories


def fit_omnibus(
    df: pd.DataFrame,
    sex: int,
    gene: str,
    model_name: str,
    trait: str,
    exclude_strict: bool,
) -> Dict:
    y, age, cats = collect_omnibus(df, sex, gene, model_name, trait, exclude_strict)
    levels = sorted(set(cats))
    counts = {level: cats.count(level) for level in levels}
    min_group = min(counts.values()) if counts else 0

    if len(levels) < 2:
        return {
            "N": len(y),
            "counts": counts,
            "min_group": min_group,
            "p": np.nan,
            "partial_R2": np.nan,
        }

    dummy_columns = [
        np.asarray([1.0 if cat == level else 0.0 for cat in cats])
        for level in levels[1:]
    ]
    X = np.column_stack([np.ones(len(y)), age] + dummy_columns)
    robust = sm.OLS(y, X).fit(cov_type="HC3")

    k = len(levels) - 1
    R = np.zeros((k, X.shape[1]))
    for j in range(k):
        R[j, 2 + j] = 1.0

    try:
        test = robust.wald_test(R, scalar=True)
        p_value = float(test.pvalue)
    except Exception:
        p_value = np.nan

    full_ols = sm.OLS(y, X).fit()
    reduced_ols = sm.OLS(y, np.column_stack([np.ones(len(y)), age])).fit()
    pr2 = float(
        (full_ols.rsquared - reduced_ols.rsquared) / (1.0 - reduced_ols.rsquared)
    )

    return {
        "N": int(len(y)),
        "counts": counts,
        "min_group": int(min_group),
        "p": p_value,
        "partial_R2": pr2,
    }


def format_counts(counts: Dict[str, int]) -> str:
    return "; ".join(f"{k}: {v}" for k, v in sorted(counts.items()))


def run_s2b(df: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict] = []
    for sex, gene, model_name in OMNIBUS_MODELS:
        for trait in TRAITS:
            full = fit_omnibus(df, sex, gene, model_name, trait, exclude_strict=False)
            strict = fit_omnibus(df, sex, gene, model_name, trait, exclude_strict=True)
            admissible = full["min_group"] >= 5 and math.isfinite(float(full["p"]))
            rows.append(
                {
                    "Sex": SEX_LABEL[sex],
                    "Gene": gene,
                    "Category model": model_name,
                    "Trait": trait,
                    "N": full["N"],
                    "Group sizes": format_counts(full["counts"]),
                    "Minimum group": full["min_group"],
                    "Admissible for inference": "yes" if admissible else "no",
                    "Omnibus p": full["p"],
                    "partial R2": full["partial_R2"],
                    "p without strict phenodeviants": strict["p"],
                }
            )

    out = pd.DataFrame(rows)
    out["q_BH"] = np.nan
    mask = out["Admissible for inference"].eq("yes") & out["Omnibus p"].notna()
    out.loc[mask, "q_BH"] = bh(out.loc[mask, "Omnibus p"])
    return out.sort_values(
        ["Omnibus p", "Sex", "Gene", "Trait"], na_position="last"
    ).reset_index(drop=True)


def run_s2c_s2d(s2a: pd.DataFrame):
    gene_rows: List[Dict] = []

    for sex in ("men", "women"):
        for gene in GENES:
            sub = s2a[(s2a["Sex"] == sex) & (s2a["Gene"] == gene)].copy()
            ps = sub["p without strict phenodeviants"].astype(float).tolist()
            global_p = simes(ps)
            best_idx = sub["p without strict phenodeviants"].astype(float).idxmin()
            best = sub.loc[best_idx]
            gene_rows.append(
                {
                    "Sex": sex,
                    "Gene": gene,
                    "Global Simes p": global_p,
                    "Minimum trait p": float(best["p without strict phenodeviants"]),
                    "Trait with minimum p": best["Trait"],
                }
            )

    genes = pd.DataFrame(gene_rows)
    genes["BH q among 7 genes"] = np.nan
    for sex in ("men", "women"):
        mask = genes["Sex"].eq(sex)
        genes.loc[mask, "BH q among 7 genes"] = bh(
            genes.loc[mask, "Global Simes p"]
        )

    genes["Selected at level 1"] = np.where(
        genes["BH q among 7 genes"] < 0.05, "yes", "no"
    )
    genes["Decision"] = np.where(
        genes["Selected at level 1"].eq("yes"),
        "proceed to within-gene level",
        "stopped at gene level",
    )
    genes = genes.sort_values(["Sex", "Global Simes p"]).reset_index(drop=True)

    trait_rows: List[Dict] = []
    for _, row in s2a.iterrows():
        gene_passed = genes.loc[
            (genes["Sex"] == row["Sex"]) & (genes["Gene"] == row["Gene"]),
            "Selected at level 1",
        ].iloc[0]

        trait_rows.append(
            {
                "Sex": row["Sex"],
                "Gene": row["Gene"],
                "Trait": row["Trait"],
                "Coding": row["Coding"],
                "N": row["N without strict phenodeviants"],
                "strict beta": row["beta without strict phenodeviants"],
                "strict p": row["p without strict phenodeviants"],
                "Gene passed level 1": gene_passed,
                "BB alpha within family": 0.0 if gene_passed == "no" else np.nan,
                "within-family q": np.nan,
                "Selected at level 2": "no",
            }
        )

    return genes, pd.DataFrame(trait_rows)


def validate(s2a: pd.DataFrame, s2b: pd.DataFrame, s2c: pd.DataFrame):
    def one(df, **conds):
        mask = np.ones(len(df), dtype=bool)
        for key, value in conds.items():
            mask &= df[key].eq(value).to_numpy()
        matches = df.loc[mask]
        if len(matches) != 1:
            raise RuntimeError(f"Expected one row for {conds}, found {len(matches)}")
        return matches.iloc[0]

    ar = one(s2a, Sex="men", Gene="AR_CAG", Trait="PhA")
    dr = one(s2a, Sex="men", Gene="DRD4eIII", Trait="PhA")
    bd = one(s2a, Sex="women", Gene="BDNF rs6265", Trait="VerbA")
    ox = one(
        s2b,
        Sex="women",
        Gene="OXTRrs53576",
        **{"Category model": "original genotypes", "Trait": "PhA"},
    )
    htr_gene = one(s2c, Sex="men", Gene="HTR2Ars6311")
    ar_gene = one(s2c, Sex="men", Gene="AR_CAG")

    checks = {
        "AR_men_PhA_p": float(ar["p"]),
        "AR_men_PhA_q84": float(ar["q_BH_84"]),
        "AR_men_PhA_strict_p": float(ar["p without strict phenodeviants"]),
        "DRD4_men_PhA_p": float(dr["p"]),
        "BDNF_women_VerbA_p": float(bd["p"]),
        "Omnibus_OXTR_women_PhA_p": float(ox["Omnibus p"]),
        "Omnibus_OXTR_women_PhA_q": float(ox["q_BH"]),
        "Simes_HTR2A_men": float(htr_gene["Global Simes p"]),
        "Simes_AR_men": float(ar_gene["Global Simes p"]),
        "Gene_q_AR_men": float(ar_gene["BH q among 7 genes"]),
    }

    failed = []
    for name, observed in checks.items():
        expected = EXPECTED[name]
        if abs(observed - expected) > 5e-10:
            failed.append((name, observed, expected))

    if failed:
        message = "\n".join(
            f"{name}: observed={obs!r}, expected={exp!r}"
            for name, obs, exp in failed
        )
        raise RuntimeError("Validation against frozen S2 targets failed:\n" + message)

    if int((s2a["q_BH_84"] < 0.05).sum()) != 0:
        raise RuntimeError("Unexpected FDR-confirmed single-locus result in S2a.")
    if int((s2c["BH q among 7 genes"] < 0.05).sum()) != 0:
        raise RuntimeError("Unexpected level-1 hierarchical-FDR gene selection.")

    print("Validation PASSED: frozen S2a-S2d targets reproduced.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--armenian",
        type=Path,
        default=Path("data/Armenian_full_analysis.xlsx"),
        help="Restricted Armenian workbook containing sheet 'первичные данные'.",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()

    df = pd.read_excel(args.armenian, sheet_name="первичные данные")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    s2a = run_s2a(df)
    s2b = run_s2b(df)
    s2c, s2d = run_s2c_s2d(s2a)

    s2a.to_csv(args.output_dir / "S2a_Single_Locus.csv", index=False)
    s2b.to_csv(args.output_dir / "S2b_Genotype_Tests.csv", index=False)
    s2c.to_csv(args.output_dir / "S2c_FDR_Genes.csv", index=False)
    s2d.to_csv(args.output_dir / "S2d_FDR_Traits.csv", index=False)

    validate(s2a, s2b, s2c)

    print(f"Outputs written to: {args.output_dir.resolve()}")
    print("No primary single-locus association survives BH-FDR across the 84-test family.")


if __name__ == "__main__":
    main()
