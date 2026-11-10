#!/usr/bin/env python3
"""
Recompute TMPRSS5 (Spinesin) pharmacodynamic statistics from GSE103799
sciatic-nerve TPM values in CMT1A (C22) and wild-type mice.
"""

from __future__ import annotations

import csv
import math
import os
from collections import defaultdict
from decimal import ROUND_HALF_EVEN, Decimal

from scipy import stats


CSV_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data",
    "gse103799_tmprss5_tpm.csv",
)

GROUPS = ("WT_PBS", "C22_PBS", "C22_ASO")


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def sample_sd(values: list[float]) -> float:
    n = len(values)
    if n < 2:
        return float("nan")
    m = mean(values)
    return math.sqrt(sum((x - m) ** 2 for x in values) / (n - 1))


def round_half_even(value: float, places: int = 2) -> float:
    """Manuscript-style half-even rounding via a 3-dp intermediate."""
    quant_fine = Decimal("0.001")
    quant = Decimal("1").scaleb(-places)
    stepped = Decimal(str(value)).quantize(quant_fine, rounding=ROUND_HALF_EVEN)
    return float(stepped.quantize(quant, rounding=ROUND_HALF_EVEN))


def load_groups(path: str) -> dict[str, list[float]]:
    groups: dict[str, list[float]] = defaultdict(list)
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            groups[row["group"]].append(float(row["tmprss5_tpm"]))
    return groups


def main() -> None:
    groups = load_groups(CSV_PATH)

    print("=" * 64)
    print("GSE103799 TMPRSS5 Sciatic Nerve Reanalysis")
    print("=" * 64)
    print()

    print("Group summary (n, Mean TPM, SD):")
    summary = {}
    for group in GROUPS:
        values = groups[group]
        n = len(values)
        m = mean(values)
        sd = sample_sd(values)
        summary[group] = {"n": n, "mean": m, "sd": sd, "values": values}
        print(f"  {group:8s}  n={n}  Mean={m:6.2f}  SD={sd:5.2f}")

    wt = summary["WT_PBS"]["mean"]
    c22 = summary["C22_PBS"]["mean"]
    aso = summary["C22_ASO"]["mean"]

    disease_delta = c22 - wt
    aso_delta = c22 - aso
    rescue_pct = 100.0 * aso_delta / disease_delta

    disease_lfc = math.log2(c22 / wt)
    treatment_lfc = math.log2(aso / c22)

    print()
    print("Pharmacodynamic contrasts:")
    print(f"  Disease Delta (C22_PBS - WT_PBS):     {disease_delta:6.2f} TPM")
    print(f"  ASO Delta (C22_PBS - C22_ASO):         {aso_delta:6.2f} TPM")
    print(f"  Rescue Percentage:                    {rescue_pct:6.2f}%")
    print()
    print("log2FoldChange:")
    print(f"  Disease  (C22_PBS / WT_PBS):           {round_half_even(disease_lfc):+.2f}")
    print(f"  Treatment (C22_ASO / C22_PBS):         {round_half_even(treatment_lfc):+.2f}")

    t_disease = stats.ttest_ind(
        summary["C22_PBS"]["values"],
        summary["WT_PBS"]["values"],
        equal_var=False,
    )
    t_treatment = stats.ttest_ind(
        summary["C22_ASO"]["values"],
        summary["C22_PBS"]["values"],
        equal_var=False,
    )

    print()
    print("Welch's two-sample t-tests (unequal variance):")
    print(
        f"  C22_PBS vs WT_PBS:   t={t_disease.statistic:7.3f}  "
        f"p={t_disease.pvalue:.5f}"
    )
    print(
        f"  C22_ASO vs C22_PBS:  t={t_treatment.statistic:7.3f}  "
        f"p={t_treatment.pvalue:.5f}"
    )
    print()


if __name__ == "__main__":
    main()
