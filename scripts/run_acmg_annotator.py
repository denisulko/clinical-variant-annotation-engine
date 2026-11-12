"""Automated ACMG/AMP 2015 & ClinGen Bayesian Point-Scale Variant Classification Engine.

Evaluates germline clinical exome variants against ACMG/AMP evidence rules:
  - Pathogenic rules: PVS1 (+8), PS1 (+4), PM1 (+2), PM2 (+2), PP3_Mod (+2), PP3 (+1), PP5 (+1)
  - Benign rules:     BA1 (-8), BS1 (-4), BP4 (-1), BP7 (-1)

Computes the net ClinGen Bayesian point score (Tavtigian et al., Hum Mutat 2020;
Pejaver et al., Am J Hum Genet 2022) and assigns each variant to one of the five
standard clinical reporting tiers:
  Pathogenic | Likely_Pathogenic | VUS | Likely_Benign | Benign
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

NULL_CONSEQUENCES = {
    "stop_gained",
    "frameshift_variant",
    "splice_donor_variant",
    "splice_acceptor_variant",
}

RULE_POINTS: dict[str, int] = {
    "PVS1": 8,
    "PS1": 4,
    "PM1": 2,
    "PM2": 2,
    "PP3_Mod": 2,
    "PP3": 1,
    "PP5": 1,
    "BA1": -8,
    "BS1": -4,
    "BP4": -1,
    "BP7": -1,
}


def classify_single_variant(row: pd.Series | dict[str, object]) -> dict[str, object]:
    """Evaluate ACMG/AMP evidence codes and compute Bayesian point classification for one variant."""
    consequence = str(row.get("consequence", ""))
    lof_intolerant = bool(row.get("lof_intolerant_gene", False))
    in_hotspot = bool(row.get("in_hotspot_domain", False))
    gnomad_af = float(row.get("gnomad_af", 0.0))
    revel = float(row.get("revel_score", 0.0))
    cadd = float(row.get("cadd_phred", 0.0))
    spliceai = float(row.get("spliceai_ds", 0.0))
    phylop = float(row.get("phylop_score", 0.0))
    clinvar_sig = str(row.get("clinvar_assertion", "Uncertain_significance"))
    clinvar_stars = int(row.get("clinvar_stars", 0))

    triggered_rules: list[str] = []

    # --- Benign Criteria ---
    if gnomad_af > 0.05:
        triggered_rules.append("BA1")
    elif gnomad_af >= 0.0008:
        triggered_rules.append("BS1")

    if revel <= 0.18 and cadd <= 15.0 and spliceai <= 0.10:
        triggered_rules.append("BP4")

    if consequence in {"synonymous_variant", "intron_variant"} and spliceai < 0.10 and phylop < 2.0:
        triggered_rules.append("BP7")

    # --- Pathogenic Criteria (only evaluated when stand-alone BA1 is not present) ---
    if "BA1" not in triggered_rules:
        if consequence in NULL_CONSEQUENCES and lof_intolerant and gnomad_af <= 0.001:
            triggered_rules.append("PVS1")

        if clinvar_sig == "Pathogenic" and clinvar_stars >= 2 and gnomad_af <= 0.001:
            triggered_rules.append("PS1")
        elif clinvar_sig in {"Pathogenic", "Likely_pathogenic"} and clinvar_stars >= 1 and gnomad_af <= 0.001:
            triggered_rules.append("PP5")

        if in_hotspot and gnomad_af <= 0.0005:
            triggered_rules.append("PM1")

        if gnomad_af <= 0.0001:
            triggered_rules.append("PM2")

        # ClinGen SVI calibrated in silico thresholds
        if revel >= 0.773 or (revel >= 0.70 and cadd >= 28.0) or spliceai >= 0.80:
            triggered_rules.append("PP3_Mod")
        elif revel >= 0.644 or cadd >= 24.5 or spliceai >= 0.50:
            triggered_rules.append("PP3")

    bayesian_points = sum(RULE_POINTS[r] for r in triggered_rules)

    # --- 5-Tier ACMG/ClinGen Point Scale Assignment ---
    if "BA1" in triggered_rules or bayesian_points <= -6:
        acmg_tier = "Benign"
    elif bayesian_points <= -1:
        acmg_tier = "Likely_Benign"
    elif bayesian_points >= 10:
        acmg_tier = "Pathogenic"
    elif bayesian_points >= 6:
        acmg_tier = "Likely_Pathogenic"
    else:
        acmg_tier = "VUS"

    # Compute posterior probability approximation under ClinGen Bayesian framework
    # Prior P_0 = 0.10, Odds_Path = 2.0^(bayesian_points * 0.85)
    prior_p = 0.10
    odds_p = float(2.0 ** (bayesian_points * 0.85))
    posterior_prob = (odds_p * prior_p) / ((odds_p - 1.0) * prior_p + 1.0)
    if "BA1" in triggered_rules:
        posterior_prob = 0.0001

    return {
        "acmg_rules": ";".join(triggered_rules) if triggered_rules else "NONE",
        "rule_count": len(triggered_rules),
        "bayesian_points": int(bayesian_points),
        "posterior_pathogenicity": round(float(posterior_prob), 4),
        "acmg_classification": acmg_tier,
        "clinical_actionable": acmg_tier in {"Pathogenic", "Likely_Pathogenic"},
    }


def annotate_variant_cohort(df: pd.DataFrame) -> pd.DataFrame:
    """Apply ACMG/AMP classification across all rows of a clinical variant DataFrame."""
    annotations = [classify_single_variant(row) for _, row in df.iterrows()]
    ann_df = pd.DataFrame(annotations, index=df.index)
    out = pd.concat([df.copy(), ann_df], axis=1)
    # Sort by Bayesian points descending, then CADD score descending
    out = out.sort_values(["bayesian_points", "cadd_phred"], ascending=[False, False]).reset_index(drop=True)
    return out


def generate_clinical_report(annotated_df: pd.DataFrame, output_path: Path) -> None:
    """Write a structured Markdown Clinical Molecular Diagnostic Report for actionable variants."""
    actionable = annotated_df[annotated_df["clinical_actionable"]].copy()
    tier_counts = annotated_df["acmg_classification"].value_counts().to_dict()

    lines = [
        "# Clinical Molecular Genetics Diagnostic Summary Report",
        "",
        "- **Assay Panel:** ACMG Secondary Findings (SF v3.2) Hereditary Cancer & Cardiomyopathy Panel (12 Genes)",
        "- **Reference Genome:** GRCh38.p14",
        f"- **Total Variants Evaluated:** {len(annotated_df)}",
        f"- **Clinically Actionable Findings (`Pathogenic` / `Likely_Pathogenic`):** {len(actionable)}",
        "",
        "## 1. Cohort Classification Tier Breakdown",
        "",
        "| ACMG/AMP Classification Tier | Variant Count | Percentage (%) | Clinical Reporting Action |",
        "| :--- | :---: | :---: | :--- |",
    ]
    tier_order = [
        ("Pathogenic", "Reportable — Tier 1 Clinical Actionable"),
        ("Likely_Pathogenic", "Reportable — Tier 1 Clinical Actionable"),
        ("VUS", "Monitor — Periodic Re-evaluation Recommended"),
        ("Likely_Benign", "Filtered — Not Routinely Reported"),
        ("Benign", "Filtered — Stand-alone Polymorphism"),
    ]
    for tier, action_desc in tier_order:
        cnt = int(tier_counts.get(tier, 0))
        pct = 100.0 * cnt / max(len(annotated_df), 1)
        lines.append(f"| `{tier}` | {cnt} | {pct:.1f}% | {action_desc} |")

    lines.extend(
        [
            "",
            "## 2. Prioritized Actionable Variants (`Pathogenic` & `Likely_Pathogenic`)",
            "",
            "| Variant ID | Gene | HGVS Coding / Protein | Consequence | gnomAD AF | REVEL | CADD | ACMG Evidence Codes | Points | Classification |",
            "| :--- | :--- | :--- | :--- | :---: | :---: | :---: | :--- | :---: | :--- |",
        ]
    )
    for _, r in actionable.head(20).iterrows():
        lines.append(
            f"| `{r['variant_id']}` | **{r['gene']}** | `{r['hgvs_c']}` ({r['hgvs_p']}) | `{r['consequence']}` | "
            f"{r['gnomad_af']:.6f} | {r['revel_score']:.2f} | {r['cadd_phred']:.1f} | "
            f"`{r['acmg_rules']}` | **{r['bayesian_points']:+d}** | **`{r['acmg_classification']}`** |"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "data/raw/clinical_variants_input.csv")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results")
    args = parser.parse_args(argv)

    df = pd.read_csv(args.input)
    annotated = annotate_variant_cohort(df)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    annotated.to_csv(args.output_dir / "annotated_variants.csv", index=False)

    actionable = annotated[annotated["clinical_actionable"]].copy()
    actionable.to_csv(args.output_dir / "actionable_clinical_variants.csv", index=False)

    generate_clinical_report(annotated, args.output_dir / "clinical_case_report.md")

    # Compute individual rule trigger frequencies
    rule_freq: dict[str, int] = {r: 0 for r in RULE_POINTS}
    for rules_str in annotated["acmg_rules"]:
        if rules_str and rules_str != "NONE":
            for code in str(rules_str).split(";"):
                if code in rule_freq:
                    rule_freq[code] = rule_freq.get(code, 0) + 1

    tier_counts = {k: int(v) for k, v in annotated["acmg_classification"].value_counts().items()}
    summary = {
        "total_variants": int(len(annotated)),
        "actionable_variants": int(len(actionable)),
        "tier_counts": tier_counts,
        "rule_trigger_counts": rule_freq,
        "mean_cadd_pathogenic": round(
            float(annotated[annotated["acmg_classification"] == "Pathogenic"]["cadd_phred"].mean()), 2
        ),
        "mean_revel_pathogenic": round(
            float(annotated[annotated["acmg_classification"] == "Pathogenic"]["revel_score"].mean()), 3
        ),
    }
    (args.output_dir / "acmg_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"Annotated {len(annotated)} clinical variants -> {args.output_dir}")
    print(f"Actionable (Pathogenic / Likely_Pathogenic): {len(actionable)}")
    print(f"Tier breakdown: {tier_counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
