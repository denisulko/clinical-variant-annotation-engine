"""Generate publication-quality 300 DPI clinical genomics & ACMG visualization figures.

Produces three multi-panel diagnostic figures in results/:
  1. acmg_classification_breakdown.png — 5-Tier ACMG classification distribution & rule trigger frequency
  2. insilico_revel_cadd_landscape.png — REVEL vs. CADD Phred & gnomAD AF vs. ClinGen Bayesian Score
  3. consequence_gene_burden.png       — Molecular consequence composition & 12-gene panel burden
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

TIER_ORDER = ["Pathogenic", "Likely_Pathogenic", "VUS", "Likely_Benign", "Benign"]
TIER_LABELS = {
    "Pathogenic": "Pathogenic (P)",
    "Likely_Pathogenic": "Likely Pathogenic (LP)",
    "VUS": "Uncertain Significance (VUS)",
    "Likely_Benign": "Likely Benign (LB)",
    "Benign": "Benign (B)",
}
TIER_COLORS = {
    "Pathogenic": "#dc2626",
    "Likely_Pathogenic": "#f97316",
    "VUS": "#eab308",
    "Likely_Benign": "#38bdf8",
    "Benign": "#2563eb",
}


def plot_acmg_breakdown(df: pd.DataFrame, summary: dict, output_path: Path) -> None:
    """Plot 2-panel 5-tier ACMG classification counts and evidence code trigger frequencies."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.5, 6.2), dpi=300)

    # Panel A: 5-Tier Classification Horizontal Bar Chart
    counts = [int((df["acmg_classification"] == t).sum()) for t in TIER_ORDER]
    y_pos = np.arange(len(TIER_ORDER))
    colors = [TIER_COLORS[t] for t in TIER_ORDER]
    labels = [TIER_LABELS[t] for t in TIER_ORDER]

    bars = ax1.barh(y_pos, counts, color=colors, edgecolor="#1e293b", height=0.62, alpha=0.92)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(labels, fontsize=10.5, fontweight="bold")
    ax1.invert_yaxis()
    ax1.set_xlabel("Number of Germline Variants", fontsize=11, fontweight="bold")
    ax1.set_title("A. 5-Tier ACMG/AMP Clinical Classification Breakdown", fontsize=12.5, fontweight="bold", pad=12)
    ax1.set_xlim(0, max(counts) * 1.24)
    ax1.grid(axis="x", linestyle="--", alpha=0.35)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    total = len(df)
    for bar, cnt in zip(bars, counts):
        pct = 100.0 * cnt / max(total, 1)
        ax1.text(
            bar.get_width() + max(counts) * 0.02,
            bar.get_y() + bar.get_height() / 2,
            f"{cnt} ({pct:.1f}%)",
            va="center",
            ha="left",
            fontsize=10,
            fontweight="bold",
            color="#0f172a",
        )

    # Panel B: ACMG Evidence Rule Trigger Frequencies
    rule_counts: dict[str, int] = summary.get("rule_trigger_counts", {})
    rule_order = ["PVS1", "PS1", "PM1", "PM2", "PP3_Mod", "PP3", "PP5", "BA1", "BS1", "BP4", "BP7"]
    r_vals = [int(rule_counts.get(r, 0)) for r in rule_order]
    r_colors = [
        "#b91c1c" if r in {"PVS1", "PS1"}
        else "#ea580c" if r in {"PM1", "PM2", "PP3_Mod"}
        else "#f59e0b" if r in {"PP3", "PP5"}
        else "#1d4ed8" if r in {"BA1", "BS1"}
        else "#0284c7"
        for r in rule_order
    ]

    x_pos = np.arange(len(rule_order))
    r_bars = ax2.bar(x_pos, r_vals, color=r_colors, edgecolor="#1e293b", width=0.65, alpha=0.92)
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(rule_order, rotation=30, ha="right", fontsize=10, fontweight="bold")
    ax2.set_ylabel("Variants Triggering Criterion", fontsize=11, fontweight="bold")
    ax2.set_title("B. ACMG/ClinGen Evidence Code Activation Profile", fontsize=12.5, fontweight="bold", pad=12)
    ax2.set_ylim(0, max(r_vals) * 1.18)
    ax2.grid(axis="y", linestyle="--", alpha=0.35)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    for bar, val in zip(r_bars, r_vals):
        ax2.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(r_vals) * 0.02,
            str(val),
            ha="center",
            va="bottom",
            fontsize=9.5,
            fontweight="bold",
            color="#0f172a",
        )

    fig.tight_layout(pad=2.0)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_insilico_landscape(df: pd.DataFrame, output_path: Path) -> None:
    """Plot 2-panel REVEL vs. CADD scatter and gnomAD AF vs. ClinGen Bayesian Points."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.5, 6.4), dpi=300)

    # Panel A: REVEL vs. CADD Phred Scatter
    for tier in reversed(TIER_ORDER):
        sub = df[df["acmg_classification"] == tier]
        if sub.empty:
            continue
        ax1.scatter(
            sub["revel_score"],
            sub["cadd_phred"],
            label=f"{TIER_LABELS[tier]} (n={len(sub)})",
            color=TIER_COLORS[tier],
            s=58 if tier in {"Pathogenic", "Likely_Pathogenic"} else 42,
            alpha=0.88,
            edgecolors="#0f172a",
            linewidths=0.5,
        )

    # ClinGen SVI calibration thresholds
    ax1.axvline(0.773, color="#dc2626", linestyle="--", linewidth=1.2, alpha=0.8, label="PP3_Mod REVEL >= 0.773")
    ax1.axvline(0.180, color="#2563eb", linestyle=":", linewidth=1.2, alpha=0.8, label="BP4 REVEL <= 0.180")
    ax1.axhline(28.0, color="#ea580c", linestyle="--", linewidth=1.1, alpha=0.7)

    # Annotate 4 distinct benchmark pathogenic variants without label overlap
    benchmark_targets = [
        ("CLINVAR_006", (-95, 10)),   # LDLR c.682G>T (0.95, 41.0)
        ("CLINVAR_005", (-128, 2)),   # MYBPC3 c.1504C>T (0.93, 38.2)
        ("CLINVAR_001", (-132, -18)), # BRCA1 c.5266dupC (0.94, 36.4)
        ("CLINVAR_003", (-112, -36)), # TP53 c.743G>A (0.97, 34.8)
    ]
    for vid, (dx, dy) in benchmark_targets:
        match = df[df["variant_id"] == vid]
        if match.empty:
            continue
        r = match.iloc[0]
        ax1.annotate(
            f"{r['gene']} {r['hgvs_c']}",
            xy=(float(r["revel_score"]), float(r["cadd_phred"])),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=8.2,
            fontweight="bold",
            color="#7f1d1d",
            bbox=dict(boxstyle="round,pad=0.2", fc="#fef2f2", ec="#f87171", lw=0.7),
            arrowprops=dict(arrowstyle="->", color="#991b1b", lw=0.7),
        )

    ax1.set_xlabel("REVEL Ensemble Pathogenicity Score (0.0 - 1.0)", fontsize=11, fontweight="bold")
    ax1.set_ylabel("CADD Phred Score (-10 log10 Rank)", fontsize=11, fontweight="bold")
    ax1.set_title("A. In Silico Deleteriousness Landscape (REVEL vs. CADD)", fontsize=12.5, fontweight="bold", pad=12)
    ax1.set_xlim(-0.04, 1.05)
    ax1.set_ylim(-1.0, 45.5)
    ax1.grid(True, linestyle="--", alpha=0.3)
    ax1.legend(loc="upper left", fontsize=8.5, frameon=True, facecolor="white", framealpha=0.92)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Panel B: log10(gnomAD AF + 1e-6) vs. ClinGen Bayesian Points
    log_af = np.log10(df["gnomad_af"].clip(lower=1e-6))
    rng = np.random.default_rng(42)
    jitter_x = log_af + rng.normal(0, 0.06, size=len(df))
    jitter_y = df["bayesian_points"] + rng.normal(0, 0.22, size=len(df))

    for tier in reversed(TIER_ORDER):
        mask = df["acmg_classification"] == tier
        if not mask.any():
            continue
        ax2.scatter(
            jitter_x[mask],
            jitter_y[mask],
            label=TIER_LABELS[tier],
            color=TIER_COLORS[tier],
            s=52,
            alpha=0.88,
            edgecolors="#0f172a",
            linewidths=0.5,
        )

    # AF & Point Thresholds
    ax2.axvline(np.log10(0.0001), color="#ea580c", linestyle="--", lw=1.1, label="PM2 Cutoff (AF <= 1e-4)")
    ax2.axvline(np.log10(0.05), color="#1d4ed8", linestyle="-.", lw=1.2, label="BA1 Stand-Alone (AF > 5%)")
    ax2.axhline(10, color="#dc2626", linestyle=":", lw=1.1, alpha=0.75)
    ax2.axhline(6, color="#f97316", linestyle=":", lw=1.1, alpha=0.75)
    ax2.axhline(-1, color="#38bdf8", linestyle=":", lw=1.1, alpha=0.75)
    ax2.axhline(-6, color="#2563eb", linestyle=":", lw=1.1, alpha=0.75)

    ax2.set_xlabel("gnomAD v4.1 Allele Frequency (log10 scale, Absent = -6)", fontsize=11, fontweight="bold")
    ax2.set_ylabel("ClinGen Net Bayesian Point Score", fontsize=11, fontweight="bold")
    ax2.set_title("B. Population Allele Frequency vs. Bayesian Pathogenicity Score", fontsize=12.5, fontweight="bold", pad=12)
    ax2.grid(True, linestyle="--", alpha=0.3)
    ax2.legend(loc="upper right", fontsize=8.5, frameon=True, facecolor="white", framealpha=0.92)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    fig.tight_layout(pad=2.0)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_consequence_and_gene_burden(df: pd.DataFrame, output_path: Path) -> None:
    """Plot 2-panel Molecular Consequence stacked bar chart and 12-Gene Panel clinical burden."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.8, 6.4), dpi=300)

    # Panel A: Consequence Breakdown by ACMG Tier
    csq_counts = df.groupby(["consequence", "acmg_classification"]).size().unstack(fill_value=0)
    csq_order = csq_counts.sum(axis=1).sort_values(ascending=True).index.tolist()
    csq_counts = csq_counts.loc[csq_order]

    y_pos = np.arange(len(csq_order))
    left = np.zeros(len(csq_order))
    for tier in TIER_ORDER:
        vals = csq_counts[tier].values if tier in csq_counts.columns else np.zeros(len(csq_order))
        ax1.barh(
            y_pos,
            vals,
            left=left,
            color=TIER_COLORS[tier],
            edgecolor="#1e293b",
            height=0.62,
            label=TIER_LABELS[tier],
            alpha=0.92,
        )
        left += vals

    clean_csq = [c.replace("_variant", "").replace("_", " ").title() for c in csq_order]
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(clean_csq, fontsize=10, fontweight="bold")
    ax1.set_xlabel("Number of Variants", fontsize=11, fontweight="bold")
    ax1.set_title("A. Sequence Ontology Consequence by ACMG Tier", fontsize=12.5, fontweight="bold", pad=12)
    ax1.grid(axis="x", linestyle="--", alpha=0.35)
    ax1.legend(loc="lower right", fontsize=8.5, frameon=True, facecolor="white")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Panel B: 12-Gene Hereditary Cancer & Cardiomyopathy Burden
    gene_counts = df.groupby(["gene", "acmg_classification"]).size().unstack(fill_value=0)
    # Sort genes by actionable count + total count
    actionable_Load = (
        gene_counts.get("Pathogenic", 0) + gene_counts.get("Likely_Pathogenic", 0)
    )
    gene_order = (
        pd.DataFrame({"act": actionable_Load, "tot": gene_counts.sum(axis=1)})
        .sort_values(["act", "tot"], ascending=[False, False])
        .index.tolist()
    )
    gene_counts = gene_counts.loc[gene_order]

    x_pos = np.arange(len(gene_order))
    bottom = np.zeros(len(gene_order))
    for tier in TIER_ORDER:
        vals = gene_counts[tier].values if tier in gene_counts.columns else np.zeros(len(gene_order))
        ax2.bar(
            x_pos,
            vals,
            bottom=bottom,
            color=TIER_COLORS[tier],
            edgecolor="#1e293b",
            width=0.62,
            label=TIER_LABELS[tier],
            alpha=0.92,
        )
        bottom += vals

    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(gene_order, rotation=32, ha="right", fontsize=10, fontweight="bold")
    ax2.set_ylabel("Variant Count per Gene", fontsize=11, fontweight="bold")
    ax2.set_title("B. Variant Spectrum Across 12 ACMG SF v3.2 Panel Genes", fontsize=12.5, fontweight="bold", pad=12)
    ax2.set_ylim(0, max(bottom) * 1.16)
    ax2.grid(axis="y", linestyle="--", alpha=0.35)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    # Annotate actionable count above each gene bar
    for idx, g in enumerate(gene_order):
        act_n = int(gene_counts.loc[g, ["Pathogenic", "Likely_Pathogenic"]].sum()) if "Pathogenic" in gene_counts.columns else 0
        if act_n > 0:
            ax2.text(
                idx,
                bottom[idx] + max(bottom) * 0.02,
                f"{act_n} P/LP",
                ha="center",
                va="bottom",
                fontsize=8.2,
                fontweight="bold",
                color="#b91c1c",
            )

    fig.tight_layout(pad=2.0)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotated-csv", type=Path, default=ROOT / "results/annotated_variants.csv")
    parser.add_argument("--summary-json", type=Path, default=ROOT / "results/acmg_summary.json")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results")
    args = parser.parse_args(argv)

    df = pd.read_csv(args.annotated_csv)
    summary = json.loads(args.summary_json.read_text(encoding="utf-8"))

    p1 = args.output_dir / "acmg_classification_breakdown.png"
    p2 = args.output_dir / "insilico_revel_cadd_landscape.png"
    p3 = args.output_dir / "consequence_gene_burden.png"

    plot_acmg_breakdown(df, summary, p1)
    plot_insilico_landscape(df, p2)
    plot_consequence_and_gene_burden(df, p3)

    print(f"Saved 300 DPI plots:\n  - {p1}\n  - {p2}\n  - {p3}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
