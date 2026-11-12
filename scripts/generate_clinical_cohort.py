"""Generate a realistic Hereditary Cancer & Cardiomyopathy Clinical Exome VCF and Annotation DB.

Simulates 180 germline variants across 12 ACMG Secondary Findings (SF v3.2) actionable
genes (BRCA1, BRCA2, TP53, PALB2, MLH1, MSH2, APC, ATM, MYBPC3, MYH7, LDLR, SCN5A)
with gnomAD v4.1 population allele frequencies, ClinVar significance & review stars,
functional domain coordinates, and in silico pathogenicity scores (REVEL, CADD, SpliceAI).
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

# Gene panel: (symbol, chrom, start_pos, lof_intolerant_hi, syndrome, critical_domain)
GENE_PANEL: list[tuple[str, str, int, bool, str, str]] = [
    ("BRCA1", "chr17", 43044295, True, "Hereditary Breast & Ovarian Cancer", "BRCT_RING_Domain"),
    ("BRCA2", "chr13", 32315474, True, "Hereditary Breast & Ovarian Cancer", "RAD51_Binding_BBRC"),
    ("TP53", "chr17", 7668402, True, "Li-Fraumeni Syndrome", "p53_DNA_Binding_Domain"),
    ("PALB2", "chr16", 23603160, True, "Hereditary Breast & Pancreatic Cancer", "WD40_Repeat_Domain"),
    ("MLH1", "chr3", 36993350, True, "Lynch Syndrome", "ATPase_Mismatch_Repair"),
    ("MSH2", "chr2", 47403067, True, "Lynch Syndrome", "DNA_Mismatch_Binding"),
    ("APC", "chr5", 112707498, True, "Familial Adenomatous Polyposis", "Beta_Catenin_Binding"),
    ("ATM", "chr11", 108222484, True, "Ataxia-Telangiectasia / HBOC", "PI3K_PI4K_Kinase"),
    ("MYBPC3", "chr11", 47331406, True, "Hypertrophic Cardiomyopathy", "Myosin_Binding_C0_C10"),
    ("MYH7", "chr14", 23412740, False, "Hypertrophic Cardiomyopathy", "Myosin_Head_Motor_Domain"),
    ("LDLR", "chr19", 11089362, True, "Familial Hypercholesterolemia", "EGF_Precursor_LDLR_ClassB"),
    ("SCN5A", "chr3", 38548062, False, "Brugada / Long QT Syndrome 3", "Voltage_Sensor_Transmembrane"),
]

CURATED_PATHOGENIC_BENCHMARKS: list[dict[str, object]] = [
    {
        "gene": "BRCA1", "chrom": "chr17", "pos": 43091983, "ref": "G", "alt": "A",
        "hgvs_c": "c.5266dupC", "hgvs_p": "p.Gln1756ProfsTer74",
        "consequence": "frameshift_variant", "gnomad_af": 0.000018,
        "revel": 0.94, "cadd": 36.4, "spliceai": 0.04, "phylop": 7.8,
        "in_hotspot_domain": True, "clinvar_sig": "Pathogenic", "clinvar_stars": 4,
        "zygo": "0/1", "dp": 142, "vaf": 0.49,
    },
    {
        "gene": "BRCA2", "chrom": "chr13", "pos": 32340301, "ref": "GT", "alt": "G",
        "hgvs_c": "c.5946delT", "hgvs_p": "p.Ser1982ArgfsTer22",
        "consequence": "frameshift_variant", "gnomad_af": 0.000024,
        "revel": 0.91, "cadd": 35.1, "spliceai": 0.02, "phylop": 6.9,
        "in_hotspot_domain": True, "clinvar_sig": "Pathogenic", "clinvar_stars": 4,
        "zygo": "0/1", "dp": 128, "vaf": 0.48,
    },
    {
        "gene": "TP53", "chrom": "chr17", "pos": 7674220, "ref": "C", "alt": "T",
        "hgvs_c": "c.743G>A", "hgvs_p": "p.Arg248Gln",
        "consequence": "missense_variant", "gnomad_af": 0.0,
        "revel": 0.97, "cadd": 34.8, "spliceai": 0.01, "phylop": 9.1,
        "in_hotspot_domain": True, "clinvar_sig": "Pathogenic", "clinvar_stars": 3,
        "zygo": "0/1", "dp": 165, "vaf": 0.51,
    },
    {
        "gene": "MLH1", "chrom": "chr3", "pos": 37012077, "ref": "G", "alt": "A",
        "hgvs_c": "c.350C>T", "hgvs_p": "p.Thr117Met",
        "consequence": "splice_donor_variant", "gnomad_af": 0.0,
        "revel": 0.89, "cadd": 33.5, "spliceai": 0.92, "phylop": 8.4,
        "in_hotspot_domain": True, "clinvar_sig": "Pathogenic", "clinvar_stars": 3,
        "zygo": "0/1", "dp": 118, "vaf": 0.47,
    },
    {
        "gene": "MYBPC3", "chrom": "chr11", "pos": 47359280, "ref": "G", "alt": "A",
        "hgvs_c": "c.1504C>T", "hgvs_p": "p.Arg502Trp",
        "consequence": "stop_gained", "gnomad_af": 0.000012,
        "revel": 0.93, "cadd": 38.2, "spliceai": 0.05, "phylop": 7.5,
        "in_hotspot_domain": True, "clinvar_sig": "Pathogenic", "clinvar_stars": 3,
        "zygo": "0/1", "dp": 134, "vaf": 0.50,
    },
    {
        "gene": "LDLR", "chrom": "chr19", "pos": 11105588, "ref": "G", "alt": "A",
        "hgvs_c": "c.682G>T", "hgvs_p": "p.Glu228Ter",
        "consequence": "stop_gained", "gnomad_af": 0.0,
        "revel": 0.95, "cadd": 41.0, "spliceai": 0.03, "phylop": 8.8,
        "in_hotspot_domain": True, "clinvar_sig": "Pathogenic", "clinvar_stars": 4,
        "zygo": "0/1", "dp": 152, "vaf": 0.52,
    },
]


def generate_clinical_dataset(
    output_dir: Path,
    total_variants: int = 180,
    seed: int = 42,
) -> pd.DataFrame:
    """Create VCF file, gene panel manifest, and variant annotation database."""
    rng = random.Random(seed)
    gene_lookup = {g[0]: g for g in GENE_PANEL}

    records: list[dict[str, object]] = []
    for idx, bm in enumerate(CURATED_PATHOGENIC_BENCHMARKS, start=1):
        g_info = gene_lookup[str(bm["gene"])]
        records.append(
            {
                "variant_id": f"CLINVAR_{idx:03d}",
                "chrom": bm["chrom"],
                "pos": bm["pos"],
                "ref": bm["ref"],
                "alt": bm["alt"],
                "gene": bm["gene"],
                "hgvs_c": bm["hgvs_c"],
                "hgvs_p": bm["hgvs_p"],
                "consequence": bm["consequence"],
                "lof_intolerant_gene": g_info[3],
                "syndrome": g_info[4],
                "domain_name": g_info[5] if bm["in_hotspot_domain"] else "None",
                "in_hotspot_domain": bm["in_hotspot_domain"],
                "gnomad_af": bm["gnomad_af"],
                "revel_score": bm["revel"],
                "cadd_phred": bm["cadd"],
                "spliceai_ds": bm["spliceai"],
                "phylop_score": bm["phylop"],
                "clinvar_assertion": bm["clinvar_sig"],
                "clinvar_stars": bm["clinvar_stars"],
                "genotype": bm["zygo"],
                "read_depth": bm["dp"],
                "vaf": bm["vaf"],
            }
        )

    bases = ["A", "C", "G", "T"]
    while len(records) < total_variants:
        vid_num = len(records) + 1
        gene_sym, chrom, base_pos, is_hi, syndrome, dom_name = rng.choice(GENE_PANEL)
        pos = base_pos + rng.randint(500, 48_000)
        ref = rng.choice(bases)
        alt = rng.choice([b for b in bases if b != ref])

        # Simulate realistic clinical exome distribution:
        # ~9% additional Pathogenic / Likely Pathogenic, ~25% VUS, ~33% Likely Benign, ~33% Benign
        r = rng.random()
        if r < 0.08:
            # High-impact LoF or hotspot pathogenic missense
            if is_hi and rng.random() < 0.60:
                consequence = rng.choice(["stop_gained", "frameshift_variant", "splice_acceptor_variant"])
                if consequence == "frameshift_variant":
                    alt = ref + rng.choice(bases)
                revel = round(rng.uniform(0.78, 0.96), 3)
                cadd = round(rng.uniform(31.0, 42.0), 1)
                spliceai = round(rng.uniform(0.75, 0.98), 2) if "splice" in consequence else round(rng.uniform(0.01, 0.12), 2)
                in_dom = rng.random() < 0.70
                clinvar_sig = rng.choice(["Pathogenic", "Likely_pathogenic"])
                stars = rng.choice([1, 2, 3])
            else:
                consequence = "missense_variant"
                revel = round(rng.uniform(0.80, 0.97), 3)
                cadd = round(rng.uniform(27.0, 35.5), 1)
                spliceai = round(rng.uniform(0.01, 0.15), 2)
                in_dom = True
                clinvar_sig = "Likely_pathogenic"
                stars = rng.choice([2, 3])
            gnomad_af = 0.0 if rng.random() < 0.75 else round(rng.uniform(1e-6, 5e-5), 6)
            phylop = round(rng.uniform(5.2, 9.4), 2)

        elif r < 0.34:
            # Variant of Uncertain Significance (VUS) — rare missense or intronic with borderline scores
            consequence = rng.choice(["missense_variant", "missense_variant", "inframe_deletion", "intron_variant"])
            gnomad_af = 0.0 if rng.random() < 0.45 else round(rng.uniform(2e-5, 0.0004), 6)
            revel = round(rng.uniform(0.22, 0.68), 3)
            cadd = round(rng.uniform(14.0, 24.2), 1)
            spliceai = round(rng.uniform(0.02, 0.32), 2)
            phylop = round(rng.uniform(0.5, 4.8), 2)
            in_dom = rng.random() < 0.25
            clinvar_sig = "Uncertain_significance"
            stars = rng.choice([0, 1])

        elif r < 0.67:
            # Likely Benign — synonymous, intronic, or low-REVEL missense with moderate AF
            consequence = rng.choice(["synonymous_variant", "intron_variant", "missense_variant"])
            gnomad_af = round(rng.uniform(0.0008, 0.035), 5)
            revel = round(rng.uniform(0.02, 0.14), 3)
            cadd = round(rng.uniform(2.0, 13.5), 1)
            spliceai = round(rng.uniform(0.0, 0.06), 2)
            phylop = round(rng.uniform(-2.0, 1.4), 2)
            in_dom = False
            clinvar_sig = "Likely_benign"
            stars = rng.choice([1, 2])

        else:
            # Stand-alone Benign (BA1: gnomAD AF > 5%)
            consequence = rng.choice(["synonymous_variant", "intron_variant", "missense_variant"])
            gnomad_af = round(rng.uniform(0.055, 0.42), 4)
            revel = round(rng.uniform(0.01, 0.12), 3)
            cadd = round(rng.uniform(0.5, 11.0), 1)
            spliceai = round(rng.uniform(0.0, 0.04), 2)
            phylop = round(rng.uniform(-3.2, 0.8), 2)
            in_dom = False
            clinvar_sig = "Benign"
            stars = rng.choice([2, 3, 4])

        zygo = "1/1" if gnomad_af > 0.12 and rng.random() < 0.4 else "0/1"
        dp = rng.randint(85, 210)
        vaf = round(rng.uniform(0.94, 1.0), 2) if zygo == "1/1" else round(rng.uniform(0.44, 0.56), 2)

        records.append(
            {
                "variant_id": f"CLINVAR_{vid_num:03d}",
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "gene": gene_sym,
                "hgvs_c": f"c.{pos % 3500 + 100}{ref}>{alt}",
                "hgvs_p": "p.(=)" if consequence == "synonymous_variant" else f"p.Val{(pos % 900) + 1}Ala",
                "consequence": consequence,
                "lof_intolerant_gene": is_hi,
                "syndrome": syndrome,
                "domain_name": dom_name if in_dom else "None",
                "in_hotspot_domain": in_dom,
                "gnomad_af": gnomad_af,
                "revel_score": revel,
                "cadd_phred": cadd,
                "spliceai_ds": spliceai,
                "phylop_score": phylop,
                "clinvar_assertion": clinvar_sig,
                "clinvar_stars": stars,
                "genotype": zygo,
                "read_depth": dp,
                "vaf": vaf,
            }
        )

    df = pd.DataFrame(records)
    output_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_dir / "clinical_variants_input.csv", index=False)

    # Also write standard VCF v4.2 file with INFO annotations
    vcf_lines = [
        "##fileformat=VCFv4.2",
        "##reference=GRCh38.p14",
        '##INFO=<ID=GENE,Number=1,Type=String,Description="HGNC Gene Symbol">',
        '##INFO=<ID=CSQ,Number=1,Type=String,Description="Sequence Ontology Consequence">',
        '##INFO=<ID=GNOMAD_AF,Number=1,Type=Float,Description="gnomAD v4.1 Global Allele Frequency">',
        '##INFO=<ID=REVEL,Number=1,Type=Float,Description="REVEL ensemble pathogenicity score">',
        '##INFO=<ID=CADD,Number=1,Type=Float,Description="CADD Phred score">',
        '##INFO=<ID=SPLICEAI,Number=1,Type=Float,Description="SpliceAI max delta score">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
        '##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Fraction">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPROBAND_001",
    ]
    for _, r in df.iterrows():
        info_str = (
            f"GENE={r['gene']};CSQ={r['consequence']};GNOMAD_AF={r['gnomad_af']};"
            f"REVEL={r['revel_score']};CADD={r['cadd_phred']};SPLICEAI={r['spliceai_ds']}"
        )
        fmt_val = f"{r['genotype']}:{r['read_depth']}:{r['vaf']}"
        vcf_lines.append(
            f"{r['chrom']}\t{r['pos']}\t{r['variant_id']}\t{r['ref']}\t{r['alt']}\t99\tPASS\t{info_str}\tGT:DP:VAF\t{fmt_val}"
        )
    (output_dir / "clinical_panel.vcf").write_text("\n".join(vcf_lines) + "\n", encoding="utf-8")
    return df


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data/raw")
    parser.add_argument("--count", type=int, default=180)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args(argv)

    df = generate_clinical_dataset(args.output_dir, total_variants=args.count, seed=args.seed)
    print(f"Generated {len(df)} clinical panel variants across {df['gene'].nunique()} ACMG SF v3.2 genes")
    print(f"Saved CSV -> {args.output_dir / 'clinical_variants_input.csv'}")
    print(f"Saved VCF -> {args.output_dir / 'clinical_panel.vcf'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
