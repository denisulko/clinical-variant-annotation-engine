"""Unit tests for the Automated ACMG/AMP & ClinGen Bayesian Point Scale Engine."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_acmg_annotator import annotate_variant_cohort, classify_single_variant  # noqa: E402


class TestACMGEngine(unittest.TestCase):
    """Verify deterministic ACMG/AMP rule evaluation and ClinGen Bayesian point scoring."""

    def test_pvs1_plus_pm2_triggers_pathogenic(self) -> None:
        """A rare frameshift variant in a LoF-intolerant gene with high CADD/REVEL must be Pathogenic."""
        variant = {
            "consequence": "frameshift_variant",
            "lof_intolerant_gene": True,
            "in_hotspot_domain": True,
            "gnomad_af": 0.00001,
            "revel_score": 0.92,
            "cadd_phred": 35.0,
            "spliceai_ds": 0.02,
            "phylop_score": 7.5,
            "clinvar_assertion": "Pathogenic",
            "clinvar_stars": 3,
        }
        res = classify_single_variant(variant)
        self.assertEqual(res["acmg_classification"], "Pathogenic")
        self.assertGreaterEqual(int(res["bayesian_points"]), 10)
        self.assertIn("PVS1", str(res["acmg_rules"]))
        self.assertIn("PM2", str(res["acmg_rules"]))
        self.assertTrue(bool(res["clinical_actionable"]))

    def test_ba1_standalone_benign_overrides_other_signals(self) -> None:
        """Any polymorphism with gnomAD allele frequency > 5% (BA1) must be classified as Benign."""
        variant = {
            "consequence": "missense_variant",
            "lof_intolerant_gene": True,
            "in_hotspot_domain": False,
            "gnomad_af": 0.14,
            "revel_score": 0.05,
            "cadd_phred": 4.2,
            "spliceai_ds": 0.01,
            "phylop_score": -1.2,
            "clinvar_assertion": "Benign",
            "clinvar_stars": 3,
        }
        res = classify_single_variant(variant)
        self.assertEqual(res["acmg_classification"], "Benign")
        self.assertIn("BA1", str(res["acmg_rules"]))
        self.assertLessEqual(int(res["bayesian_points"]), -8)
        self.assertFalse(bool(res["clinical_actionable"]))

    def test_hotspot_missense_triggers_likely_pathogenic(self) -> None:
        """Rare hotspot missense variant (PM1 + PM2 + PP3_Mod + PP5 = 7 pts) must be Likely_Pathogenic."""
        variant = {
            "consequence": "missense_variant",
            "lof_intolerant_gene": True,
            "in_hotspot_domain": True,
            "gnomad_af": 0.0,
            "revel_score": 0.85,
            "cadd_phred": 29.5,
            "spliceai_ds": 0.03,
            "phylop_score": 6.4,
            "clinvar_assertion": "Likely_pathogenic",
            "clinvar_stars": 2,
        }
        res = classify_single_variant(variant)
        self.assertEqual(res["acmg_classification"], "Likely_Pathogenic")
        self.assertEqual(int(res["bayesian_points"]), 7)
        self.assertTrue(bool(res["clinical_actionable"]))

    def test_borderline_missense_classified_as_vus(self) -> None:
        """Rare missense variant outside hotspot with intermediate REVEL (0.42) must remain VUS."""
        variant = {
            "consequence": "missense_variant",
            "lof_intolerant_gene": False,
            "in_hotspot_domain": False,
            "gnomad_af": 0.00005,
            "revel_score": 0.42,
            "cadd_phred": 19.0,
            "spliceai_ds": 0.05,
            "phylop_score": 2.1,
            "clinvar_assertion": "Uncertain_significance",
            "clinvar_stars": 1,
        }
        res = classify_single_variant(variant)
        self.assertEqual(res["acmg_classification"], "VUS")
        self.assertEqual(int(res["bayesian_points"]), 2)  # PM2 only
        self.assertFalse(bool(res["clinical_actionable"]))

    def test_cohort_annotation_preserves_row_count_and_sorts_by_points(self) -> None:
        """Batch cohort annotation should preserve all rows and sort highest pathogenicity first."""
        df = pd.DataFrame(
            [
                {
                    "variant_id": "V_BENIGN",
                    "consequence": "synonymous_variant",
                    "lof_intolerant_gene": True,
                    "in_hotspot_domain": False,
                    "gnomad_af": 0.12,
                    "revel_score": 0.02,
                    "cadd_phred": 1.5,
                    "spliceai_ds": 0.0,
                    "phylop_score": -0.5,
                    "clinvar_assertion": "Benign",
                    "clinvar_stars": 2,
                },
                {
                    "variant_id": "V_PATH",
                    "consequence": "stop_gained",
                    "lof_intolerant_gene": True,
                    "in_hotspot_domain": True,
                    "gnomad_af": 0.0,
                    "revel_score": 0.94,
                    "cadd_phred": 38.0,
                    "spliceai_ds": 0.01,
                    "phylop_score": 8.0,
                    "clinvar_assertion": "Pathogenic",
                    "clinvar_stars": 4,
                },
            ]
        )
        out = annotate_variant_cohort(df)
        self.assertEqual(len(out), 2)
        self.assertEqual(out.iloc[0]["variant_id"], "V_PATH")
        self.assertEqual(out.iloc[0]["acmg_classification"], "Pathogenic")
        self.assertEqual(out.iloc[1]["acmg_classification"], "Benign")


if __name__ == "__main__":
    unittest.main()
