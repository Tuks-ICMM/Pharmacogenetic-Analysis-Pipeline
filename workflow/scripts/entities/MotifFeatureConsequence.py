from typing import List

from workflow.scripts.entities.MotifFeaturePhenotype import MotifFeaturePhenotype


class MotifFeatureConsequence:
    def __init__(
        self,
        cadd_raw: float | None = None,
        impact: str | None = None,
        transcription_factors: List[str] | None = None,
        motif_score_change: float | None = None,
        motif_feature_id: str | None = None,
        phenotypes: List[MotifFeaturePhenotype] | None = None,
        motif_name: str | None = None,
        cadd_phred: float | None = None,
        high_inf_pos: str | None = None,
        variant_allele: str | None = None,
        consequence_terms: List[str] | None = None,
        motif_pos: int | None = None,
        strand: int | None = None,
    ) -> None:
        self.cadd_raw = (cadd_raw,)
        self.impact = (impact,)
        self.transcription_factors = (transcription_factors,)
        self.motif_score_change = (motif_score_change,)
        self.motif_feature_id = (motif_feature_id,)
        self.phenotypes = (phenotypes,)
        self.motif_name = (motif_name,)
        self.cadd_phred = (cadd_phred,)
        self.high_inf_pos = (high_inf_pos,)
        self.variant_allele = (variant_allele,)
        self.consequence_terms = (consequence_terms,)
        self.motif_pos = (motif_pos,)
        self.strand = (strand,)
