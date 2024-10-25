from typing import List

from workflow.scripts.entities.GA4GHvariantRepresentationSpec import (
    GA4GHvariantRepresentationSpec,
)

from workflow.scripts.entities.TranscriptPhenotype import TranscriptPhenotype


class RegulatoryFeatureConsequence:
    def __init__(
        self,
        phenotypes: TranscriptPhenotype | None = None,
        regulatory_feature_id: str | None = None,
        biotype: str | None = None,
        cadd_phred: int | None = None,
        cadd_raw: int | None = None,
        impact: str | None = None,
        consequence_terms: List[str] = None,
        variant_allele: str | None = None,
        source: str | None = None,
        review_status: str | None = None,
        datelastevaluated: str | None = None,
        submitter_name: str | None = None,
        ga4gh_vrs: GA4GHvariantRepresentationSpec | None = None,
        aa: str | None = None,
        **other
    ):
        self.id = regulatory_feature_id
        self.biotype = biotype
        self.cadd_phred = cadd_phred
        self.cadd_raw = cadd_raw
        self.impact = impact
        self.consequence_terms = consequence_terms
        self.variant_allele = variant_allele
        self.source = source
        self.review_status = review_status
        self.datelastevaluated = datelastevaluated
        self.submitter_name = submitter_name
        self.aa = aa

        # TODO: Add custom processor here
        if phenotypes is not None:
            self.phenotypes = map(
                lambda phenotype: TranscriptPhenotype(**phenotype), phenotypes
            )

        # TODO: Add custom processor here
        if ga4gh_vrs is not None:
            self.ga4gh_vrs = GA4GHvariantRepresentationSpec(**ga4gh_vrs)
