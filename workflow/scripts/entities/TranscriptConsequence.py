from typing import Dict, List, Union

from workflow.scripts.entities.GA4GHvariantRepresentationSpec import (
    GA4GHvariantRepresentationSpec,
)
from workflow.scripts.entities.TranscriptConsequenceDomain import (
    TranscriptConsequenceDomain,
)
from workflow.scripts.entities.TranscriptPhenotype import TranscriptPhenotype


class TranscriptConsequence:
    def __init__(
        self,
        cds_start: str | None = None,
        given_ref: str | None = None,
        gene_symbol_source: str | None = None,
        transcript_id: str | None = None,
        cds_end: int | None = None,
        gene_symbol: str | None = None,
        hgvsp: str | None = None,
        used_ref: str | None = None,
        biotype: str | None = None,
        cdna_start: int | None = None,
        consequence_terms: List[str] | None = None,
        cadd_raw: float = 0.0,
        protein_start: int | None = None,
        strand: int | None = None,
        protein_end: int | None = None,
        gene_id: str | None = None,
        cdna_end: int | None = None,
        codons: str | None = None,
        cadd_phred: float = 0.0,
        phenotypes: List[TranscriptPhenotype] | None = None,
        impact: str | None = None,
        hgvsc: str | None = None,
        aa: str | None = None,
        variant_allele: str | None = None,
        mane_select: str | None = None,
        sift_score: float | None = None,
        sift_prediction: str | None = None,
        polyphen_score: float | None = None,
        polyphen_prediction: str | None = None,
        sift4g_score: str | None = None,
        sift4g_pred: str | None = None,
        polyphen2_hvar_score: str | None = None,
        polyphen2_hvar_pred: str | None = None,
        lof_info: str | None = None,
        lof: str | None = None,
        distance: int | None = None,
        mastermind_mmid3: str | None = None,
        protein_id: str | None = None,
        ga4gh_vrs: Dict[str, Dict[str, str]] | None = None,
        tsl: int | None = None,
        exon: str | None = None,
        amino_acids: str | None = None,
        ccds: str | None = None,
        appris: str | None = None,
        geno2mp_url: str | None = None,
        go: str | None = None,
        mastermind_counts: str | None = None,
        swissprot: List[str] | None = None,
        uniprot_isoform: List[str] | None = None,
        geno2mp_hpo_count: int | None = None,
        hgnc_id: str | None = None,
        domains: List[TranscriptConsequenceDomain] | None = None,
        spliceai: Dict[str, Union[int, str]] | None = None,
        uniparc: List[str] | None = None,
        canonical: int | None = None,
    ):
        self.cds_start = cds_start
        self.given_ref = given_ref
        self.gene_symbol_source = gene_symbol_source
        self.id = transcript_id
        self.cds_end = cds_end
        self.gene_symbol = gene_symbol
        self.hgvsp = hgvsp
        self.used_ref = used_ref
        self.biotype = biotype
        self.cdna_start = cdna_start
        self.consequence_terms = consequence_terms
        self.cadd_raw = cadd_raw
        self.protein_start = protein_start
        self.strand = strand
        self.protein_end = protein_end
        self.gene_id = gene_id
        self.cdna_end = cdna_end
        self.codons = codons
        self.cadd_phred = cadd_phred
        self.lof_info = lof_info
        self.lof = lof
        self.distance = distance

        if phenotypes is not None:
            mapped = map(lambda phenotype: TranscriptPhenotype(**phenotype), phenotypes)
            self.phenotypes = list(
                filter(
                    lambda consequence: consequence.source is not None
                    or consequence.review_status is not None,
                    mapped,
                )
            )
        else:
            self.phenotypes = None
        self.impact = impact
        self.hgvsc = hgvsc
        self.amino_acids = aa
        self.variant_allele = variant_allele
        self.mane_select = mane_select
        self.sift_score = sift_score
        self.sift_prediction = sift_prediction
        self.polyphen_score = polyphen_score
        self.polyphen_prediction = polyphen_prediction

        # Unused
        self.sift4g_pred = sift4g_pred
        self.sift4g_score = sift4g_score
        self.polyphen2_hvar_score = polyphen2_hvar_score
        self.polyphen2_hvar_pred = polyphen2_hvar_pred

        self.master_mmid3 = mastermind_mmid3
        self.protein_id = protein_id
        self.ga4gh_vrs = (ga4gh_vrs,)

        self.tsl = tsl
        self.exon = exon
        self.amino_acids = amino_acids
        self.ccds = ccds
        self.appris = appris
        self.geno2mp_url = geno2mp_url
        self.go = go
        self.mastermind_counts = mastermind_counts
        self.swissprot = swissprot
        self.uniprot_isoform = uniprot_isoform
        self.geno2mp_hpo_count = geno2mp_hpo_count
        self.hgnc_id = hgnc_id

        if domains is not None:
            self.domains = list(
                map(lambda domain: TranscriptConsequenceDomain(**domain), domains)
            )

        self.spliceai = spliceai
        self.uniparc = uniparc
        self.canonical = canonical
