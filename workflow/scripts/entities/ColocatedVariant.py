from json import dumps
from typing import Dict, List

from workflow.scripts.entities.VariantFrequencies import VariantFrequencies
from workflow.scripts.entities.VariantSynonyms import VariantSynonyms


class ColocatedVariant:
    """This object represents a co-located variant at a location that was provided to E! Ensembl for variant-effect-prediction."""

    def __init__(
        self,
        start: int | None = None,
        strand: str | None = None,
        allele_string: str | None = None,
        seq_region_name: str | None = None,
        phenotype_or_disease: int | None = None,
        somatic: int | None = None,
        end: int | None = None,
        id: str | None = None,
        frequencies: Dict[str, str] | None = None,
        pubmed: List[int] = None,
        var_synonyms: str = None,
        clin_sig_allele: str | None = None,
        clin_sig: str | None = None,
    ):
        self.start = start
        self.strand = strand
        self.allele_string = allele_string
        self.seq_region_name = seq_region_name
        self.phenotype_or_disease = phenotype_or_disease
        self.somatic = somatic
        self.end = end
        self.id = id
        self.clin_sig_allele = clin_sig_allele
        self.clin_sig = clin_sig

        if var_synonyms is not None:
            variant_synonyms = VariantSynonyms(**var_synonyms)
            self.var_synonyms = variant_synonyms.getFlattenedList()
        else:
            self.var_synonyms = None

        if pubmed is not None:
            self.pubmed = pubmed
        else:
            self.pubmed = None

        self.frequencies = dict()

        if frequencies != None:
            for allele, freq_results in frequencies.items():
                self.frequencies[allele] = str(
                    VariantFrequencies(**freq_results).toDict()
                )
        else:
            self.frequencies = None

    def toJSON(self):
        return dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)
