from itertools import repeat
from json import dumps
from typing import Dict, List, Union

from workflow.scripts.entities.ColocatedVariant import ColocatedVariant
from workflow.scripts.entities.MotifFeatureConsequence import MotifFeatureConsequence
from workflow.scripts.entities.RegulatoryFeatureConsequence import (
    RegulatoryFeatureConsequence,
)
from workflow.scripts.entities.TranscriptConsequence import TranscriptConsequence
from workflow.scripts.entities.TranscriptPhenotype import TranscriptPhenotype


class VariantConsequenceResult:
    """This class represents a Variant Effect Prediction consequence as returned by E! Ensembl."""

    # def __init__(
    #     self,
    #     # payload,
    #     # assembly_name: str | None = None,
    #     # id: str | None = None,
    #     # strand: int | None = None,
    #     # most_severe_consequence: str | None = None,
    #     # start: str | None = None,
    #     # allele_string: str | None = None,
    #     # end: int | None = None,
    #     # colocated_variants: List[ColocatedVariant] | None = None,
    #     # input: str | None = None,
    #     # seq_region_name: str | None = None,
    #     # motif_feature_consequences: List[MotifFeatureConsequence] | None = None,
    #     # transcript_consequences: List[TranscriptConsequence] | None = None,
    #     # regulatory_feature_consequences: List[RegulatoryFeatureConsequence]
    #     # | None = None,
    #     # **other,
    # ):

    @classmethod
    def from_dict(cls, payload):
        obj = cls()
        obj.__dict__.update(payload)
        return obj

    #     self.assembly_name = assembly_name

    #     self.id = id
    #     self.strand = strand
    #     self.most_severe_consequence = most_severe_consequence
    #     self.start = start
    #     self.allele_string = allele_string
    #     self.end = end

    #     # TODO: Convert to custom setter for 'colocated_variants' property

    #     if colocated_variants is not None:
    #         self.colocated_variants = list(
    #             map(lambda variant: ColocatedVariant(**variant), colocated_variants)
    #         )
    #     else:
    #         self.colocated_variants = None

    #     # TODO: Create custom setter for 'regulatory_feature_consequences' property
    #     if regulatory_feature_consequences is not None:
    #         mapped = map(
    #             lambda consequence: RegulatoryFeatureConsequence(**consequence),
    #             regulatory_feature_consequences,
    #         )
    #         self.regulatory_feature_consequences = list(
    #             filter(
    #                 lambda consequence: consequence.source is not None
    #                 or consequence.review_status is not None,
    #                 mapped,
    #             )
    #         )
    #     else:
    #         self.regulatory_feature_consequences = None

    #     self.input = input
    #     self.seq_region_name = seq_region_name

    #     # TODO: Create custom setter for 'motif_feature_consequences' property
    #     if motif_feature_consequences is not None:
    #         self.motif_feature_consequences = list(
    #             map(
    #                 lambda consequence: MotifFeatureConsequence(**consequence),
    #                 motif_feature_consequences,
    #             )
    #         )
    #     else:
    #         self.motif_feature_consequences = None

    #     if transcript_consequences is not None:
    #         self.transcript_consequences = list(
    #             map(
    #                 lambda consequence: TranscriptConsequence(**consequence),
    #                 transcript_consequences,
    #             )
    #         )
    #     else:
    #         self.transcript_consequences = None

    # def getBiotypes(self) -> str | None:
    #     if self.regulatory_feature_consequences is not None:
    #         return "|".join(
    #             {
    #                 consequence.biotype
    #                 for consequence in self.regulatory_feature_consequences
    #             }
    #         )
    #     else:
    #         return None

    # def getPhenotypes(self) -> str | None:
    #     if self.regulatory_feature_consequences is not None:
    #         return "|".join(
    #             {
    #                 consequence.phenotypes
    #                 for consequence in self.regulatory_feature_consequences
    #             }
    #         )
    #     else:
    #         return None

    # def getImpact(self) -> str | None:
    #     if self.regulatory_feature_consequences is not None:
    #         return "|".join(
    #             {
    #                 consequence.impact
    #                 for consequence in self.regulatory_feature_consequences
    #             }
    #         )
    #     else:
    #         return None

    # def matchTranscript(self, transcript_id: str) -> TranscriptConsequence | None:
    #     """Within a given variant consequence report, more than one transcript can match to the given variation. In this case, this function will identify and return the matching transcript based on ID."""
    #     # TODO: Make transcript fallback in case transcript not found
    #     transcript: TranscriptConsequence | None = next(
    #         (
    #             transcript
    #             for transcript in self.transcript_consequences
    #             if transcript.mane_select == transcript_id
    #         ),
    #         None,
    #     )
    #     return transcript

    # def fetchCanonicalTranscript(self) -> TranscriptConsequence | None:
    #     return next(
    #         transcript
    #         for transcript in self.transcript_consequences
    #         if transcript.canonical
    #     )

    # def matchPhenotypeOverlap(
    #     self, phenotypes: TranscriptPhenotype | RegulatoryFeatureConsequence
    # ) -> List[TranscriptPhenotype] | None:
    #     filtered = filter(
    #         lambda phenotype: self.start <= phenotype.start <= self.end
    #         or self.start <= phenotype.end <= self.end,
    #         phenotypes,
    #     )
    #     if filtered is not None:
    #         phenotype = list(filtered)
    #         return phenotype
    #     else:
    #         return filtered

    # def fetchAttributsInList(self, transcript_phenotypes, attribute):
    #     return list(
    #         filter(
    #             lambda pheno: pheno is not None,
    #             map(
    #                 lambda phenotype: getattr(phenotype, f"{attribute}"),
    #                 transcript_phenotypes,
    #             ),
    #         )
    #     )

    # def generateDataFrameRow(
    #     self,
    # ) -> List[Union[int, str, float, None]]:
    #     # First, start by building a blank list for the row that represents the values to be added into the DataFrame:
    #     row = list()
    #     if self.transcript_consequences is not None:
    #         transcript = self.fetchCanonicalTranscript()
    #     else:
    #         transcript = None

    #     if transcript is not None and transcript.phenotypes is not None:
    #         transcript_phenotypes = self.matchPhenotypeOverlap(transcript.phenotypes)
    #     else:
    #         transcript_phenotypes = None

    #     row.extend(
    #         [
    #             self.id,
    #             self.start,
    #             self.end,
    #             self.allele_string,
    #             self.strand,
    #             self.assembly_name,
    #             self.most_severe_consequence,
    #             self.seq_region_name,
    #             # Matching transcript consequence fields
    #         ]
    #     )
    #     # Colocated variants
    #     if self.colocated_variants is not None:
    #         id = self.fetchAttributsInList(self.colocated_variants, "id")
    #         if len(id) != 0:
    #             row.append("|".join(id))
    #         else:
    #             row.append(None)

    #         start = self.fetchAttributsInList(self.colocated_variants, "start")
    #         if len(start) != 0:
    #             row.append("|".join([str(elem) for elem in start]))
    #         else:
    #             row.append(None)

    #         end = self.fetchAttributsInList(self.colocated_variants, "end")
    #         if len(end) != 0:
    #             row.append("|".join([str(elem) for elem in end]))
    #         else:
    #             row.append(None)

    #         allele_string = self.fetchAttributsInList(
    #             self.colocated_variants, "allele_string"
    #         )
    #         if len(allele_string) != 0:
    #             row.append("|".join(allele_string))
    #         else:
    #             row.append(None)

    #         frequencies = self.fetchAttributsInList(
    #             self.colocated_variants, "frequencies"
    #         )
    #         if len(frequencies) != 0:
    #             row.append(dumps(frequencies))
    #         else:
    #             row.append(None)

    #         var_synonyms = self.fetchAttributsInList(
    #             self.colocated_variants, "var_synonyms"
    #         )
    #         if len(var_synonyms) != 0:
    #             result = "|".join(var_synonyms[0])
    #             row.append(result)
    #         else:
    #             row.append(None)

    #         pubmed = self.fetchAttributsInList(self.colocated_variants, "pubmed")
    #         if len(pubmed) != 0:
    #             print("PUBMED: ", [str(pubmedId) for pubmedId in pubmed[0]])
    #             row.append("|".join([str(pubmedId) for pubmedId in pubmed[0]]))
    #         else:
    #             row.append(None)
    #     else:
    #         row.extend(list(repeat(None, 7)))

    #     if transcript is not None:
    #         row.extend(
    #             [
    #                 transcript.id,
    #                 transcript.given_ref,
    #                 transcript.used_ref,
    #                 transcript.amino_acids,
    #                 transcript.biotype,
    #                 transcript.strand,
    #                 transcript.variant_allele,
    #                 transcript.gene_symbol,
    #                 transcript.hgvsc,
    #                 transcript.hgvsp,
    #                 transcript.cadd_phred,
    #                 transcript.cadd_raw,
    #                 transcript.sift_score,
    #                 transcript.sift_prediction,
    #                 transcript.polyphen_score,
    #                 transcript.polyphen_prediction,
    #             ]
    #         )
    #     else:
    #         row.extend(list(repeat(None, 16)))

    #     if transcript_phenotypes is not None:
    #         # Matching fields for transcript phenotypes with overlap
    #         transcript_phenotypes_phenotype = self.fetchAttributsInList(
    #             transcript_phenotypes, "phenotype"
    #         )
    #         if (
    #             len(transcript_phenotypes_phenotype) != 0
    #             and transcript_phenotypes_phenotype[0] is not None
    #         ):
    #             row.append("|".join(transcript_phenotypes_phenotype))
    #         else:
    #             row.append(None)

    #         source = self.fetchAttributsInList(transcript_phenotypes, "source")
    #         if len(source) != 0 and source[0] is not None:
    #             row.append("|".join(source))
    #         else:
    #             row.append(None)

    #         review_status = self.fetchAttributsInList(
    #             transcript_phenotypes, "review_status"
    #         )
    #         if len(review_status) != 0:
    #             row.append("|".join(review_status))
    #         else:
    #             row.append(None)

    #         datelastevaluated = self.fetchAttributsInList(
    #             transcript_phenotypes, "datelastevaluated"
    #         )
    #         if len(datelastevaluated) != 0 and datelastevaluated[0] is not None:
    #             row.append("|".join(datelastevaluated))
    #         else:
    #             row.append(None)

    #         submitter_names = self.fetchAttributsInList(
    #             transcript_phenotypes, "submitter_name"
    #         )
    #         if len(submitter_names) != 0 and submitter_names[0] is not None:
    #             row.append("|".join(submitter_names))
    #         else:
    #             row.append(None)
    #     else:
    #         row.append(None)
    #         row.append(None)
    #         row.append(None)
    #         row.append(None)
    #         row.append(None)

    #     # Matching fields for regulatory feature with overlap

    #     if self.regulatory_feature_consequences is not None:
    #         row.append(
    #             "|".join(
    #                 {
    #                     regfeature.id
    #                     for regfeature in self.regulatory_feature_consequences
    #                 }
    #             )
    #         )
    #         row.append(
    #             "|".join(
    #                 {
    #                     regfeature.biotype
    #                     for regfeature in self.regulatory_feature_consequences
    #                 }
    #             )
    #         )
    #         row.append(
    #             "|".join(
    #                 {
    #                     regfeature.impact
    #                     for regfeature in self.regulatory_feature_consequences
    #                 }
    #             )
    #         )
    #         row.append(
    #             "|".join(
    #                 [
    #                     str(regfeature.cadd_phred)
    #                     for regfeature in self.regulatory_feature_consequences
    #                 ]
    #             )
    #         )
    #         row.append(
    #             "|".join(
    #                 [
    #                     str(regfeature.cadd_raw)
    #                     for regfeature in self.regulatory_feature_consequences
    #                 ]
    #             )
    #         )
    #     else:
    #         row.append(None)
    #         row.append(None)
    #         row.append(None)
    #         row.append(None)
    #         row.append(None)

    #     if self.regulatory_feature_consequences is not None:
    #         row.append(
    #             "|".join(
    #                 {
    #                     phenotype.source
    #                     for phenotype in self.regulatory_feature_consequences
    #                 }
    #             )
    #         )
    #         row.append(
    #             "|".join(
    #                 [
    #                     phenotype.review_status
    #                     for phenotype in self.regulatory_feature_consequences
    #                 ]
    #             )
    #         )
    #         row.append(
    #             "|".join(
    #                 [
    #                     phenotype.datelastevaluated
    #                     for phenotype in self.regulatory_feature_consequences
    #                 ]
    #             )
    #         )
    #         row.append(
    #             "|".join(
    #                 [
    #                     phenotype.submitter_name
    #                     for phenotype in self.regulatory_feature_consequences
    #                 ]
    #             )
    #         )
    #     else:
    #         row.append(None)
    #         row.append(None)
    #         row.append(None)
    #         row.append(None)

    #     return row

    # def colocated_pubmed_ids(self) -> List[str] | None:
    #     # Here we use list-comprehension, which is an excellent way to comb through a list of complex objects,
    #     # and pull out specific values in a list. This one is more complicated since we have to nest this, and go two layers down.
    #     if self.colocated_variants:
    #         ID_LISTS = [
    #             colocated_variant.pubmed
    #             for colocated_variant in self.colocated_variants
    #             if colocated_variant.pubmed
    #         ]

    #         LIST_OF_IDS = [
    #             str(pubmedId)
    #             for VARIANT_CONSEQUENCE in ID_LISTS
    #             for pubmedId in VARIANT_CONSEQUENCE
    #         ]

    #         if LIST_OF_IDS:
    #             return ",".join(LIST_OF_IDS)
    #         else:
    #             return None
    #     else:
    #         return str()
