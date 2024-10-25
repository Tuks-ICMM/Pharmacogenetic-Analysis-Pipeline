from typing import List


class VariantSynonyms:
    """This object contains a record of all this variants synonyms on otehr databses"""

    def __init__(
        self,
        PharmGKB: List[str] | None = None,
        UniProt: List[str] | None = None,
        COSMIC: List[str] | None = None,
        ClinVar: List[str] | None = None,
        OMIM: List[float] | None = None,
    ):
        self.PharmGKB = PharmGKB
        self.UniProt = UniProt
        self.COSMIC = COSMIC
        self.ClinVar = ClinVar
        self.OMIM = OMIM

    def getFlattenedList(self):
        filtered = list(
            filter(
                lambda synonym: synonym is not None,
                [self.PharmGKB, self.UniProt, self.COSMIC, self.ClinVar, self.OMIM],
            )
        )
        return [item for row in filtered for item in row]
