class TranscriptConsequencePhenotype:
    def __init__(
        self,
        end: int | None = None,
        start: int | None = None,
        source: str | None = None,
        phenotype: str | None = None,
        attrib_type: str | None = None,
        external_id: int | None = None,
        id: str | None = None,
        seq_region_name: str | None = None,
        strand: int | None = None,
        type: str | None = None,
    ):
        self.end = end
        self.start = start
        self.source = source
        self.phenotype = phenotype
        self.attrib_type = attrib_type
        self.external_id = external_id
        self.id = id
        self.seq_region_name = seq_region_name
        self.strand = strand
        self.type = type
