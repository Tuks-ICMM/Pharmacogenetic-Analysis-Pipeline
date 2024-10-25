class MotifFeaturePhenotype:
    def __init__(
        self,
        phenotype: str | None = None,
        attrib_type: str | None = None,
        start: int | None = None,
        strand: str | None = None,
        id: str | None = None,
        seq_region_name: str | None = None,
        external_id: int | None = None,
        source: str | None = None,
        end: int | None = None,
        type: str | None = None,
    ) -> None:
        self.phenotype = (phenotype,)
        self.attrib_type = (attrib_type,)
        self.start = (start,)
        self.strand = (strand,)
        self.id = (id,)
        self.seq_region_name = (seq_region_name,)
        self.external_id = (external_id,)
        self.source = (source,)
        self.end = (end,)
        self.type = (type,)
