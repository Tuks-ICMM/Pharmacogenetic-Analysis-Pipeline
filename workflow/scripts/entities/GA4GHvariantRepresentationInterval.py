class GA4GHvariantRepresentationInterval:
    def __init__(
        self, type: str | None = None, start: int | None = None, end: str | None = None
    ) -> None:
        self.type = type
        self.start = start
        self.end = end
