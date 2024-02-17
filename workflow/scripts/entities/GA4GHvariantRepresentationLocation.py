from workflow.scripts.entities.GA4GHvariantRepresentationInterval import (
    GA4GHvariantRepresentationInterval,
)


class GA4GHvariantRepresentationLocation:
    def __init__(
        self,
        type: str | None = None,
        sequence_id: str | None = None,
        interval: GA4GHvariantRepresentationInterval | None = None,
    ) -> None:
        self.type = type
        self.sequence_id = sequence_id
        if interval is not None:
            self.interval = interval
