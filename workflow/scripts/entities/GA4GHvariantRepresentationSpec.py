from workflow.scripts.entities.GA4GHvariantRepresentationLocation import (
    GA4GHvariantRepresentationLocation,
)
from workflow.scripts.entities.GA4GHvariantRepresentationState import (
    GA4GHvariantRepresentationState,
)


class GA4GHvariantRepresentationSpec:
    def __init__(
        self,
        type: str | None = None,
        location: GA4GHvariantRepresentationLocation | None = None,
        state: GA4GHvariantRepresentationState | None = None,
    ) -> None:
        self.type = type

        if location is not None:
            self.location = GA4GHvariantRepresentationLocation(**location)

        if state is not None:
            self.state = GA4GHvariantRepresentationState(**state)
