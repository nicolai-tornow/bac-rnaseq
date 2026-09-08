from __future__ import annotations
from pathlib import Path
from typing import Literal, Optional
import yaml
from pydantic import BaseModel, Field, model_validator

Strand = Literal["reverse", "forward", "unstranded"]
Species = Literal["mabs", "mtb", "custom"]


class Reference(BaseModel):
    species: Species
    fasta: Optional[str] = None
    gff: Optional[str] = None
    exclude_seqids: list[str] = Field(default_factory=list)
    strandedness: Strand = "reverse"
    feature_types: list[str] = Field(default_factory=lambda: ["gene"])
    id_attribute: str = "locus_tag"
    seqid_map: dict[str, str] = Field(default_factory=dict)

    @model_validator(mode="after")
    def _check(self):
        if self.species == "custom":
            if not self.fasta or not self.gff:
                raise ValueError("custom reference requires 'fasta' and 'gff'")
        if self.species == "mabs" and not self.exclude_seqids:
            self.exclude_seqids = ["NC_010394.1"]
        if self.species == "mtb" and not self.seqid_map:
            self.seqid_map = {"H37RvBD": "NC_018143.2"}
        return self


class Resources(BaseModel):
    threads: Optional[int] = None


class Config(BaseModel):
    run_name: str
    reference: Reference
    resources: Resources = Field(default_factory=Resources)


def load_config(src) -> Config:
    if isinstance(src, (str, Path)):
        src = yaml.safe_load(Path(src).read_text())
    return Config.model_validate(src)
