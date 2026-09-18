from __future__ import annotations
from pathlib import Path
from typing import Literal, Optional
import yaml
from pydantic import BaseModel, ConfigDict, Field, model_validator

Strand = Literal["reverse", "forward", "unstranded"]
Species = Literal["mabs", "mtb", "custom"]
Layout = Literal["auto", "mate1_only"]


class _Strict(BaseModel):
    # A typo such as `strandness:` must fail validation, not silently fall back
    # to the default.
    model_config = ConfigDict(extra="forbid")


class Reference(_Strict):
    species: Species
    fasta: Optional[str] = None
    gff: Optional[str] = None
    exclude_seqids: list[str] = Field(default_factory=list)
    strandedness: Strand = "reverse"
    feature_types: list[str] = Field(default_factory=lambda: ["gene"])
    id_attribute: str = "locus_tag"
    seqid_map: dict[str, str] = Field(default_factory=dict)
    # Optional table of structural RNAs (see refs/README.md) for a custom genome.
    # The mabs/mtb bundles ship their own.
    structural_rna: Optional[str] = None

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


class Resources(_Strict):
    threads: Optional[int] = None


class Reads(_Strict):
    # auto: every sample is single-end or every sample is paired-end.
    # mate1_only: run every sample single-end from its fastq_r1 (use this when the
    # sample sheet mixes layouts, so layout is not confounded with condition).
    layout: Layout = "auto"


class Contrast(_Strict):
    name: str
    numerator: str
    denominator: str


class Design(_Strict):
    variable: Literal["condition"] = "condition"
    batch_variable: Optional[Literal["batch"]] = None


class Contrasts(_Strict):
    explicit: list[Contrast] = Field(default_factory=list)
    all_vs_all: bool = False


class Thresholds(_Strict):
    padj: float = 0.05
    log2fc: float = 1.0


class Config(_Strict):
    run_name: str
    reference: Reference
    resources: Resources = Field(default_factory=Resources)
    reads: Reads = Field(default_factory=Reads)
    design: Design = Field(default_factory=Design)
    contrasts: Contrasts = Field(default_factory=Contrasts)
    thresholds: Thresholds = Field(default_factory=Thresholds)


def load_config(src) -> Config:
    if isinstance(src, (str, Path)):
        src = yaml.safe_load(Path(src).read_text())
    return Config.model_validate(src)
