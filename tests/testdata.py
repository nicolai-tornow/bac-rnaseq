"""Where the tests find data.

Bundled files (refs/, tests/fixtures/) always exist. Lab data that cannot ship
with the plugin is found under $BAC_RNASEQ_TESTDATA (a folder that contains e.g.
rnaseq_mabs_media/); tests that need it skip when the variable is unset.
"""
import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
REFS = REPO / "refs"
FIXTURES = Path(__file__).resolve().parent / "fixtures"
MABS_GFF = REFS / "mabs" / "NC_010397.1.gff3"
MABS_CATEGORIES = REFS / "mabs" / "categories.xlsx"
READS_R1 = FIXTURES / "SRR10958838_20k_R1.fastq.gz"
READS_R2 = FIXTURES / "SRR10958838_20k_R2.fastq.gz"

_LAB = os.environ.get("BAC_RNASEQ_TESTDATA")
LAB = Path(_LAB) if _LAB else None


def lab(rel: str) -> Path:
    """Path under $BAC_RNASEQ_TESTDATA, or a path that does not exist when unset."""
    return (LAB / rel) if LAB else Path("/nonexistent") / rel
