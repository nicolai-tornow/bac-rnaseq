import pytest
from pathlib import Path
from engine.enrichment.genesets import load_genesets, align_hits
from tests.testdata import MABS_CATEGORIES

SHEET = MABS_CATEGORIES
pytestmark = pytest.mark.skipif(not SHEET.exists(), reason="no category sheet")


def test_load_module_sets():
    sets = load_genesets(str(SHEET), "2 Gene Annotations", "Gene ID", "Module")
    assert len(sets) > 20
    assert all(isinstance(v, set) and v for v in sets.values())


def test_align_hits_case_insensitive():
    sets = load_genesets(str(SHEET), "2 Gene Annotations", "Gene ID", "Module")
    all_ids = set().union(*sets.values())
    # lower-cased DE IDs (as in rnaseq_suppressor) must still align
    lowered = [g.lower() for g in list(all_ids)[:200]]
    mapping = align_hits(lowered, all_ids, max_unmatched=0.05)
    assert len(mapping) >= int(0.95 * len(lowered))
