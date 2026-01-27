"""Tests for orthology transform."""

from unittest.mock import MagicMock

import pytest
from biolink_model.datamodel.pydanticmodel_v2 import (
    AgentTypeEnum,
    GeneToGeneHomologyAssociation,
    KnowledgeLevelEnum,
)

from orthology import EVIDENCE_MAP, transform_record


@pytest.fixture
def mock_koza_transform():
    """Create a mock koza transform."""
    return MagicMock()


@pytest.fixture
def fly_ortholog_row():
    """Create a sample row for a fly ortholog."""
    return {
        "zfin_gene": "ZFIN:ZDB-GENE-000112-47",
        "ortholog_gene": "FB:FBgn0000017",
        "evidence": "AA",
        "publications": "ZDB-PUB-020723-7|ZDB-PUB-030508-1",
    }


@pytest.fixture
def mouse_ortholog_row():
    """Create a sample row for a mouse ortholog."""
    return {
        "zfin_gene": "ZFIN:ZDB-GENE-000112-47",
        "ortholog_gene": "MGI:97490",
        "evidence": "PT",
        "publications": "ZDB-PUB-020723-7",
    }


@pytest.fixture
def human_ortholog_row():
    """Create a sample row for a human ortholog."""
    return {
        "zfin_gene": "ZFIN:ZDB-GENE-000112-47",
        "ortholog_gene": "HGNC:7",
        "evidence": "NT",
        "publications": "ZDB-PUB-020723-7",
    }


@pytest.fixture
def row_no_publications():
    """Create a sample row with no publications."""
    return {
        "zfin_gene": "ZFIN:ZDB-GENE-000112-47",
        "ortholog_gene": "FB:FBgn0000017",
        "evidence": "AA",
        "publications": "",
    }


@pytest.fixture
def row_unknown_evidence():
    """Create a sample row with unknown evidence code."""
    return {
        "zfin_gene": "ZFIN:ZDB-GENE-000112-47",
        "ortholog_gene": "FB:FBgn0000017",
        "evidence": "XX",
        "publications": "ZDB-PUB-020723-7",
    }


class TestOrthology:
    """Tests for the orthology transform."""

    def test_fly_ortholog_creates_association(self, mock_koza_transform, fly_ortholog_row):
        """Test that a fly ortholog row creates a GeneToGeneHomologyAssociation."""
        result = transform_record(mock_koza_transform, fly_ortholog_row)

        assert len(result) == 1
        association = result[0]

        assert isinstance(association, GeneToGeneHomologyAssociation)
        assert association.subject == "ZFIN:ZDB-GENE-000112-47"
        assert association.predicate == "biolink:orthologous_to"
        assert association.object == "FB:FBgn0000017"
        assert association.has_evidence == ["ECO:0000031"]  # AA = Amino acid sequence comparison
        assert association.publications == ["ZFIN:ZDB-PUB-020723-7", "ZFIN:ZDB-PUB-030508-1"]
        assert association.aggregator_knowledge_source == ["infores:monarchinitiative"]
        assert association.primary_knowledge_source == "infores:zfin"
        assert association.knowledge_level == KnowledgeLevelEnum.knowledge_assertion
        assert association.agent_type == AgentTypeEnum.manual_agent
        assert association.id.startswith("uuid:")

    def test_mouse_ortholog_creates_association(self, mock_koza_transform, mouse_ortholog_row):
        """Test that a mouse ortholog row creates a GeneToGeneHomologyAssociation."""
        result = transform_record(mock_koza_transform, mouse_ortholog_row)

        assert len(result) == 1
        association = result[0]

        assert isinstance(association, GeneToGeneHomologyAssociation)
        assert association.subject == "ZFIN:ZDB-GENE-000112-47"
        assert association.object == "MGI:97490"
        assert association.has_evidence == ["ECO:0007750"]  # PT = Phylogenetic tree

    def test_human_ortholog_creates_association(self, mock_koza_transform, human_ortholog_row):
        """Test that a human ortholog row creates a GeneToGeneHomologyAssociation."""
        result = transform_record(mock_koza_transform, human_ortholog_row)

        assert len(result) == 1
        association = result[0]

        assert isinstance(association, GeneToGeneHomologyAssociation)
        assert association.subject == "ZFIN:ZDB-GENE-000112-47"
        assert association.object == "HGNC:7"
        assert association.has_evidence == ["ECO:0000032"]  # NT = Nucleotide sequence comparison

    def test_multiple_publications_split_correctly(self, mock_koza_transform, fly_ortholog_row):
        """Test that multiple publications are correctly parsed."""
        result = transform_record(mock_koza_transform, fly_ortholog_row)

        assert len(result) == 1
        assert result[0].publications == ["ZFIN:ZDB-PUB-020723-7", "ZFIN:ZDB-PUB-030508-1"]

    def test_single_publication(self, mock_koza_transform, mouse_ortholog_row):
        """Test that single publication is correctly parsed."""
        result = transform_record(mock_koza_transform, mouse_ortholog_row)

        assert len(result) == 1
        assert result[0].publications == ["ZFIN:ZDB-PUB-020723-7"]

    def test_empty_publications(self, mock_koza_transform, row_no_publications):
        """Test that empty publications field results in None."""
        result = transform_record(mock_koza_transform, row_no_publications)

        assert len(result) == 1
        assert result[0].publications is None

    def test_unknown_evidence_code(self, mock_koza_transform, row_unknown_evidence):
        """Test that unknown evidence code results in None has_evidence."""
        result = transform_record(mock_koza_transform, row_unknown_evidence)

        assert len(result) == 1
        assert result[0].has_evidence is None

    def test_all_evidence_codes_mapped(self):
        """Test that all expected evidence codes are in the mapping."""
        expected_codes = ["AA", "CE", "CL", "FC", "NT", "PT", "OT"]
        for code in expected_codes:
            assert code in EVIDENCE_MAP

    def test_evidence_code_aa(self, mock_koza_transform, fly_ortholog_row):
        """Test AA evidence code mapping."""
        fly_ortholog_row["evidence"] = "AA"
        result = transform_record(mock_koza_transform, fly_ortholog_row)
        assert result[0].has_evidence == ["ECO:0000031"]

    def test_evidence_code_ce(self, mock_koza_transform, fly_ortholog_row):
        """Test CE evidence code mapping."""
        fly_ortholog_row["evidence"] = "CE"
        result = transform_record(mock_koza_transform, fly_ortholog_row)
        assert result[0].has_evidence == ["ECO:0001163"]

    def test_evidence_code_cl(self, mock_koza_transform, fly_ortholog_row):
        """Test CL evidence code mapping."""
        fly_ortholog_row["evidence"] = "CL"
        result = transform_record(mock_koza_transform, fly_ortholog_row)
        assert result[0].has_evidence == ["ECO:0000354"]

    def test_evidence_code_fc(self, mock_koza_transform, fly_ortholog_row):
        """Test FC evidence code mapping."""
        fly_ortholog_row["evidence"] = "FC"
        result = transform_record(mock_koza_transform, fly_ortholog_row)
        assert result[0].has_evidence == ["ECO:0006091"]

    def test_evidence_code_nt(self, mock_koza_transform, fly_ortholog_row):
        """Test NT evidence code mapping."""
        fly_ortholog_row["evidence"] = "NT"
        result = transform_record(mock_koza_transform, fly_ortholog_row)
        assert result[0].has_evidence == ["ECO:0000032"]

    def test_evidence_code_pt(self, mock_koza_transform, fly_ortholog_row):
        """Test PT evidence code mapping."""
        fly_ortholog_row["evidence"] = "PT"
        result = transform_record(mock_koza_transform, fly_ortholog_row)
        assert result[0].has_evidence == ["ECO:0007750"]

    def test_evidence_code_ot(self, mock_koza_transform, fly_ortholog_row):
        """Test OT evidence code mapping."""
        fly_ortholog_row["evidence"] = "OT"
        result = transform_record(mock_koza_transform, fly_ortholog_row)
        assert result[0].has_evidence == ["ECO:0000352"]
