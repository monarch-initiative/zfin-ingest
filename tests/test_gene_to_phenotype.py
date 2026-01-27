"""Tests for gene_to_phenotype transform."""

from unittest.mock import MagicMock

import pytest
from biolink_model.datamodel.pydanticmodel_v2 import (
    AgentTypeEnum,
    GeneToPhenotypicFeatureAssociation,
    KnowledgeLevelEnum,
)

from gene_to_phenotype import transform_record


@pytest.fixture
def mock_koza_transform():
    """Create a mock koza transform with lookup method."""
    mock = MagicMock()
    # Mock the lookup method to return a ZP term for a specific key
    mock.lookup.return_value = "ZP:0000001"
    return mock


@pytest.fixture
def abnormal_row():
    """Create a sample row with abnormal phenotype tag."""
    return {
        "ID": "1",
        "Gene Symbol": "brca2",
        "Gene ID": "ZDB-GENE-060510-3",
        "Affected Structure or Process 1 subterm ID": "GO:0006281",
        "Affected Structure or Process 1 subterm Name": "DNA repair",
        "Post-composed Relationship ID": "BFO:0000050",
        "Post-composed Relationship Name": "part of",
        "Affected Structure or Process 1 superterm ID": "ZFA:0001439",
        "Affected Structure or Process 1 superterm Name": "anatomical structure",
        "Phenotype Keyword ID": "PATO:0000460",
        "Phenotype Keyword Name": "abnormal",
        "Phenotype Tag": "abnormal",
        "Affected Structure or Process 2 subterm ID": "",
        "Affected Structure or Process 2 subterm name": "",
        "Post-composed Relationship (rel) ID": "",
        "Post-composed Relationship (rel) Name": "",
        "Affected Structure or Process 2 superterm ID": "",
        "Affected Structure or Process 2 superterm name": "",
        "Fish ID": "ZDB-FISH-150901-1",
        "Fish Display Name": "brca2<sup>hg5/hg5</sup>",
        "Start Stage ID": "ZFS:0000001",
        "End Stage ID": "ZFS:0000002",
        "Fish Environment ID": "ZDB-GENOX-041102-1",
        "Publication ID": "ZDB-PUB-170214-55",
        "Figure ID": "ZDB-FIG-170214-10",
    }


@pytest.fixture
def normal_row():
    """Create a sample row with normal phenotype tag (should be skipped)."""
    return {
        "ID": "2",
        "Gene Symbol": "tp53",
        "Gene ID": "ZDB-GENE-990415-270",
        "Affected Structure or Process 1 subterm ID": "GO:0006915",
        "Affected Structure or Process 1 subterm Name": "apoptosis",
        "Post-composed Relationship ID": "",
        "Post-composed Relationship Name": "",
        "Affected Structure or Process 1 superterm ID": "ZFA:0001439",
        "Affected Structure or Process 1 superterm Name": "anatomical structure",
        "Phenotype Keyword ID": "PATO:0000461",
        "Phenotype Keyword Name": "normal",
        "Phenotype Tag": "normal",
        "Affected Structure or Process 2 subterm ID": "",
        "Affected Structure or Process 2 subterm name": "",
        "Post-composed Relationship (rel) ID": "",
        "Post-composed Relationship (rel) Name": "",
        "Affected Structure or Process 2 superterm ID": "",
        "Affected Structure or Process 2 superterm name": "",
        "Fish ID": "ZDB-FISH-150901-2",
        "Fish Display Name": "tp53<sup>zdf1/zdf1</sup>",
        "Start Stage ID": "ZFS:0000001",
        "End Stage ID": "ZFS:0000002",
        "Fish Environment ID": "ZDB-GENOX-041102-1",
        "Publication ID": "ZDB-PUB-170214-56",
        "Figure ID": "ZDB-FIG-170214-11",
    }


class TestGeneToPhenotype:
    """Tests for the gene_to_phenotype transform."""

    def test_abnormal_phenotype_creates_association(self, mock_koza_transform, abnormal_row):
        """Test that abnormal phenotype tag creates a GeneToPhenotypicFeatureAssociation."""
        result = transform_record(mock_koza_transform, abnormal_row)

        assert len(result) == 1
        association = result[0]

        assert isinstance(association, GeneToPhenotypicFeatureAssociation)
        assert association.subject == "ZFIN:ZDB-GENE-060510-3"
        assert association.predicate == "biolink:has_phenotype"
        assert association.object == "ZP:0000001"
        assert association.publications == ["ZFIN:ZDB-PUB-170214-55"]
        assert association.aggregator_knowledge_source == ["infores:monarchinitiative"]
        assert association.primary_knowledge_source == "infores:zfin"
        assert association.knowledge_level == KnowledgeLevelEnum.knowledge_assertion
        assert association.agent_type == AgentTypeEnum.manual_agent
        assert association.id.startswith("uuid:")

    def test_normal_phenotype_skipped(self, mock_koza_transform, normal_row):
        """Test that normal phenotype tag is skipped."""
        result = transform_record(mock_koza_transform, normal_row)

        assert result == []

    def test_lookup_key_construction(self, mock_koza_transform, abnormal_row):
        """Test that the ZP lookup key is constructed correctly."""
        transform_record(mock_koza_transform, abnormal_row)

        # The key should be constructed from the 7 elements, with empty values replaced by "0"
        expected_key = "GO:0006281-BFO:0000050-ZFA:0001439-PATO:0000460-0-0-0"
        mock_koza_transform.lookup.assert_called_once_with(expected_key, "iri", "eqe2zp")

    def test_lookup_failure_returns_empty(self, mock_koza_transform, abnormal_row):
        """Test that lookup failure returns empty list."""
        mock_koza_transform.lookup.side_effect = KeyError("Key not found")

        result = transform_record(mock_koza_transform, abnormal_row)

        assert result == []

    def test_lookup_returns_none(self, mock_koza_transform, abnormal_row):
        """Test that None lookup result returns empty list."""
        mock_koza_transform.lookup.return_value = None

        result = transform_record(mock_koza_transform, abnormal_row)

        assert result == []

    def test_attribute_error_returns_empty(self, mock_koza_transform, abnormal_row):
        """Test that AttributeError during lookup returns empty list."""
        mock_koza_transform.lookup.side_effect = AttributeError("No lookup method")

        result = transform_record(mock_koza_transform, abnormal_row)

        assert result == []

    def test_all_fields_populated_in_key(self, mock_koza_transform):
        """Test key construction when all fields are populated."""
        row = {
            "ID": "3",
            "Gene Symbol": "fgf8a",
            "Gene ID": "ZDB-GENE-990415-72",
            "Affected Structure or Process 1 subterm ID": "GO:0001",
            "Affected Structure or Process 1 subterm Name": "subterm1",
            "Post-composed Relationship ID": "REL:0001",
            "Post-composed Relationship Name": "relationship1",
            "Affected Structure or Process 1 superterm ID": "ZFA:0001",
            "Affected Structure or Process 1 superterm Name": "superterm1",
            "Phenotype Keyword ID": "PATO:0001",
            "Phenotype Keyword Name": "keyword1",
            "Phenotype Tag": "abnormal",
            "Affected Structure or Process 2 subterm ID": "GO:0002",
            "Affected Structure or Process 2 subterm name": "subterm2",
            "Post-composed Relationship (rel) ID": "REL:0002",
            "Post-composed Relationship (rel) Name": "relationship2",
            "Affected Structure or Process 2 superterm ID": "ZFA:0002",
            "Affected Structure or Process 2 superterm name": "superterm2",
            "Fish ID": "ZDB-FISH-150901-3",
            "Fish Display Name": "fgf8a<sup>ti282/ti282</sup>",
            "Start Stage ID": "ZFS:0000001",
            "End Stage ID": "ZFS:0000002",
            "Fish Environment ID": "ZDB-GENOX-041102-1",
            "Publication ID": "ZDB-PUB-170214-57",
            "Figure ID": "ZDB-FIG-170214-12",
        }

        transform_record(mock_koza_transform, row)

        expected_key = "GO:0001-REL:0001-ZFA:0001-PATO:0001-GO:0002-REL:0002-ZFA:0002"
        mock_koza_transform.lookup.assert_called_once_with(expected_key, "iri", "eqe2zp")
