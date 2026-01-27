"""Tests for genotype_to_phenotype transform."""

from unittest.mock import MagicMock

import pytest
from biolink_model.datamodel.pydanticmodel_v2 import (
    AgentTypeEnum,
    GenotypeToPhenotypicFeatureAssociation,
    KnowledgeLevelEnum,
)

import genotype_to_phenotype
from genotype_to_phenotype import transform_record


@pytest.fixture(autouse=True)
def reset_seen_records():
    """Reset the seen_records dict before each test."""
    genotype_to_phenotype.seen_records = {}
    yield
    genotype_to_phenotype.seen_records = {}


@pytest.fixture
def mock_koza_transform():
    """Create a mock koza transform with lookup method."""
    mock = MagicMock()

    def mock_lookup(key, field, mapping_name):
        if mapping_name == "eqe2zp":
            return "ZP:0000001"
        elif mapping_name == "pheno_environment_fish":
            return "ZECO:0000103"  # Standard condition
        elif mapping_name == "pub2pubmed":
            return "12345678"
        raise KeyError(f"Unknown mapping: {mapping_name}")

    mock.lookup.side_effect = mock_lookup
    return mock


@pytest.fixture
def abnormal_row():
    """Create a sample row with abnormal phenotype tag and standard condition."""
    return {
        "Fish ID": "ZDB-FISH-150901-1",
        "Fish Name": "brca2<sup>hg5/hg5</sup>",
        "Start Stage ID": "ZFS:0000001",
        "Start Stage Name": "Zygote:1-cell",
        "End Stage ID": "ZFS:0000002",
        "End Stage Name": "Cleavage:2-cell",
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
        "Publication ID": "ZDB-PUB-170214-55",
        "Environment ID": "ZDB-EXP-041102-1",
    }


@pytest.fixture
def normal_row():
    """Create a sample row with normal phenotype tag (should be skipped)."""
    return {
        "Fish ID": "ZDB-FISH-150901-2",
        "Fish Name": "tp53<sup>zdf1/zdf1</sup>",
        "Start Stage ID": "ZFS:0000001",
        "Start Stage Name": "Zygote:1-cell",
        "End Stage ID": "ZFS:0000002",
        "End Stage Name": "Cleavage:2-cell",
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
        "Publication ID": "ZDB-PUB-170214-56",
        "Environment ID": "ZDB-EXP-041102-1",
    }


class TestGenotypeToPhenotype:
    """Tests for the genotype_to_phenotype transform."""

    def test_abnormal_phenotype_creates_association(self, mock_koza_transform, abnormal_row):
        """Test that abnormal phenotype tag creates a GenotypeToPhenotypicFeatureAssociation."""
        result = transform_record(mock_koza_transform, abnormal_row)

        assert len(result) == 1
        association = result[0]

        assert isinstance(association, GenotypeToPhenotypicFeatureAssociation)
        assert association.subject == "ZFIN:ZDB-FISH-150901-1"
        assert association.predicate == "biolink:has_phenotype"
        assert association.object == "ZP:0000001"
        assert association.publications == ["PMID:12345678"]
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
        calls = mock_koza_transform.lookup.call_args_list
        eqe2zp_call = [c for c in calls if c[0][2] == "eqe2zp"][0]
        assert eqe2zp_call[0][0] == expected_key

    def test_zp_lookup_failure_returns_empty(self, mock_koza_transform, abnormal_row):
        """Test that ZP lookup failure returns empty list."""

        def mock_lookup(key, field, mapping_name):
            if mapping_name == "eqe2zp":
                raise KeyError("Key not found")
            return "ZECO:0000103"

        mock_koza_transform.lookup.side_effect = mock_lookup

        result = transform_record(mock_koza_transform, abnormal_row)

        assert result == []

    def test_zp_lookup_returns_none(self, mock_koza_transform, abnormal_row):
        """Test that None ZP lookup result returns empty list."""

        def mock_lookup(key, field, mapping_name):
            if mapping_name == "eqe2zp":
                return None
            return "ZECO:0000103"

        mock_koza_transform.lookup.side_effect = mock_lookup

        result = transform_record(mock_koza_transform, abnormal_row)

        assert result == []

    def test_non_standard_condition_skipped(self, mock_koza_transform, abnormal_row):
        """Test that non-standard environment conditions are skipped."""

        def mock_lookup(key, field, mapping_name):
            if mapping_name == "eqe2zp":
                return "ZP:0000001"
            elif mapping_name == "pheno_environment_fish":
                return "ZECO:0000104"  # Non-standard condition
            return "12345678"

        mock_koza_transform.lookup.side_effect = mock_lookup

        result = transform_record(mock_koza_transform, abnormal_row)

        assert result == []

    def test_environment_lookup_failure_returns_empty(self, mock_koza_transform, abnormal_row):
        """Test that environment lookup failure returns empty list."""

        def mock_lookup(key, field, mapping_name):
            if mapping_name == "eqe2zp":
                return "ZP:0000001"
            elif mapping_name == "pheno_environment_fish":
                raise KeyError("Key not found")
            return "12345678"

        mock_koza_transform.lookup.side_effect = mock_lookup

        result = transform_record(mock_koza_transform, abnormal_row)

        assert result == []

    def test_deduplication(self, mock_koza_transform, abnormal_row):
        """Test that duplicate records are skipped."""
        # First call should return an association
        result1 = transform_record(mock_koza_transform, abnormal_row)
        assert len(result1) == 1

        # Second call with same Fish ID, Publication ID, and ZP term should be skipped
        result2 = transform_record(mock_koza_transform, abnormal_row)
        assert result2 == []

    def test_deduplication_different_life_stages(self, mock_koza_transform, abnormal_row):
        """Test that records with different life stages but same key are deduplicated."""
        # First call
        result1 = transform_record(mock_koza_transform, abnormal_row)
        assert len(result1) == 1

        # Second call with different life stage
        abnormal_row["End Stage Name"] = "Larval:Day 5"
        abnormal_row["End Stage ID"] = "ZFS:0000010"
        result2 = transform_record(mock_koza_transform, abnormal_row)
        assert result2 == []

    def test_pubmed_fallback_to_zfin(self, mock_koza_transform, abnormal_row):
        """Test that publication falls back to ZFIN ID when PubMed lookup fails."""

        def mock_lookup(key, field, mapping_name):
            if mapping_name == "eqe2zp":
                return "ZP:0000001"
            elif mapping_name == "pheno_environment_fish":
                return "ZECO:0000103"
            elif mapping_name == "pub2pubmed":
                raise KeyError("Key not found")
            return None

        mock_koza_transform.lookup.side_effect = mock_lookup

        result = transform_record(mock_koza_transform, abnormal_row)

        assert len(result) == 1
        assert result[0].publications == ["ZFIN:ZDB-PUB-170214-55"]

    def test_pubmed_fallback_when_none(self, mock_koza_transform, abnormal_row):
        """Test that publication falls back to ZFIN ID when PubMed lookup returns None."""

        def mock_lookup(key, field, mapping_name):
            if mapping_name == "eqe2zp":
                return "ZP:0000001"
            elif mapping_name == "pheno_environment_fish":
                return "ZECO:0000103"
            elif mapping_name == "pub2pubmed":
                return None
            return None

        mock_koza_transform.lookup.side_effect = mock_lookup

        result = transform_record(mock_koza_transform, abnormal_row)

        assert len(result) == 1
        assert result[0].publications == ["ZFIN:ZDB-PUB-170214-55"]

    def test_all_fields_populated_in_key(self, mock_koza_transform):
        """Test key construction when all fields are populated."""
        row = {
            "Fish ID": "ZDB-FISH-150901-3",
            "Fish Name": "fgf8a<sup>ti282/ti282</sup>",
            "Start Stage ID": "ZFS:0000001",
            "Start Stage Name": "Zygote:1-cell",
            "End Stage ID": "ZFS:0000002",
            "End Stage Name": "Cleavage:2-cell",
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
            "Publication ID": "ZDB-PUB-170214-57",
            "Environment ID": "ZDB-EXP-041102-1",
        }

        transform_record(mock_koza_transform, row)

        expected_key = "GO:0001-REL:0001-ZFA:0001-PATO:0001-GO:0002-REL:0002-ZFA:0002"
        calls = mock_koza_transform.lookup.call_args_list
        eqe2zp_call = [c for c in calls if c[0][2] == "eqe2zp"][0]
        assert eqe2zp_call[0][0] == expected_key

    def test_attribute_error_on_zp_lookup_returns_empty(self, mock_koza_transform, abnormal_row):
        """Test that AttributeError during ZP lookup returns empty list."""

        def mock_lookup(key, field, mapping_name):
            if mapping_name == "eqe2zp":
                raise AttributeError("No lookup method")
            return "ZECO:0000103"

        mock_koza_transform.lookup.side_effect = mock_lookup

        result = transform_record(mock_koza_transform, abnormal_row)

        assert result == []

    def test_attribute_error_on_environment_lookup_returns_empty(self, mock_koza_transform, abnormal_row):
        """Test that AttributeError during environment lookup returns empty list."""

        def mock_lookup(key, field, mapping_name):
            if mapping_name == "eqe2zp":
                return "ZP:0000001"
            elif mapping_name == "pheno_environment_fish":
                raise AttributeError("No lookup method")
            return "12345678"

        mock_koza_transform.lookup.side_effect = mock_lookup

        result = transform_record(mock_koza_transform, abnormal_row)

        assert result == []
