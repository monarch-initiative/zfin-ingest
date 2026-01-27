import uuid

import koza
from biolink_model.datamodel.pydanticmodel_v2 import (
    AgentTypeEnum,
    GenotypeToPhenotypicFeatureAssociation,
    KnowledgeLevelEnum,
)
from loguru import logger

# ZECO ID for standard condition - only process records with this condition
STANDARD_CONDITION = "ZECO:0000103"

# Module-level dictionary for deduplication across transform calls
# Records can have multiple life stages for the same Fish/Publication/ZP combination
seen_records = {}


@koza.transform_record()
def transform_record(koza_transform, row):
    global seen_records

    # Skip normal phenotypes - only process abnormal
    if row["Phenotype Tag"] == "normal":
        return []

    # Build ZP key from 7 phenotype fields
    zp_key_elements = [
        row["Affected Structure or Process 1 subterm ID"],
        row["Post-composed Relationship ID"],
        row["Affected Structure or Process 1 superterm ID"],
        row["Phenotype Keyword ID"],
        row["Affected Structure or Process 2 subterm ID"],
        row["Post-composed Relationship (rel) ID"],
        row["Affected Structure or Process 2 superterm ID"],
    ]

    zp_key = "-".join([element or "0" for element in zp_key_elements])

    # Look up ZP term from eqe2zp mapping
    try:
        zp_term = koza_transform.lookup(zp_key, "iri", "eqe2zp")
    except (KeyError, AttributeError):
        logger.debug("ZP concatenation " + zp_key + " did not match a ZP term")
        return []

    if not zp_term:
        logger.debug("ZP concatenation " + zp_key + " did not match a ZP term")
        return []

    # Look up ZECO environment from pheno_environment_fish mapping
    try:
        zeco_term = koza_transform.lookup(row["Environment ID"], "ZECO Term ID (ZECO:ID)", "pheno_environment_fish")
    except (KeyError, AttributeError):
        logger.debug("Environment ID " + row["Environment ID"] + " did not match a ZECO term")
        return []

    # Skip non-standard conditions
    if zeco_term != STANDARD_CONDITION:
        logger.debug("ZP Environment not standard condition")
        return []

    # Deduplication - avoid duplicate records from multiple life stages
    key = "-".join([row["Fish ID"], row["Publication ID"], zp_term])
    if key in seen_records:
        logger.debug(
            "Duplicate record found presumably for differences in life stages, Record={}, LifeStage={}".format(
                key, row["End Stage Name"]
            )
        )
        return []
    seen_records[key] = True

    # Look up PubMed ID from pub2pubmed mapping (fall back to ZFIN ID)
    zdb_pub_id = row["Publication ID"]
    try:
        pubmed_id = koza_transform.lookup(zdb_pub_id, "pubmed", "pub2pubmed")
        if pubmed_id:
            publication_id = "PMID:" + pubmed_id
        else:
            publication_id = "ZFIN:" + zdb_pub_id
    except (KeyError, AttributeError):
        publication_id = "ZFIN:" + zdb_pub_id

    association = GenotypeToPhenotypicFeatureAssociation(
        id="uuid:" + str(uuid.uuid1()),
        subject="ZFIN:" + row["Fish ID"],
        predicate="biolink:has_phenotype",
        object=zp_term,
        publications=[publication_id],
        aggregator_knowledge_source=["infores:monarchinitiative"],
        primary_knowledge_source="infores:zfin",
        knowledge_level=KnowledgeLevelEnum.knowledge_assertion,
        agent_type=AgentTypeEnum.manual_agent,
    )

    return [association]
