import uuid

import koza
from biolink_model.datamodel.pydanticmodel_v2 import (
    AgentTypeEnum,
    GeneToGeneHomologyAssociation,
    KnowledgeLevelEnum,
)

# Evidence code mapping from ZFIN evidence codes to ECO terms
EVIDENCE_MAP = {
    "AA": "ECO:0000031",  # Amino acid sequence comparison
    "CE": "ECO:0001163",  # Coincident expression
    "CL": "ECO:0000354",  # Conserved map location
    "FC": "ECO:0006091",  # Functional complementation
    "NT": "ECO:0000032",  # Nucleotide sequence comparison
    "PT": "ECO:0007750",  # Phylogenetic tree
    "OT": "ECO:0000352",  # Other
}


@koza.transform_record()
def transform_record(koza_transform, row):
    # Parse publications - split on pipe and prefix with ZFIN:
    publications = None
    if row["publications"]:
        publications = [f"ZFIN:{pub}" for pub in row["publications"].split("|")]

    # Map evidence code to ECO term
    evidence = EVIDENCE_MAP.get(row["evidence"], None)
    has_evidence = [evidence] if evidence else None

    association = GeneToGeneHomologyAssociation(
        id="uuid:" + str(uuid.uuid1()),
        subject=row["zfin_gene"],
        predicate="biolink:orthologous_to",
        object=row["ortholog_gene"],
        has_evidence=has_evidence,
        publications=publications,
        primary_knowledge_source="infores:zfin",
        aggregator_knowledge_source=["infores:monarchinitiative"],
        knowledge_level=KnowledgeLevelEnum.knowledge_assertion,
        agent_type=AgentTypeEnum.manual_agent,
    )

    return [association]
