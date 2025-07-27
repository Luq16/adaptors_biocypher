from __future__ import annotations

import random
import string
from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class InterproAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    DOMAIN = auto()
    FAMILY = auto()
    SITE = auto()


class InterproAdapterDomainField(Enum):
    """
    Define possible fields the adapter can provide for protein domains.
    """
    NAME = "name"
    DESCRIPTION = "description"
    TYPE = "type"
    INTERPRO_ID = "interpro_id"
    PFAM_ID = "pfam_id"
    PROSITE_ID = "prosite_id"
    SMART_ID = "smart_id"
    SUPFAM_ID = "supfam_id"
    LENGTH = "length"
    ABSTRACT = "abstract"


class InterproAdapterEdgeType(Enum):
    """
    Enum for the types of edges the adapter can provide.
    """
    PROTEIN_HAS_DOMAIN = "protein_has_domain"
    DOMAIN_IS_PART_OF_FAMILY = "domain_is_part_of_family"
    DOMAIN_CONTAINS_SITE = "domain_contains_site"
    DOMAIN_RELATED_TO_DOMAIN = "domain_related_to_domain"


class InterproAdapterProteinDomainEdgeField(Enum):
    """
    Define possible fields for protein-domain annotation edges.
    """
    START_POSITION = "start_position"
    END_POSITION = "end_position"
    E_VALUE = "e_value"
    SCORE = "score"
    SOURCE = "source"
    MATCH_STATUS = "match_status"
    SIGNATURE_LIBRARY = "signature_library"
    MODEL_ACCESSION = "model_accession"
    INTERPRO_ACCESSION = "interpro_accession"
    EVIDENCE = "evidence"


class InterproAdapter:
    """
    InterPro BioCypher adapter. Generates protein domain/family nodes and 
    protein-domain annotation edges for protein functional classification.

    Args:
        node_types: List of node types to include in the result.
        node_fields: List of node fields to include in the result.
        edge_types: List of edge types to include in the result.
        edge_fields: List of edge fields to include in the result.
        test_mode: If True, limits the number of entries for testing.
    """

    def __init__(
        self,
        node_types: Optional[list] = None,
        node_fields: Optional[list] = None,
        edge_types: Optional[list] = None,
        edge_fields: Optional[list] = None,
        test_mode: Optional[bool] = False,
    ):
        self.test_mode = test_mode
        self._set_types_and_fields(node_types, node_fields, edge_types, edge_fields)

    def get_nodes(self):
        """
        Returns a generator of node tuples for node types specified in the
        adapter constructor.
        """
        logger.info("Generating InterPro domain/family nodes.")

        self.nodes = []
        self.domains = []
        self.families = []
        self.sites = []

        # Generate domains
        if InterproAdapterNodeType.DOMAIN in self.node_types:
            num_domains = 25 if self.test_mode else 250
            for _ in range(num_domains):
                domain_node = InterproDomainNode(fields=self.node_fields, node_type="domain")
                self.domains.append(domain_node)
                self.nodes.append(domain_node)

        # Generate families
        if InterproAdapterNodeType.FAMILY in self.node_types:
            num_families = 15 if self.test_mode else 150
            for _ in range(num_families):
                family_node = InterproDomainNode(fields=self.node_fields, node_type="family")
                self.families.append(family_node)
                self.nodes.append(family_node)

        # Generate sites
        if InterproAdapterNodeType.SITE in self.node_types:
            num_sites = 10 if self.test_mode else 100
            for _ in range(num_sites):
                site_node = InterproDomainNode(fields=self.node_fields, node_type="site")
                self.sites.append(site_node)
                self.nodes.append(site_node)

        for node in self.nodes:
            yield (node.get_id(), node.get_label(), node.get_properties())

    def get_edges(self, target_proteins=None):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.

        Args:
            target_proteins: List of protein nodes to create domain annotations for.
        """
        logger.info("Generating InterPro annotation edges.")

        if not self.nodes:
            list(self.get_nodes())  # Generate nodes first

        # Protein-domain annotation relationships
        if (InterproAdapterEdgeType.PROTEIN_HAS_DOMAIN in self.edge_types 
            and target_proteins):
            
            for protein in target_proteins:
                # Each protein has a 40% chance of having domain annotations
                if random.random() < 0.4:
                    # Select 1-4 random domains for this protein
                    num_domains = random.randint(1, 4)
                    available_domains = [node for node in self.nodes if node.get_label() in ["domain", "family", "site"]]
                    
                    if available_domains:
                        selected_domains = random.sample(available_domains, min(num_domains, len(available_domains)))
                        
                        for domain in selected_domains:
                            protein_id = protein[0] if isinstance(protein, tuple) else str(protein)
                            
                            relationship_id = "".join(
                                random.choice(string.ascii_letters + string.digits)
                                for _ in range(12)
                            )
                            
                            properties = {}
                            for field in self.edge_fields:
                                if field == InterproAdapterProteinDomainEdgeField.START_POSITION:
                                    properties["start_position"] = random.randint(1, 500)
                                elif field == InterproAdapterProteinDomainEdgeField.END_POSITION:
                                    start = properties.get("start_position", random.randint(1, 500))
                                    properties["end_position"] = start + random.randint(20, 200)
                                elif field == InterproAdapterProteinDomainEdgeField.E_VALUE:
                                    properties["e_value"] = f"{random.uniform(1e-10, 1e-3):.2e}"
                                elif field == InterproAdapterProteinDomainEdgeField.SCORE:
                                    properties["score"] = round(random.uniform(20.0, 500.0), 1)
                                elif field == InterproAdapterProteinDomainEdgeField.SOURCE:
                                    properties["source"] = random.choice([
                                        "InterPro", "Pfam", "PROSITE", "SMART", "SUPERFAMILY",
                                        "PRINTS", "PANTHER", "TIGRFAMs", "PIRSF", "HAMAP"
                                    ])
                                elif field == InterproAdapterProteinDomainEdgeField.MATCH_STATUS:
                                    properties["match_status"] = random.choice([
                                        "T", "F", "?"  # True, False, Unknown
                                    ])
                                elif field == InterproAdapterProteinDomainEdgeField.SIGNATURE_LIBRARY:
                                    properties["signature_library"] = random.choice([
                                        "Pfam", "PROSITE", "SMART", "SUPERFAMILY", "PRINTS",
                                        "PANTHER", "TIGRFAMs", "PIRSF", "HAMAP", "CDD"
                                    ])
                                elif field == InterproAdapterProteinDomainEdgeField.MODEL_ACCESSION:
                                    libraries = {
                                        "Pfam": f"PF{random.randint(10000, 99999):05d}",
                                        "PROSITE": f"PS{random.randint(10000, 99999):05d}",
                                        "SMART": f"SM{random.randint(10000, 99999):05d}",
                                        "SUPERFAMILY": f"SSF{random.randint(10000, 99999)}"
                                    }
                                    library = random.choice(list(libraries.keys()))
                                    properties["model_accession"] = libraries[library]
                                elif field == InterproAdapterProteinDomainEdgeField.INTERPRO_ACCESSION:
                                    properties["interpro_accession"] = f"IPR{random.randint(100000, 999999):06d}"
                                elif field == InterproAdapterProteinDomainEdgeField.EVIDENCE:
                                    properties["evidence"] = random.choice([
                                        "computational", "experimental", "manual assertion",
                                        "automatic assertion", "combinatorial evidence"
                                    ])

                            yield (
                                relationship_id,
                                protein_id,
                                domain.get_id(),
                                InterproAdapterEdgeType.PROTEIN_HAS_DOMAIN.value,
                                properties,
                            )

        # Domain-Family relationships
        if InterproAdapterEdgeType.DOMAIN_IS_PART_OF_FAMILY in self.edge_types:
            for domain in self.domains:
                # 25% chance of belonging to a family
                if random.random() < 0.25 and self.families:
                    family = random.choice(self.families)
                    
                    relationship_id = "".join(
                        random.choice(string.ascii_letters + string.digits)
                        for _ in range(12)
                    )
                    
                    yield (
                        relationship_id,
                        domain.get_id(),
                        family.get_id(),
                        InterproAdapterEdgeType.DOMAIN_IS_PART_OF_FAMILY.value,
                        {"relationship_type": "is_part_of"},
                    )

        # Domain-Site relationships
        if InterproAdapterEdgeType.DOMAIN_CONTAINS_SITE in self.edge_types:
            for domain in self.domains:
                # 20% chance of containing sites
                if random.random() < 0.2 and self.sites:
                    num_sites = random.randint(1, 2)
                    selected_sites = random.sample(self.sites, min(num_sites, len(self.sites)))
                    
                    for site in selected_sites:
                        relationship_id = "".join(
                            random.choice(string.ascii_letters + string.digits)
                            for _ in range(12)
                        )
                        
                        yield (
                            relationship_id,
                            domain.get_id(),
                            site.get_id(),
                            InterproAdapterEdgeType.DOMAIN_CONTAINS_SITE.value,
                            {"relationship_type": "contains"},
                        )

        # Domain-Domain relationships
        if InterproAdapterEdgeType.DOMAIN_RELATED_TO_DOMAIN in self.edge_types:
            for i, domain_a in enumerate(self.domains):
                # 10% chance of being related to another domain
                if random.random() < 0.1 and i < len(self.domains) - 1:
                    domain_b = random.choice(self.domains[i+1:])
                    
                    relationship_id = "".join(
                        random.choice(string.ascii_letters + string.digits)
                        for _ in range(12)
                    )
                    
                    relationship_type = random.choice([
                        "homologous_to", "similar_to", "related_to", "paralogous_to"
                    ])
                    
                    yield (
                        relationship_id,
                        domain_a.get_id(),
                        domain_b.get_id(),
                        InterproAdapterEdgeType.DOMAIN_RELATED_TO_DOMAIN.value,
                        {"relationship_type": relationship_type},
                    )

    def get_node_count(self):
        """
        Returns the number of nodes generated by the adapter.
        """
        return len(list(self.get_nodes()))

    def _set_types_and_fields(self, node_types, node_fields, edge_types, edge_fields):
        if node_types:
            self.node_types = node_types
        else:
            self.node_types = [type for type in InterproAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [field for field in InterproAdapterDomainField]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in InterproAdapterEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in InterproAdapterProteinDomainEdgeField]


class Node:
    """
    Base class for nodes.
    """

    def __init__(self):
        self.id = None
        self.label = None
        self.properties = {}

    def get_id(self):
        """
        Returns the node id.
        """
        return self.id

    def get_label(self):
        """
        Returns the node label.
        """
        return self.label

    def get_properties(self):
        """
        Returns the node properties.
        """
        return self.properties


class InterproDomainNode(Node):
    """
    Generates instances of InterPro domain/family/site nodes.
    """

    def __init__(self, fields: Optional[list] = None, node_type: str = "domain"):
        self.fields = fields
        self.node_type = node_type
        self.id = self._generate_id()
        self.label = node_type
        self.properties = self._generate_properties()

    def _generate_id(self):
        """
        Generate an InterPro-style identifier.
        """
        # InterPro format: IPR000001 to IPR999999
        return f"IPR{random.randint(100000, 999999):06d}"

    def _generate_properties(self):
        properties = {}

        # Generate name based on node type
        if InterproAdapterDomainField.NAME in self.fields:
            if self.node_type == "domain":
                domain_names = [
                    "Immunoglobulin-like fold", "EGF-like domain", "SH3 domain",
                    "Kinase domain", "DNA-binding domain", "Zinc finger domain",
                    "Helix-turn-helix motif", "Leucine zipper", "WD40 repeat",
                    "Ankyrin repeat", "Armadillo repeat", "RING finger domain",
                    "Bromodomain", "Chromodomain", "Tudor domain", "PDZ domain",
                    "Pleckstrin homology domain", "Src homology 2 domain",
                    "Death domain", "Caspase domain", "Thioredoxin fold",
                    "Rossmann fold", "TIM barrel", "Beta barrel"
                ]
                properties["name"] = random.choice(domain_names)
            elif self.node_type == "family":
                family_names = [
                    "Protein kinase family", "Transcription factor family",
                    "Immunoglobulin family", "G-protein coupled receptor family",
                    "Ion channel family", "Protease family", "Oxidase family",
                    "Transferase family", "Hydrolase family", "Ligase family",
                    "Ribosomal protein family", "Histone family", "Globin family",
                    "Cytochrome family", "Hemoglobin family"
                ]
                properties["name"] = random.choice(family_names)
            else:  # site
                site_names = [
                    "Active site", "Binding site", "Catalytic site", "Cleavage site",
                    "Metal binding site", "Phosphorylation site", "Glycosylation site",
                    "Ubiquitination site", "Sumoylation site", "Acetylation site",
                    "Methylation site", "Nuclear localization signal",
                    "Signal peptide", "Transmembrane region"
                ]
                properties["name"] = random.choice(site_names)

        # Generate description
        if InterproAdapterDomainField.DESCRIPTION in self.fields:
            if self.node_type == "domain":
                properties["description"] = "Conserved protein domain with specific structural and functional characteristics"
            elif self.node_type == "family":
                properties["description"] = "Group of related proteins sharing common evolutionary origin and function"
            else:  # site
                properties["description"] = "Functionally important region within protein sequence"

        # Generate type
        if InterproAdapterDomainField.TYPE in self.fields:
            properties["type"] = self.node_type

        # Generate InterPro ID
        if InterproAdapterDomainField.INTERPRO_ID in self.fields:
            properties["interpro_id"] = self.id

        # Generate cross-references
        if InterproAdapterDomainField.PFAM_ID in self.fields:
            if random.random() < 0.7:  # 70% chance
                properties["pfam_id"] = f"PF{random.randint(10000, 99999):05d}"

        if InterproAdapterDomainField.PROSITE_ID in self.fields:
            if random.random() < 0.4:  # 40% chance
                properties["prosite_id"] = f"PS{random.randint(10000, 99999):05d}"

        if InterproAdapterDomainField.SMART_ID in self.fields:
            if random.random() < 0.3:  # 30% chance
                properties["smart_id"] = f"SM{random.randint(10000, 99999):05d}"

        if InterproAdapterDomainField.SUPFAM_ID in self.fields:
            if random.random() < 0.2:  # 20% chance
                properties["supfam_id"] = f"SSF{random.randint(10000, 99999)}"

        # Generate length (for domains)
        if InterproAdapterDomainField.LENGTH in self.fields and self.node_type == "domain":
            properties["length"] = random.randint(20, 300)

        # Generate abstract
        if InterproAdapterDomainField.ABSTRACT in self.fields:
            abstracts = [
                "This entry represents a conserved domain found in various proteins",
                "Members of this family are involved in protein-protein interactions",
                "This domain is characterized by a specific fold and binding properties",
                "Proteins containing this domain are typically involved in signal transduction",
                "This family includes enzymes with similar catalytic mechanisms"
            ]
            properties["abstract"] = random.choice(abstracts)

        return properties