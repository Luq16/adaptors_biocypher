from __future__ import annotations

import random
import string
from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class PPIAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    PROTEIN = auto()


class PPIAdapterProteinField(Enum):
    """
    Define possible fields the adapter can provide for proteins.
    """
    ID = "id"
    NAME = "name"
    ORGANISM = "organism"
    GENE_NAME = "gene_name"


class PPIAdapterEdgeType(Enum):
    """
    Enum for the types of edges the adapter can provide.
    """
    PROTEIN_PROTEIN_INTERACTION = "protein_protein_interaction"


class PPIAdapterInteractionEdgeField(Enum):
    """
    Define possible fields for protein-protein interaction edges.
    """
    INTERACTION_TYPE = "interaction_type"
    DETECTION_METHOD = "detection_method"
    CONFIDENCE_SCORE = "confidence_score"
    SOURCE_DATABASE = "source_database"
    PUBMED_IDS = "pubmed_ids"
    EXPERIMENTAL_SYSTEM = "experimental_system"
    INTERACTION_IDENTIFIER = "interaction_identifier"
    THROUGHPUT = "throughput"


class PPIAdapter:
    """
    Protein-Protein Interaction BioCypher adapter. Generates protein nodes 
    and protein-protein interaction edges for creating a knowledge graph.

    Args:
        node_types: List of node types to include in the result.
        node_fields: List of node fields to include in the result.
        edge_types: List of edge types to include in the result.
        edge_fields: List of edge fields to include in the result.
        organism: NCBI taxonomy ID for organism filter (default: 9606 for human).
        test_mode: If True, limits the number of entries for testing.
        interaction_probability: Probability of creating interaction between proteins.
    """

    def __init__(
        self,
        node_types: Optional[list] = None,
        node_fields: Optional[list] = None,
        edge_types: Optional[list] = None,
        edge_fields: Optional[list] = None,
        organism: Optional[str] = "9606",
        test_mode: Optional[bool] = False,
        interaction_probability: Optional[float] = 0.05,
    ):
        self.organism = organism
        self.test_mode = test_mode
        self.interaction_probability = interaction_probability
        self._set_types_and_fields(node_types, node_fields, edge_types, edge_fields)

    def get_nodes(self):
        """
        Returns a generator of node tuples for node types specified in the
        adapter constructor. 
        
        Note: In a real implementation, this might not be needed if proteins
        are already provided by a UniProt adapter.
        """
        logger.info("Generating protein nodes for PPI.")

        self.nodes = []

        if PPIAdapterNodeType.PROTEIN in self.node_types:
            # Generate proteins for PPI network
            num_proteins = 50 if self.test_mode else 500
            [self.nodes.append(ProteinForPPI(fields=self.node_fields, organism=self.organism)) 
             for _ in range(num_proteins)]

        for node in self.nodes:
            yield (node.get_id(), node.get_label(), node.get_properties())

    def get_edges(self, protein_nodes=None):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.

        Args:
            protein_nodes: List of protein nodes to create interactions between.
                          If None, uses nodes from get_nodes().
        """
        logger.info("Generating protein-protein interaction edges.")

        if protein_nodes is None:
            if not self.nodes:
                raise ValueError("No protein nodes found. Please run get_nodes() first or provide protein_nodes.")
            protein_nodes = self.nodes
        
        # Convert to list if it's a generator
        if not isinstance(protein_nodes, list):
            protein_nodes = list(protein_nodes)

        if PPIAdapterEdgeType.PROTEIN_PROTEIN_INTERACTION in self.edge_types:
            interactions_created = set()  # To avoid duplicate interactions
            
            for i, protein_a in enumerate(protein_nodes):
                for protein_b in protein_nodes[i+1:]:  # Avoid self-interactions and duplicates
                    
                    if random.random() < self.interaction_probability:
                        # Create unique interaction identifier
                        protein_a_id = protein_a.get_id() if hasattr(protein_a, 'get_id') else str(protein_a)
                        protein_b_id = protein_b.get_id() if hasattr(protein_b, 'get_id') else str(protein_b)
                        
                        # Ensure consistent ordering for interaction pairs
                        pair = tuple(sorted([protein_a_id, protein_b_id]))
                        
                        if pair not in interactions_created:
                            interactions_created.add(pair)
                            
                            relationship_id = "".join(
                                random.choice(string.ascii_letters + string.digits)
                                for _ in range(12)
                            )
                            
                            properties = {}
                            for field in self.edge_fields:
                                if field == PPIAdapterInteractionEdgeField.INTERACTION_TYPE:
                                    properties["interaction_type"] = random.choice([
                                        "physical association", "direct interaction", 
                                        "association", "colocalization", "genetic interaction"
                                    ])
                                elif field == PPIAdapterInteractionEdgeField.DETECTION_METHOD:
                                    properties["detection_method"] = random.choice([
                                        "yeast two-hybrid", "affinity chromatography",
                                        "co-immunoprecipitation", "pull down",
                                        "fluorescence microscopy", "mass spectrometry",
                                        "cross-linking", "surface plasmon resonance"
                                    ])
                                elif field == PPIAdapterInteractionEdgeField.CONFIDENCE_SCORE:
                                    properties["confidence_score"] = round(random.uniform(0.1, 1.0), 3)
                                elif field == PPIAdapterInteractionEdgeField.SOURCE_DATABASE:
                                    properties["source_database"] = random.choice([
                                        "BioGRID", "IntAct", "MINT", "DIP", "HPRD", 
                                        "STRING", "APID", "MIPS", "Reactome"
                                    ])
                                elif field == PPIAdapterInteractionEdgeField.PUBMED_IDS:
                                    # Generate 1-3 PubMed IDs
                                    num_pubmed = random.randint(1, 3)
                                    pubmed_ids = [str(random.randint(10000000, 35000000)) for _ in range(num_pubmed)]
                                    properties["pubmed_ids"] = "|".join(pubmed_ids)
                                elif field == PPIAdapterInteractionEdgeField.EXPERIMENTAL_SYSTEM:
                                    properties["experimental_system"] = random.choice([
                                        "Yeast Two-Hybrid", "Affinity Capture-MS",
                                        "Co-fractionation", "Reconstituted Complex",
                                        "Co-crystal Structure", "FRET", "PCA",
                                        "Far Western", "Co-localization"
                                    ])
                                elif field == PPIAdapterInteractionEdgeField.INTERACTION_IDENTIFIER:
                                    properties["interaction_identifier"] = f"EBI-{random.randint(1000000, 9999999)}"
                                elif field == PPIAdapterInteractionEdgeField.THROUGHPUT:
                                    properties["throughput"] = random.choice([
                                        "High Throughput", "Low Throughput"
                                    ])

                            yield (
                                relationship_id,
                                protein_a_id,
                                protein_b_id,
                                PPIAdapterEdgeType.PROTEIN_PROTEIN_INTERACTION.value,
                                properties,
                            )

    def get_node_count(self):
        """
        Returns the number of nodes generated by the adapter.
        """
        return len(list(self.get_nodes()))

    def get_interaction_count(self, protein_nodes=None):
        """
        Returns the estimated number of interactions that would be generated.
        """
        if protein_nodes is None:
            if not self.nodes:
                self.get_nodes()  # Generate nodes if not already done
            protein_nodes = self.nodes
        
        n_proteins = len(list(protein_nodes))
        max_interactions = (n_proteins * (n_proteins - 1)) // 2
        estimated_interactions = int(max_interactions * self.interaction_probability)
        return estimated_interactions

    def _set_types_and_fields(self, node_types, node_fields, edge_types, edge_fields):
        if node_types:
            self.node_types = node_types
        else:
            self.node_types = [type for type in PPIAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [field for field in PPIAdapterProteinField]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in PPIAdapterEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in PPIAdapterInteractionEdgeField]


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


class ProteinForPPI(Node):
    """
    Generates instances of proteins for PPI networks.
    """

    def __init__(self, fields: Optional[list] = None, organism: str = "9606"):
        self.fields = fields
        self.organism = organism
        self.id = self._generate_id()
        self.label = "protein"
        self.properties = self._generate_properties()

    def _generate_id(self):
        """
        Generate a random UniProt-style accession ID.
        """
        # UniProt format: 1-2 letters + 4-5 alphanumeric characters
        prefix = "".join(random.choices(string.ascii_uppercase, k=random.choice([1, 2])))
        suffix = "".join(random.choices(string.ascii_uppercase + string.digits, k=random.choice([4, 5])))
        return f"{prefix}{suffix}"

    def _generate_properties(self):
        properties = {}

        # Add organism information
        if PPIAdapterProteinField.ORGANISM in self.fields:
            properties["organism"] = self.organism

        # Generate protein name
        if PPIAdapterProteinField.NAME in self.fields:
            protein_names = [
                "Tumor protein p53", "Epidermal growth factor receptor",
                "Insulin receptor", "Hemoglobin subunit alpha",
                "Cytochrome c oxidase subunit 1", "Heat shock protein 90",
                "Ubiquitin", "Histone H3", "Beta-actin", "Albumin",
                "BRCA1", "BRCA2", "ATM kinase", "p21", "MDM2",
                "EGFR", "HER2", "VEGFR", "PDGFR", "c-MYC",
                "p16", "Rb protein", "E2F1", "Cyclin D1", "CDK4"
            ]
            properties["name"] = random.choice(protein_names)

        # Generate gene name
        if PPIAdapterProteinField.GENE_NAME in self.fields:
            gene_names = [
                "TP53", "EGFR", "INSR", "HBA1", "COX1", "HSP90AA1", 
                "UBB", "H3C1", "ACTB", "ALB", "BRCA1", "BRCA2", 
                "ATM", "CDKN1A", "MDM2", "ERBB2", "VEGFR1", "PDGFRA", 
                "MYC", "CDKN2A", "RB1", "E2F1", "CCND1", "CDK4"
            ]
            properties["gene_name"] = random.choice(gene_names)

        return properties