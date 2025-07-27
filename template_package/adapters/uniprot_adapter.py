from __future__ import annotations

import random
import string
from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class UniProtAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    PROTEIN = auto()


class UniProtAdapterProteinField(Enum):
    """
    Define possible fields the adapter can provide for proteins.
    """
    ID = "id"
    NAME = "name" 
    SEQUENCE = "sequence"
    ORGANISM = "organism"
    GENE_NAME = "gene_name"
    FUNCTION = "function"
    LENGTH = "length"
    MASS = "mass"
    SUBCELLULAR_LOCATION = "subcellular_location"


class UniProtAdapterEdgeType(Enum):
    """
    Enum for the types of edges the adapter can provide.
    """
    PROTEIN_PROTEIN_INTERACTION = "protein_protein_interaction"


class UniProtAdapterProteinProteinEdgeField(Enum):
    """
    Define possible fields the adapter can provide for protein-protein edges.
    """
    INTERACTION_TYPE = "interaction_type"
    INTERACTION_SOURCE = "interaction_source"
    CONFIDENCE_SCORE = "confidence_score"


class UniProtAdapter:
    """
    UniProt BioCypher adapter. Generates protein nodes and protein-protein 
    interaction edges for creating a knowledge graph.

    Args:
        node_types: List of node types to include in the result.
        node_fields: List of node fields to include in the result.
        edge_types: List of edge types to include in the result.
        edge_fields: List of edge fields to include in the result.
        organism: NCBI taxonomy ID for organism filter (default: 9606 for human).
        test_mode: If True, limits the number of entries for testing.
    """

    def __init__(
        self,
        node_types: Optional[list] = None,
        node_fields: Optional[list] = None,
        edge_types: Optional[list] = None,
        edge_fields: Optional[list] = None,
        organism: Optional[str] = "9606",
        test_mode: Optional[bool] = False,
    ):
        self.organism = organism
        self.test_mode = test_mode
        self._set_types_and_fields(node_types, node_fields, edge_types, edge_fields)

    def get_nodes(self):
        """
        Returns a generator of node tuples for node types specified in the
        adapter constructor.
        """
        logger.info("Generating UniProt protein nodes.")

        self.nodes = []

        if UniProtAdapterNodeType.PROTEIN in self.node_types:
            # In a real implementation, this would fetch from UniProt API/database
            # For now, generating example proteins
            num_proteins = 50 if self.test_mode else 1000
            [self.nodes.append(Protein(fields=self.node_fields, organism=self.organism)) 
             for _ in range(num_proteins)]

        for node in self.nodes:
            yield (node.get_id(), node.get_label(), node.get_properties())

    def get_edges(self, probability: float = 0.1):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.

        Args:
            probability: Probability of generating an edge between two proteins.
        """
        logger.info("Generating protein-protein interaction edges.")

        if not self.nodes:
            raise ValueError("No nodes found. Please run get_nodes() first.")

        for node in self.nodes:
            if random.random() < probability:
                other_node = random.choice(self.nodes)
                
                # Only create edges between different proteins
                if node.get_id() != other_node.get_id():
                    # Generate random relationship id
                    relationship_id = "".join(
                        random.choice(string.ascii_letters + string.digits)
                        for _ in range(10)
                    )

                    if (
                        isinstance(other_node, Protein)
                        and UniProtAdapterEdgeType.PROTEIN_PROTEIN_INTERACTION
                        in self.edge_types
                    ):
                        edge_type = UniProtAdapterEdgeType.PROTEIN_PROTEIN_INTERACTION.value
                        
                        # Create edge properties
                        properties = {}
                        for field in self.edge_fields:
                            if field == UniProtAdapterProteinProteinEdgeField.INTERACTION_TYPE:
                                properties["interaction_type"] = random.choice([
                                    "physical association", "direct interaction", 
                                    "association", "colocalization"
                                ])
                            elif field == UniProtAdapterProteinProteinEdgeField.INTERACTION_SOURCE:
                                properties["interaction_source"] = random.choice([
                                    "experimental", "database", "text mining"
                                ])
                            elif field == UniProtAdapterProteinProteinEdgeField.CONFIDENCE_SCORE:
                                properties["confidence_score"] = round(random.uniform(0.1, 1.0), 3)

                        yield (
                            relationship_id,
                            node.get_id(),
                            other_node.get_id(),
                            edge_type,
                            properties,
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
            self.node_types = [type for type in UniProtAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [field for field in UniProtAdapterProteinField]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in UniProtAdapterEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in UniProtAdapterProteinProteinEdgeField]


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


class Protein(Node):
    """
    Generates instances of proteins from UniProt.
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
        if UniProtAdapterProteinField.ORGANISM in self.fields:
            properties["organism"] = self.organism

        # Generate protein name
        if UniProtAdapterProteinField.NAME in self.fields:
            protein_names = [
                "Tumor protein p53", "Epidermal growth factor receptor",
                "Insulin receptor", "Hemoglobin subunit alpha",
                "Cytochrome c oxidase subunit 1", "Heat shock protein 90",
                "Ubiquitin", "Histone H3", "Beta-actin", "Albumin"
            ]
            properties["name"] = random.choice(protein_names)

        # Generate gene name
        if UniProtAdapterProteinField.GENE_NAME in self.fields:
            gene_names = ["TP53", "EGFR", "INSR", "HBA1", "COX1", "HSP90AA1", "UBB", "H3C1", "ACTB", "ALB"]
            properties["gene_name"] = random.choice(gene_names)

        # Generate sequence
        if UniProtAdapterProteinField.SEQUENCE in self.fields:
            length = random.randint(100, 1000)
            amino_acids = "ACDEFGHIKLMNPQRSTVWY"
            properties["sequence"] = "".join(
                random.choices(amino_acids, k=length)
            )

        # Generate length
        if UniProtAdapterProteinField.LENGTH in self.fields:
            properties["length"] = random.randint(100, 1000)

        # Generate mass
        if UniProtAdapterProteinField.MASS in self.fields:
            properties["mass"] = random.randint(10000, 200000)

        # Generate function
        if UniProtAdapterProteinField.FUNCTION in self.fields:
            functions = [
                "Transcription factor", "Kinase", "Phosphatase", "Transporter",
                "Receptor", "Enzyme", "Structural protein", "Signaling protein"
            ]
            properties["function"] = random.choice(functions)

        # Generate subcellular location
        if UniProtAdapterProteinField.SUBCELLULAR_LOCATION in self.fields:
            locations = [
                "Nucleus", "Cytoplasm", "Mitochondria", "Endoplasmic reticulum",
                "Golgi apparatus", "Cell membrane", "Extracellular"
            ]
            properties["subcellular_location"] = random.choice(locations)

        return properties