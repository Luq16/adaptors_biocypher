from __future__ import annotations

import random
import string
from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class OrthologyAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    ORTHOLOG_GROUP = auto()


class OrthologyAdapterOrthologGroupField(Enum):
    """
    Define possible fields the adapter can provide for ortholog groups.
    """
    ID = "id"
    NAME = "name"
    DESCRIPTION = "description"
    SIZE = "size"
    DATABASE_SOURCE = "database_source"
    FUNCTIONAL_CATEGORY = "functional_category"


class OrthologyAdapterEdgeType(Enum):
    """
    Enum for the types of edges the adapter can provide.
    """
    PROTEIN_IN_ORTHOLOG_GROUP = "protein_in_ortholog_group"
    PROTEIN_ORTHOLOGOUS_TO_PROTEIN = "protein_orthologous_to_protein"


class OrthologyAdapterProteinGroupEdgeField(Enum):
    """
    Define possible fields for protein-ortholog group edges.
    """
    ORGANISM = "organism"
    CONFIDENCE_SCORE = "confidence_score"
    SOURCE = "source"


class OrthologyAdapterOrthologyEdgeField(Enum):
    """
    Define possible fields for protein-protein orthology edges.
    """
    ORTHOLOGY_TYPE = "orthology_type"
    CONFIDENCE_SCORE = "confidence_score"
    SEQUENCE_IDENTITY = "sequence_identity"
    COVERAGE = "coverage"
    SOURCE = "source"
    METHOD = "method"


class OrthologyAdapter:
    """
    Orthology BioCypher adapter. Generates ortholog group nodes and protein-ortholog 
    association edges for creating a knowledge graph.

    Args:
        node_types: List of node types to include in the result.
        node_fields: List of node fields to include in the result.
        edge_types: List of edge types to include in the result.
        edge_fields: List of edge fields to include in the result.
        organisms: List of NCBI taxonomy IDs to include (default: human and mouse).
        test_mode: If True, limits the number of entries for testing.
    """

    def __init__(
        self,
        node_types: Optional[list] = None,
        node_fields: Optional[list] = None,
        edge_types: Optional[list] = None,
        edge_fields: Optional[list] = None,
        organisms: Optional[list] = None,
        test_mode: Optional[bool] = False,
    ):
        self.organisms = organisms or ["9606", "10090"]  # Human and mouse by default
        self.test_mode = test_mode
        self._set_types_and_fields(node_types, node_fields, edge_types, edge_fields)

    def get_nodes(self):
        """
        Returns a generator of node tuples for node types specified in the
        adapter constructor.
        """
        logger.info("Generating ortholog group nodes.")

        self.nodes = []

        if OrthologyAdapterNodeType.ORTHOLOG_GROUP in self.node_types:
            # In a real implementation, this would fetch from OrthoDB, OMA, PANTHER, etc.
            # For now, generating example ortholog groups
            num_groups = 15 if self.test_mode else 100
            [self.nodes.append(OrthologGroup(fields=self.node_fields)) 
             for _ in range(num_groups)]

        for node in self.nodes:
            yield (node.get_id(), node.get_label(), node.get_properties())

    def get_edges(self, target_proteins=None):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.

        Args:
            target_proteins: List of protein nodes to create orthology relationships with.
        """
        logger.info("Generating orthology association edges.")

        if not self.nodes:
            raise ValueError("No nodes found. Please run get_nodes() first.")

        # Create protein-ortholog group edges
        if (OrthologyAdapterEdgeType.PROTEIN_IN_ORTHOLOG_GROUP in self.edge_types 
            and target_proteins):
            
            for ortholog_group in self.nodes:
                # Each ortholog group contains 2-6 proteins from different organisms
                num_proteins = random.randint(2, 6)
                selected_proteins = random.sample(
                    target_proteins, 
                    min(num_proteins, len(target_proteins))
                )
                
                for protein in selected_proteins:
                    relationship_id = "".join(
                        random.choice(string.ascii_letters + string.digits)
                        for _ in range(10)
                    )
                    
                    properties = {}
                    for field in self.edge_fields:
                        if field == OrthologyAdapterProteinGroupEdgeField.ORGANISM:
                            properties["organism"] = random.choice(self.organisms)
                        elif field == OrthologyAdapterProteinGroupEdgeField.CONFIDENCE_SCORE:
                            properties["confidence_score"] = round(random.uniform(0.5, 1.0), 3)
                        elif field == OrthologyAdapterProteinGroupEdgeField.SOURCE:
                            properties["source"] = random.choice([
                                "OrthoDB", "OMA", "PANTHER", "EggNOG", "TreeFam", 
                                "OrthoMCL", "InParanoid", "Ensembl Compara"
                            ])

                    yield (
                        relationship_id,
                        protein.get_id() if hasattr(protein, 'get_id') else protein,
                        ortholog_group.get_id(),
                        OrthologyAdapterEdgeType.PROTEIN_IN_ORTHOLOG_GROUP.value,
                        properties,
                    )

        # Create direct protein-protein orthology edges
        if (OrthologyAdapterEdgeType.PROTEIN_ORTHOLOGOUS_TO_PROTEIN in self.edge_types 
            and target_proteins):
            
            # Create some direct orthology relationships between proteins
            proteins_list = list(target_proteins)
            num_orthology_pairs = min(50 if self.test_mode else 200, len(proteins_list) // 2)
            
            created_pairs = set()
            
            for _ in range(num_orthology_pairs):
                # Select two proteins
                protein_a, protein_b = random.sample(proteins_list, 2)
                
                # Create consistent pair identifier to avoid duplicates
                protein_a_id = protein_a.get_id() if hasattr(protein_a, 'get_id') else str(protein_a)
                protein_b_id = protein_b.get_id() if hasattr(protein_b, 'get_id') else str(protein_b)
                pair = tuple(sorted([protein_a_id, protein_b_id]))
                
                if pair not in created_pairs:
                    created_pairs.add(pair)
                    
                    relationship_id = "".join(
                        random.choice(string.ascii_letters + string.digits)
                        for _ in range(10)
                    )
                    
                    properties = {}
                    for field in self.edge_fields:
                        if field == OrthologyAdapterOrthologyEdgeField.ORTHOLOGY_TYPE:
                            properties["orthology_type"] = random.choice([
                                "ortholog", "co-ortholog", "in-paralog", "out-paralog"
                            ])
                        elif field == OrthologyAdapterOrthologyEdgeField.CONFIDENCE_SCORE:
                            properties["confidence_score"] = round(random.uniform(0.6, 1.0), 3)
                        elif field == OrthologyAdapterOrthologyEdgeField.SEQUENCE_IDENTITY:
                            properties["sequence_identity"] = round(random.uniform(0.3, 0.95), 3)
                        elif field == OrthologyAdapterOrthologyEdgeField.COVERAGE:
                            properties["coverage"] = round(random.uniform(0.5, 1.0), 3)
                        elif field == OrthologyAdapterOrthologyEdgeField.SOURCE:
                            properties["source"] = random.choice([
                                "OrthoDB", "OMA", "PANTHER", "EggNOG", "InParanoid", 
                                "Ensembl Compara", "BLAST", "OrthoMCL"
                            ])
                        elif field == OrthologyAdapterOrthologyEdgeField.METHOD:
                            properties["method"] = random.choice([
                                "reciprocal_best_hit", "phylogenetic_analysis", 
                                "synteny_analysis", "domain_architecture", 
                                "sequence_similarity", "gene_tree_reconciliation"
                            ])

                    yield (
                        relationship_id,
                        protein_a_id,
                        protein_b_id,
                        OrthologyAdapterEdgeType.PROTEIN_ORTHOLOGOUS_TO_PROTEIN.value,
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
            self.node_types = [type for type in OrthologyAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [field for field in OrthologyAdapterOrthologGroupField]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in OrthologyAdapterEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [
                field for field in chain(
                    OrthologyAdapterProteinGroupEdgeField,
                    OrthologyAdapterOrthologyEdgeField,
                )
            ]


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


class OrthologGroup(Node):
    """
    Generates instances of ortholog groups.
    """

    def __init__(self, fields: Optional[list] = None):
        self.fields = fields
        self.id = self._generate_id()
        self.label = "ortholog_group"
        self.properties = self._generate_properties()

    def _generate_id(self):
        """
        Generate a random ortholog group ID (OrthoDB-style).
        """
        return f"{random.randint(100000, 999999)}at{random.randint(2, 33208)}"

    def _generate_properties(self):
        properties = {}

        # Generate ortholog group name
        if OrthologyAdapterOrthologGroupField.NAME in self.fields:
            group_names = [
                "DNA repair protein family", "Histone family", "Ribosomal protein family",
                "Heat shock protein family", "Cytochrome oxidase family", 
                "Kinase family", "Transcription factor family", "Transporter family",
                "Enzyme family", "Structural protein family", "Signaling protein family",
                "Metabolic enzyme family", "Cell cycle protein family", 
                "Apoptosis protein family", "Immune system protein family",
                "Developmental protein family", "Neurotransmitter receptor family",
                "Ion channel family", "Growth factor family", "Hormone receptor family"
            ]
            properties["name"] = random.choice(group_names)

        # Generate description
        if OrthologyAdapterOrthologGroupField.DESCRIPTION in self.fields:
            group_name = properties.get("name", "protein family")
            properties["description"] = f"Orthologous group containing members of the {group_name.lower()}"

        # Generate group size
        if OrthologyAdapterOrthologGroupField.SIZE in self.fields:
            properties["size"] = random.randint(2, 20)

        # Generate database source
        if OrthologyAdapterOrthologGroupField.DATABASE_SOURCE in self.fields:
            properties["database_source"] = random.choice([
                "OrthoDB", "OMA", "PANTHER", "EggNOG", "TreeFam", "OrthoMCL"
            ])

        # Generate functional category
        if OrthologyAdapterOrthologGroupField.FUNCTIONAL_CATEGORY in self.fields:
            functional_categories = [
                "DNA repair", "transcription", "translation", "metabolism",
                "signal transduction", "transport", "cell cycle", "apoptosis",
                "immune response", "development", "protein folding", 
                "oxidative phosphorylation", "glycolysis", "protein degradation",
                "cell adhesion", "cytoskeleton", "vesicle transport", 
                "chromatin modification", "RNA processing", "protein synthesis"
            ]
            properties["functional_category"] = random.choice(functional_categories)

        return properties