from __future__ import annotations

import random
import string
from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class StringAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    PROTEIN = auto()


class StringAdapterProteinField(Enum):
    """
    Define possible fields the adapter can provide for proteins.
    """
    NAME = "name"
    ORGANISM = "organism"
    GENE_NAME = "gene_name"
    STRING_ID = "string_id"
    UNIPROT_ID = "uniprot_id"
    ENSEMBL_ID = "ensembl_id"


class StringAdapterEdgeType(Enum):
    """
    Enum for the types of edges the adapter can provide.
    """
    PROTEIN_PROTEIN_INTERACTION = "protein_protein_interaction"
    FUNCTIONAL_ASSOCIATION = "functional_association"


class StringAdapterInteractionEdgeField(Enum):
    """
    Define possible fields for protein-protein interaction edges.
    """
    SOURCE = "source"
    COMBINED_SCORE = "combined_score"
    NEIGHBORHOOD_SCORE = "neighborhood_score"
    FUSION_SCORE = "fusion_score"
    COOCCURRENCE_SCORE = "cooccurrence_score"
    COEXPRESSION_SCORE = "coexpression_score"
    EXPERIMENTAL_SCORE = "experimental_score"
    DATABASE_SCORE = "database_score"
    TEXTMINING_SCORE = "textmining_score"
    PHYSICAL_INTERACTION_SCORE = "physical_interaction_score"
    FUNCTIONAL_ASSOCIATION_SCORE = "functional_association_score"
    CONFIDENCE_LEVEL = "confidence_level"
    ANNOTATION_TRANSFER = "annotation_transfer"


class StringAdapter:
    """
    STRING BioCypher adapter. Generates protein nodes and functional protein 
    association edges from STRING database with confidence scoring.

    Args:
        node_types: List of node types to include in the result.
        node_fields: List of node fields to include in the result.
        edge_types: List of edge types to include in the result.
        edge_fields: List of edge fields to include in the result.
        organism: NCBI taxonomy ID for organism filter (default: 9606 for human).
        test_mode: If True, limits the number of entries for testing.
        min_combined_score: Minimum STRING combined score threshold (0-1000).
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
        min_combined_score: Optional[int] = 400,
        interaction_probability: Optional[float] = 0.12,
    ):
        self.organism = organism
        self.test_mode = test_mode
        self.min_combined_score = min_combined_score
        self.interaction_probability = interaction_probability
        self._set_types_and_fields(node_types, node_fields, edge_types, edge_fields)

    def get_nodes(self):
        """
        Returns a generator of node tuples for node types specified in the
        adapter constructor.
        """
        logger.info("Generating STRING protein nodes.")

        self.nodes = []

        if StringAdapterNodeType.PROTEIN in self.node_types:
            # Generate proteins for STRING network
            num_proteins = 70 if self.test_mode else 700
            [self.nodes.append(StringProtein(fields=self.node_fields, organism=self.organism)) 
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
        logger.info("Generating STRING protein association edges.")

        if protein_nodes is None:
            if not self.nodes:
                raise ValueError("No protein nodes found. Please run get_nodes() first or provide protein_nodes.")
            protein_nodes = self.nodes
        
        # Convert to list if it's a generator
        if not isinstance(protein_nodes, list):
            protein_nodes = list(protein_nodes)

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
                        
                        # Generate STRING scores
                        combined_score = random.randint(self.min_combined_score, 999)
                        
                        # Determine edge type based on combined score
                        if combined_score > 700:
                            edge_type = StringAdapterEdgeType.PROTEIN_PROTEIN_INTERACTION.value
                        else:
                            edge_type = StringAdapterEdgeType.FUNCTIONAL_ASSOCIATION.value
                        
                        # Check if this edge type should be included
                        edge_type_enum = None
                        for et in self.edge_types:
                            if et.value == edge_type:
                                edge_type_enum = et
                                break
                        
                        if edge_type_enum is None:
                            continue
                        
                        relationship_id = "".join(
                            random.choice(string.ascii_letters + string.digits)
                            for _ in range(12)
                        )
                        
                        properties = {}
                        for field in self.edge_fields:
                            if field == StringAdapterInteractionEdgeField.SOURCE:
                                properties["source"] = "STRING"
                            elif field == StringAdapterInteractionEdgeField.COMBINED_SCORE:
                                properties["combined_score"] = combined_score
                            elif field == StringAdapterInteractionEdgeField.NEIGHBORHOOD_SCORE:
                                properties["neighborhood_score"] = random.randint(0, 300)
                            elif field == StringAdapterInteractionEdgeField.FUSION_SCORE:
                                properties["fusion_score"] = random.randint(0, 200)
                            elif field == StringAdapterInteractionEdgeField.COOCCURRENCE_SCORE:
                                properties["cooccurrence_score"] = random.randint(0, 400)
                            elif field == StringAdapterInteractionEdgeField.COEXPRESSION_SCORE:
                                properties["coexpression_score"] = random.randint(0, 500)
                            elif field == StringAdapterInteractionEdgeField.EXPERIMENTAL_SCORE:
                                properties["experimental_score"] = random.randint(0, 600)
                            elif field == StringAdapterInteractionEdgeField.DATABASE_SCORE:
                                properties["database_score"] = random.randint(0, 400)
                            elif field == StringAdapterInteractionEdgeField.TEXTMINING_SCORE:
                                properties["textmining_score"] = random.randint(0, 800)
                            elif field == StringAdapterInteractionEdgeField.PHYSICAL_INTERACTION_SCORE:
                                # Physical score for high-confidence interactions
                                if combined_score > 700:
                                    properties["physical_interaction_score"] = random.randint(300, 900)
                                else:
                                    properties["physical_interaction_score"] = random.randint(0, 300)
                            elif field == StringAdapterInteractionEdgeField.FUNCTIONAL_ASSOCIATION_SCORE:
                                properties["functional_association_score"] = random.randint(200, 800)
                            elif field == StringAdapterInteractionEdgeField.CONFIDENCE_LEVEL:
                                if combined_score > 900:
                                    properties["confidence_level"] = "highest"
                                elif combined_score > 700:
                                    properties["confidence_level"] = "high"
                                elif combined_score > 400:
                                    properties["confidence_level"] = "medium"
                                else:
                                    properties["confidence_level"] = "low"
                            elif field == StringAdapterInteractionEdgeField.ANNOTATION_TRANSFER:
                                properties["annotation_transfer"] = random.choice([
                                    "true", "false"
                                ])

                        yield (
                            relationship_id,
                            protein_a_id,
                            protein_b_id,
                            edge_type,
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
            self.node_types = [type for type in StringAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [field for field in StringAdapterProteinField]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in StringAdapterEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in StringAdapterInteractionEdgeField]


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


class StringProtein(Node):
    """
    Generates instances of proteins for STRING networks.
    """

    def __init__(self, fields: Optional[list] = None, organism: str = "9606"):
        self.fields = fields
        self.organism = organism
        self.id = self._generate_id()
        self.label = "protein"
        self.properties = self._generate_properties()

    def _generate_id(self):
        """
        Generate a STRING-style protein identifier.
        """
        # STRING format: organism.ENSPxxxxxxxxxxxxxxx or similar
        ensembl_id = f"ENSP{random.randint(10000000000000, 99999999999999):014d}"
        return f"{self.organism}.{ensembl_id}"

    def _generate_properties(self):
        properties = {}

        # Add organism information
        if StringAdapterProteinField.ORGANISM in self.fields:
            properties["organism"] = self.organism

        # Generate protein name
        if StringAdapterProteinField.NAME in self.fields:
            protein_names = [
                "60S ribosomal protein L10", "40S ribosomal protein S3",
                "Glyceraldehyde-3-phosphate dehydrogenase", "Beta-actin",
                "Alpha-tubulin", "Beta-tubulin", "Elongation factor 1-alpha 1",
                "Heat shock protein 90", "Heat shock cognate 71 kDa protein",
                "Pyruvate kinase PKM", "Lactate dehydrogenase A",
                "Aldolase A", "Phosphoglycerate kinase 1", "Enolase 1",
                "Triose phosphate isomerase 1", "Phosphoglycerate mutase 1",
                "Fructose-bisphosphate aldolase A", "Glucose-6-phosphate isomerase",
                "Hexokinase-1", "6-phosphofructokinase", "Ubiquitin-60S ribosomal protein L40",
                "Histone H1.0", "Histone H2A", "Histone H2B", "Histone H3",
                "Histone H4", "40S ribosomal protein S6", "60S ribosomal protein L7a"
            ]
            properties["name"] = random.choice(protein_names)

        # Generate gene name
        if StringAdapterProteinField.GENE_NAME in self.fields:
            gene_names = [
                "RPL10", "RPS3", "GAPDH", "ACTB", "TUBA1A", "TUBB", "EEF1A1",
                "HSP90AA1", "HSPA8", "PKM", "LDHA", "ALDOA", "PGK1", "ENO1",
                "TPI1", "PGAM1", "ALDOA", "GPI", "HK1", "PFKM", "UBA52",
                "H1-0", "H2AX", "H2BC1", "H3C1", "H4C1", "RPS6", "RPL7A",
                "PPIA", "YWHAZ", "CALM1", "PSMA1", "PSMB1", "PSMC1"
            ]
            properties["gene_name"] = random.choice(gene_names)

        # Generate STRING ID
        if StringAdapterProteinField.STRING_ID in self.fields:
            properties["string_id"] = self.id

        # Generate UniProt ID
        if StringAdapterProteinField.UNIPROT_ID in self.fields:
            # Generate UniProt-style ID
            prefix = "".join(random.choices(string.ascii_uppercase, k=random.choice([1, 2])))
            suffix = "".join(random.choices(string.ascii_uppercase + string.digits, k=random.choice([4, 5])))
            properties["uniprot_id"] = f"{prefix}{suffix}"

        # Generate Ensembl ID
        if StringAdapterProteinField.ENSEMBL_ID in self.fields:
            ensembl_protein_id = f"ENSP{random.randint(10000000000000, 99999999999999):014d}"
            properties["ensembl_id"] = ensembl_protein_id

        return properties