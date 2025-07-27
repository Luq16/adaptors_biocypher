from __future__ import annotations

import random
import string
from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class IntactAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    PROTEIN = auto()


class IntactAdapterProteinField(Enum):
    """
    Define possible fields the adapter can provide for proteins.
    """
    NAME = "name"
    ORGANISM = "organism"
    GENE_NAME = "gene_name"
    INTACT_AC = "intact_ac"
    UNIPROT_ID = "uniprot_id"
    TAXONOMY_ID = "taxonomy_id"


class IntactAdapterEdgeType(Enum):
    """
    Enum for the types of edges the adapter can provide.
    """
    PROTEIN_PROTEIN_INTERACTION = "protein_protein_interaction"


class IntactAdapterInteractionEdgeField(Enum):
    """
    Define possible fields for protein-protein interaction edges.
    """
    SOURCE = "source"
    INTERACTION_AC = "interaction_ac"
    MI_SCORE = "mi_score"
    INTERACTION_TYPE = "interaction_type"
    DETECTION_METHOD = "detection_method"
    PUBMED_IDS = "pubmed_ids"
    INTERACTION_XREFS = "interaction_xrefs"
    CONFIDENCE_VALUE = "confidence_value"
    EXPANSION_METHOD = "expansion_method"
    BIOLOGICAL_ROLE_A = "biological_role_a"
    BIOLOGICAL_ROLE_B = "biological_role_b"
    EXPERIMENTAL_ROLE_A = "experimental_role_a"
    EXPERIMENTAL_ROLE_B = "experimental_role_b"
    INTERACTOR_TYPE_A = "interactor_type_a"
    INTERACTOR_TYPE_B = "interactor_type_b"
    SOURCE_DATABASE = "source_database"
    INTERACTION_IDENTIFIER = "interaction_identifier"
    COMPLEX_EXPANSION = "complex_expansion"
    DATASET = "dataset"


class IntactAdapter:
    """
    IntAct BioCypher adapter. Generates protein nodes and high-quality curated 
    protein-protein interaction edges from IntAct database.

    Args:
        node_types: List of node types to include in the result.
        node_fields: List of node fields to include in the result.
        edge_types: List of edge types to include in the result.
        edge_fields: List of edge fields to include in the result.
        organism: NCBI taxonomy ID for organism filter (default: 9606 for human).
        test_mode: If True, limits the number of entries for testing.
        min_mi_score: Minimum MI score threshold for interactions.
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
        min_mi_score: Optional[float] = 0.4,
        interaction_probability: Optional[float] = 0.06,
    ):
        self.organism = organism
        self.test_mode = test_mode
        self.min_mi_score = min_mi_score
        self.interaction_probability = interaction_probability
        self._set_types_and_fields(node_types, node_fields, edge_types, edge_fields)

    def get_nodes(self):
        """
        Returns a generator of node tuples for node types specified in the
        adapter constructor.
        """
        logger.info("Generating IntAct protein nodes.")

        self.nodes = []

        if IntactAdapterNodeType.PROTEIN in self.node_types:
            # Generate proteins for IntAct network
            num_proteins = 40 if self.test_mode else 400
            [self.nodes.append(IntactProtein(fields=self.node_fields, organism=self.organism)) 
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
        logger.info("Generating IntAct protein-protein interaction edges.")

        if protein_nodes is None:
            if not self.nodes:
                raise ValueError("No protein nodes found. Please run get_nodes() first or provide protein_nodes.")
            protein_nodes = self.nodes
        
        # Convert to list if it's a generator
        if not isinstance(protein_nodes, list):
            protein_nodes = list(protein_nodes)

        if IntactAdapterEdgeType.PROTEIN_PROTEIN_INTERACTION in self.edge_types:
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
                            
                            # Generate MI score (molecular interaction confidence)
                            mi_score = round(random.uniform(self.min_mi_score, 1.0), 3)
                            
                            relationship_id = "".join(
                                random.choice(string.ascii_letters + string.digits)
                                for _ in range(12)
                            )
                            
                            properties = {}
                            for field in self.edge_fields:
                                if field == IntactAdapterInteractionEdgeField.SOURCE:
                                    properties["source"] = "IntAct"
                                elif field == IntactAdapterInteractionEdgeField.INTERACTION_AC:
                                    properties["interaction_ac"] = f"EBI-{random.randint(1000000, 9999999)}"
                                elif field == IntactAdapterInteractionEdgeField.MI_SCORE:
                                    properties["mi_score"] = mi_score
                                elif field == IntactAdapterInteractionEdgeField.INTERACTION_TYPE:
                                    properties["interaction_type"] = random.choice([
                                        "direct interaction", "physical association",
                                        "association", "colocalization", "genetic interaction",
                                        "enzymatic reaction", "cleavage reaction"
                                    ])
                                elif field == IntactAdapterInteractionEdgeField.DETECTION_METHOD:
                                    properties["detection_method"] = random.choice([
                                        "two hybrid", "anti tag coimmunoprecipitation",
                                        "pull down", "affinity chromatography technology",
                                        "anti bait coimmunoprecipitation", "fluorescence microscopy",
                                        "x-ray crystallography", "nuclear magnetic resonance",
                                        "surface plasmon resonance", "isothermal titration calorimetry"
                                    ])
                                elif field == IntactAdapterInteractionEdgeField.PUBMED_IDS:
                                    # IntAct typically has well-curated references
                                    num_pubmed = random.randint(1, 3)
                                    pubmed_ids = [str(random.randint(12000000, 35000000)) for _ in range(num_pubmed)]
                                    properties["pubmed_ids"] = "|".join(pubmed_ids)
                                elif field == IntactAdapterInteractionEdgeField.INTERACTION_XREFS:
                                    # Cross-references to other databases
                                    xrefs = [
                                        f"mint:MINT-{random.randint(1000000, 9999999)}",
                                        f"dip:DIP-{random.randint(10000, 99999)}N"
                                    ]
                                    properties["interaction_xrefs"] = "|".join(random.sample(xrefs, random.randint(1, 2)))
                                elif field == IntactAdapterInteractionEdgeField.CONFIDENCE_VALUE:
                                    properties["confidence_value"] = f"intact-miscore:{mi_score}"
                                elif field == IntactAdapterInteractionEdgeField.EXPANSION_METHOD:
                                    properties["expansion_method"] = random.choice([
                                        "spoke expansion", "matrix expansion", "bipartite expansion"
                                    ])
                                elif field == IntactAdapterInteractionEdgeField.BIOLOGICAL_ROLE_A:
                                    properties["biological_role_a"] = random.choice([
                                        "unspecified role", "enzyme", "enzyme target",
                                        "transcription factor", "cofactor", "electron acceptor"
                                    ])
                                elif field == IntactAdapterInteractionEdgeField.BIOLOGICAL_ROLE_B:
                                    properties["biological_role_b"] = random.choice([
                                        "unspecified role", "enzyme", "enzyme target",
                                        "transcription factor", "cofactor", "electron donor"
                                    ])
                                elif field == IntactAdapterInteractionEdgeField.EXPERIMENTAL_ROLE_A:
                                    properties["experimental_role_a"] = random.choice([
                                        "bait", "prey", "neutral component", "self"
                                    ])
                                elif field == IntactAdapterInteractionEdgeField.EXPERIMENTAL_ROLE_B:
                                    properties["experimental_role_b"] = random.choice([
                                        "bait", "prey", "neutral component", "self"
                                    ])
                                elif field == IntactAdapterInteractionEdgeField.INTERACTOR_TYPE_A:
                                    properties["interactor_type_a"] = "protein"
                                elif field == IntactAdapterInteractionEdgeField.INTERACTOR_TYPE_B:
                                    properties["interactor_type_b"] = "protein"
                                elif field == IntactAdapterInteractionEdgeField.SOURCE_DATABASE:
                                    properties["source_database"] = random.choice([
                                        "intact", "mint", "dip", "biogrid", "string"
                                    ])
                                elif field == IntactAdapterInteractionEdgeField.INTERACTION_IDENTIFIER:
                                    properties["interaction_identifier"] = f"EBI-{random.randint(10000000, 99999999)}"
                                elif field == IntactAdapterInteractionEdgeField.COMPLEX_EXPANSION:
                                    if random.random() < 0.3:  # 30% chance
                                        properties["complex_expansion"] = random.choice([
                                            "spoke", "matrix", "bipartite"
                                        ])
                                elif field == IntactAdapterInteractionEdgeField.DATASET:
                                    properties["dataset"] = random.choice([
                                        "Krogan - Coimmunoprecipitation",
                                        "Ho - Immunoprecipitation",
                                        "Gavin - Immunoprecipitation", 
                                        "Collins - Immunoprecipitation",
                                        "Boutros - Y2H",
                                        "Vidal - Y2H"
                                    ])

                            yield (
                                relationship_id,
                                protein_a_id,
                                protein_b_id,
                                IntactAdapterEdgeType.PROTEIN_PROTEIN_INTERACTION.value,
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
            self.node_types = [type for type in IntactAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [field for field in IntactAdapterProteinField]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in IntactAdapterEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in IntactAdapterInteractionEdgeField]


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


class IntactProtein(Node):
    """
    Generates instances of proteins for IntAct networks.
    """

    def __init__(self, fields: Optional[list] = None, organism: str = "9606"):
        self.fields = fields
        self.organism = organism
        self.id = self._generate_id()
        self.label = "protein"
        self.properties = self._generate_properties()

    def _generate_id(self):
        """
        Generate an IntAct-style protein identifier.
        """
        # IntAct format: EBI-xxxxxxx
        return f"EBI-{random.randint(1000000, 9999999)}"

    def _generate_properties(self):
        properties = {}

        # Add organism information
        if IntactAdapterProteinField.ORGANISM in self.fields:
            properties["organism"] = self.organism

        # Add taxonomy ID
        if IntactAdapterProteinField.TAXONOMY_ID in self.fields:
            properties["taxonomy_id"] = self.organism

        # Generate protein name
        if IntactAdapterProteinField.NAME in self.fields:
            protein_names = [
                "Cellular tumor antigen p53",
                "Epidermal growth factor receptor",
                "Insulin receptor",
                "Catenin beta-1",
                "Proto-oncogene tyrosine-protein kinase Src",
                "Serine/threonine-protein kinase Akt",
                "Mitogen-activated protein kinase 1",
                "Transcription factor AP-1",
                "Nuclear factor NF-kappa-B p65 subunit",
                "Cyclin-dependent kinase 2",
                "Retinoblastoma-associated protein",
                "DNA topoisomerase 2-alpha",
                "Heat shock protein HSP 90-alpha",
                "Ubiquitin-protein ligase E3A",
                "26S proteasome regulatory subunit 6A",
                "Histone deacetylase 1",
                "Chromatin modification complex subunit",
                "Ribosomal protein S6 kinase alpha-1",
                "AMP-activated protein kinase catalytic subunit alpha-1",
                "Phosphatidylinositol 3-kinase regulatory subunit alpha"
            ]
            properties["name"] = random.choice(protein_names)

        # Generate gene name
        if IntactAdapterProteinField.GENE_NAME in self.fields:
            gene_names = [
                "TP53", "EGFR", "INSR", "CTNNB1", "SRC", "AKT1", "MAPK1", 
                "JUN", "RELA", "CDK2", "RB1", "TOP2A", "HSP90AA1", "UBE3A",
                "PSMC3", "HDAC1", "SMARCA4", "RPS6KA1", "PRKAA1", "PIK3R1",
                "BRCA1", "BRCA2", "ATM", "PTEN", "MYC", "RAS", "RAF1", "MEK1"
            ]
            properties["gene_name"] = random.choice(gene_names)

        # Generate IntAct AC
        if IntactAdapterProteinField.INTACT_AC in self.fields:
            properties["intact_ac"] = self.id

        # Generate UniProt ID
        if IntactAdapterProteinField.UNIPROT_ID in self.fields:
            # Generate UniProt-style ID
            prefix = "".join(random.choices(string.ascii_uppercase, k=random.choice([1, 2])))
            suffix = "".join(random.choices(string.ascii_uppercase + string.digits, k=random.choice([4, 5])))
            properties["uniprot_id"] = f"{prefix}{suffix}"

        return properties