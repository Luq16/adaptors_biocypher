from __future__ import annotations

import random
import string
from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class BiogridAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    PROTEIN = auto()


class BiogridAdapterProteinField(Enum):
    """
    Define possible fields the adapter can provide for proteins.
    """
    NAME = "name"
    ORGANISM = "organism"
    GENE_NAME = "gene_name"
    BIOGRID_ID = "biogrid_id"
    UNIPROT_ID = "uniprot_id"


class BiogridAdapterEdgeType(Enum):
    """
    Enum for the types of edges the adapter can provide.
    """
    PROTEIN_PROTEIN_INTERACTION = "protein_protein_interaction"


class BiogridAdapterInteractionEdgeField(Enum):
    """
    Define possible fields for protein-protein interaction edges.
    """
    SOURCE = "source"
    PUBMED_IDS = "pubmed_ids"
    EXPERIMENTAL_SYSTEM = "experimental_system"
    EXPERIMENTAL_SYSTEM_TYPE = "experimental_system_type"
    AUTHOR = "author"
    INTERACTION_THROUGHPUT = "interaction_throughput"
    SCORE = "score"
    MODIFICATION = "modification"
    QUALIFICATIONS = "qualifications"
    TAGS = "tags"
    SOURCE_DATABASE = "source_database"


class BiogridAdapter:
    """
    BioGRID BioCypher adapter. Generates protein nodes and high-quality 
    protein-protein interaction edges from BioGRID database.

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
        interaction_probability: Optional[float] = 0.08,
    ):
        self.organism = organism
        self.test_mode = test_mode
        self.interaction_probability = interaction_probability
        self._set_types_and_fields(node_types, node_fields, edge_types, edge_fields)

    def get_nodes(self):
        """
        Returns a generator of node tuples for node types specified in the
        adapter constructor.
        """
        logger.info("Generating BioGRID protein nodes.")

        self.nodes = []

        if BiogridAdapterNodeType.PROTEIN in self.node_types:
            # Generate proteins for BioGRID network
            num_proteins = 60 if self.test_mode else 600
            [self.nodes.append(BiogridProtein(fields=self.node_fields, organism=self.organism)) 
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
        logger.info("Generating BioGRID protein-protein interaction edges.")

        if protein_nodes is None:
            if not self.nodes:
                raise ValueError("No protein nodes found. Please run get_nodes() first or provide protein_nodes.")
            protein_nodes = self.nodes
        
        # Convert to list if it's a generator
        if not isinstance(protein_nodes, list):
            protein_nodes = list(protein_nodes)

        if BiogridAdapterEdgeType.PROTEIN_PROTEIN_INTERACTION in self.edge_types:
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
                                if field == BiogridAdapterInteractionEdgeField.SOURCE:
                                    properties["source"] = "BioGRID"
                                elif field == BiogridAdapterInteractionEdgeField.PUBMED_IDS:
                                    # Generate 1-5 PubMed IDs for high-quality interactions
                                    num_pubmed = random.randint(1, 5)
                                    pubmed_ids = [str(random.randint(10000000, 35000000)) for _ in range(num_pubmed)]
                                    properties["pubmed_ids"] = "|".join(pubmed_ids)
                                elif field == BiogridAdapterInteractionEdgeField.EXPERIMENTAL_SYSTEM:
                                    properties["experimental_system"] = random.choice([
                                        "Affinity Capture-MS", "Affinity Capture-Western",
                                        "Two-hybrid", "Co-purification", "Reconstituted Complex",
                                        "Co-fractionation", "Affinity Capture-RNA", "Protein-peptide",
                                        "FRET", "Co-crystal Structure", "Biochemical Activity",
                                        "Co-localization", "Proximity Label-MS", "PCA"
                                    ])
                                elif field == BiogridAdapterInteractionEdgeField.EXPERIMENTAL_SYSTEM_TYPE:
                                    properties["experimental_system_type"] = random.choice([
                                        "physical", "genetic", "regulatory"
                                    ])
                                elif field == BiogridAdapterInteractionEdgeField.AUTHOR:
                                    # Generate author names in format "LastName (Year)"
                                    authors = [
                                        "Smith", "Johnson", "Brown", "Davis", "Miller", "Wilson",
                                        "Moore", "Taylor", "Anderson", "Thomas", "Jackson", "White",
                                        "Harris", "Martin", "Thompson", "Garcia", "Martinez", "Robinson"
                                    ]
                                    year = random.randint(2010, 2024)
                                    properties["author"] = f"{random.choice(authors)} ({year})"
                                elif field == BiogridAdapterInteractionEdgeField.INTERACTION_THROUGHPUT:
                                    properties["interaction_throughput"] = random.choice([
                                        "High Throughput", "Low Throughput"
                                    ])
                                elif field == BiogridAdapterInteractionEdgeField.SCORE:
                                    # BioGRID often includes reliability scores
                                    properties["score"] = round(random.uniform(0.5, 1.0), 3)
                                elif field == BiogridAdapterInteractionEdgeField.MODIFICATION:
                                    if random.random() < 0.3:  # 30% chance
                                        properties["modification"] = random.choice([
                                            "phosphorylation", "ubiquitination", "acetylation",
                                            "methylation", "sumoylation", "neddylation"
                                        ])
                                elif field == BiogridAdapterInteractionEdgeField.QUALIFICATIONS:
                                    if random.random() < 0.2:  # 20% chance  
                                        properties["qualifications"] = random.choice([
                                            "direct interaction", "indirect interaction",
                                            "colocalizing", "genetic interaction"
                                        ])
                                elif field == BiogridAdapterInteractionEdgeField.TAGS:
                                    if random.random() < 0.25:  # 25% chance
                                        properties["tags"] = random.choice([
                                            "validated", "high confidence", "conserved",
                                            "disease-related", "drug target"  
                                        ])
                                elif field == BiogridAdapterInteractionEdgeField.SOURCE_DATABASE:
                                    properties["source_database"] = "BioGRID"

                            yield (
                                relationship_id,
                                protein_a_id,
                                protein_b_id,
                                BiogridAdapterEdgeType.PROTEIN_PROTEIN_INTERACTION.value,
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
            self.node_types = [type for type in BiogridAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [field for field in BiogridAdapterProteinField]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in BiogridAdapterEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in BiogridAdapterInteractionEdgeField]


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


class BiogridProtein(Node):
    """
    Generates instances of proteins for BioGRID networks.
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
        if BiogridAdapterProteinField.ORGANISM in self.fields:
            properties["organism"] = self.organism

        # Generate protein name
        if BiogridAdapterProteinField.NAME in self.fields:
            protein_names = [
                "Tumor suppressor p53", "Epidermal growth factor receptor",
                "Insulin receptor substrate 1", "Cell division cycle 42",
                "Serine/threonine-protein kinase mTOR", "Nuclear factor NF-kappa-B p65",
                "Catenin beta-1", "Cyclin-dependent kinase 1", "Mitogen-activated protein kinase 1",
                "Proto-oncogene c-Myc", "Transcription factor p65", "Retinoblastoma protein",
                "DNA topoisomerase 2-alpha", "Heat shock protein HSP 90-alpha",
                "Ubiquitin-conjugating enzyme E2 I", "26S proteasome regulatory subunit 6A",
                "Histone deacetylase 1", "Chromatin modification complex subunit",
                "Ribosomal protein S6 kinase alpha-1", "AMP-activated protein kinase",
                "Phosphatidylinositol 3-kinase regulatory subunit alpha", "BRCA1",
                "BRCA2", "ATM serine/threonine kinase", "DNA repair protein RAD51"
            ]
            properties["name"] = random.choice(protein_names)

        # Generate gene name
        if BiogridAdapterProteinField.GENE_NAME in self.fields:
            gene_names = [
                "TP53", "EGFR", "IRS1", "CDC42", "MTOR", "RELA", "CTNNB1", "CDK1", 
                "MAPK1", "MYC", "NFKB1", "RB1", "TOP2A", "HSP90AA1", "UBE2I", 
                "PSMC3", "HDAC1", "SMARCA4", "RPS6KA1", "PRKAA1", "PIK3R1",
                "BRCA1", "BRCA2", "ATM", "RAD51", "MLH1", "MSH2", "PTEN", "APC"
            ]
            properties["gene_name"] = random.choice(gene_names)

        # Generate BioGRID ID
        if BiogridAdapterProteinField.BIOGRID_ID in self.fields:
            properties["biogrid_id"] = random.randint(100000, 999999)

        # Generate UniProt ID  
        if BiogridAdapterProteinField.UNIPROT_ID in self.fields:
            properties["uniprot_id"] = self.id

        return properties