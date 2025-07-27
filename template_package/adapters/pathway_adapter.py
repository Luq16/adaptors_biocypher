from __future__ import annotations

import random
import string
from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class PathwayAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    PATHWAY = auto()


class PathwayAdapterPathwayField(Enum):
    """
    Define possible fields the adapter can provide for pathways.
    """
    ID = "id"
    NAME = "name"
    DESCRIPTION = "description"
    ORGANISM = "organism"
    PATHWAY_TYPE = "pathway_type"
    DATABASE_SOURCE = "database_source"
    EXTERNAL_IDS = "external_ids"
    GENE_COUNT = "gene_count"
    COMPOUND_COUNT = "compound_count"


class PathwayAdapterEdgeType(Enum):
    """
    Enum for the types of edges the adapter can provide.
    """
    PROTEIN_PARTICIPATES_IN_PATHWAY = "protein_participates_in_pathway"
    COMPOUND_PARTICIPATES_IN_PATHWAY = "compound_participates_in_pathway"
    PATHWAY_IS_PART_OF_PATHWAY = "pathway_is_part_of_pathway"


class PathwayAdapterProteinPathwayEdgeField(Enum):
    """
    Define possible fields for protein-pathway edges.
    """
    ROLE = "role"
    EVIDENCE_CODE = "evidence_code"
    SOURCE = "source"
    REFERENCE = "reference"


class PathwayAdapterCompoundPathwayEdgeField(Enum):
    """
    Define possible fields for compound-pathway edges.
    """
    ROLE = "role"
    REACTION_TYPE = "reaction_type"
    SOURCE = "source"


class PathwayAdapterPathwayHierarchyEdgeField(Enum):
    """
    Define possible fields for pathway hierarchy edges.
    """
    RELATIONSHIP_TYPE = "relationship_type"
    SOURCE = "source"


class PathwayAdapter:
    """
    Pathway BioCypher adapter. Generates pathway nodes and protein/compound-pathway 
    participation edges for creating a knowledge graph.

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
        logger.info("Generating pathway nodes.")

        self.nodes = []

        if PathwayAdapterNodeType.PATHWAY in self.node_types:
            # In a real implementation, this would fetch from KEGG, Reactome, WikiPathways, etc.
            # For now, generating example pathways
            num_pathways = 20 if self.test_mode else 150
            [self.nodes.append(Pathway(fields=self.node_fields, organism=self.organism)) 
             for _ in range(num_pathways)]

        for node in self.nodes:
            yield (node.get_id(), node.get_label(), node.get_properties())

    def get_edges(self, target_proteins=None, target_compounds=None):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.

        Args:
            target_proteins: List of protein nodes to create associations with.
            target_compounds: List of compound nodes to create associations with.
        """
        logger.info("Generating pathway participation edges.")

        if not self.nodes:
            raise ValueError("No nodes found. Please run get_nodes() first.")

        # Create protein-pathway participation edges
        if (PathwayAdapterEdgeType.PROTEIN_PARTICIPATES_IN_PATHWAY in self.edge_types 
            and target_proteins):
            
            for pathway in self.nodes:
                # Each pathway involves 5-20 proteins
                num_proteins = random.randint(5, 20)
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
                        if field == PathwayAdapterProteinPathwayEdgeField.ROLE:
                            properties["role"] = random.choice([
                                "enzyme", "substrate", "product", "regulator", 
                                "inhibitor", "activator", "cofactor", "transporter"
                            ])
                        elif field == PathwayAdapterProteinPathwayEdgeField.EVIDENCE_CODE:
                            properties["evidence_code"] = random.choice([
                                "IEA", "IDA", "IMP", "TAS", "NAS", "IC", "ISS", "IPI"
                            ])
                        elif field == PathwayAdapterProteinPathwayEdgeField.SOURCE:
                            properties["source"] = random.choice([
                                "KEGG", "Reactome", "WikiPathways", "BioCyc", 
                                "PANTHER", "GO", "Pathway Commons"
                            ])
                        elif field == PathwayAdapterProteinPathwayEdgeField.REFERENCE:
                            properties["reference"] = f"PMID:{random.randint(10000000, 35000000)}"

                    yield (
                        relationship_id,
                        protein.get_id() if hasattr(protein, 'get_id') else protein,
                        pathway.get_id(),
                        PathwayAdapterEdgeType.PROTEIN_PARTICIPATES_IN_PATHWAY.value,
                        properties,
                    )

        # Create compound-pathway participation edges
        if (PathwayAdapterEdgeType.COMPOUND_PARTICIPATES_IN_PATHWAY in self.edge_types 
            and target_compounds):
            
            for pathway in self.nodes:
                # Each pathway involves 3-10 compounds
                num_compounds = random.randint(3, 10)
                selected_compounds = random.sample(
                    target_compounds, 
                    min(num_compounds, len(target_compounds))
                )
                
                for compound in selected_compounds:
                    relationship_id = "".join(
                        random.choice(string.ascii_letters + string.digits)
                        for _ in range(10)
                    )
                    
                    properties = {}
                    for field in self.edge_fields:
                        if field == PathwayAdapterCompoundPathwayEdgeField.ROLE:
                            properties["role"] = random.choice([
                                "substrate", "product", "cofactor", "inhibitor", 
                                "activator", "intermediate", "precursor", "metabolite"
                            ])
                        elif field == PathwayAdapterCompoundPathwayEdgeField.REACTION_TYPE:
                            properties["reaction_type"] = random.choice([
                                "oxidation", "reduction", "hydrolysis", "condensation",
                                "phosphorylation", "methylation", "acetylation", "glycosylation"
                            ])
                        elif field == PathwayAdapterCompoundPathwayEdgeField.SOURCE:
                            properties["source"] = random.choice([
                                "KEGG", "Reactome", "WikiPathways", "BioCyc", "ChEBI"
                            ])

                    yield (
                        relationship_id,
                        compound.get_id() if hasattr(compound, 'get_id') else compound,
                        pathway.get_id(),
                        PathwayAdapterEdgeType.COMPOUND_PARTICIPATES_IN_PATHWAY.value,
                        properties,
                    )

        # Create pathway hierarchy edges
        if PathwayAdapterEdgeType.PATHWAY_IS_PART_OF_PATHWAY in self.edge_types:
            # Create some parent-child relationships between pathways
            for i, child_pathway in enumerate(self.nodes):
                if i > 0 and random.random() < 0.3:  # 30% chance of having a parent
                    parent_pathway = random.choice(self.nodes[:i])  # Choose from earlier pathways
                    
                    relationship_id = "".join(
                        random.choice(string.ascii_letters + string.digits)
                        for _ in range(10)
                    )
                    
                    properties = {}
                    for field in self.edge_fields:
                        if field == PathwayAdapterPathwayHierarchyEdgeField.RELATIONSHIP_TYPE:
                            properties["relationship_type"] = random.choice([
                                "is_part_of", "is_subpathway_of", "regulates", 
                                "positively_regulates", "negatively_regulates"
                            ])
                        elif field == PathwayAdapterPathwayHierarchyEdgeField.SOURCE:
                            properties["source"] = random.choice([
                                "KEGG", "Reactome", "GO", "WikiPathways"
                            ])

                    yield (
                        relationship_id,
                        child_pathway.get_id(),
                        parent_pathway.get_id(),
                        PathwayAdapterEdgeType.PATHWAY_IS_PART_OF_PATHWAY.value,
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
            self.node_types = [type for type in PathwayAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [field for field in PathwayAdapterPathwayField]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in PathwayAdapterEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [
                field for field in chain(
                    PathwayAdapterProteinPathwayEdgeField,
                    PathwayAdapterCompoundPathwayEdgeField,
                    PathwayAdapterPathwayHierarchyEdgeField,
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


class Pathway(Node):
    """
    Generates instances of biological pathways.
    """

    def __init__(self, fields: Optional[list] = None, organism: str = "9606"):
        self.fields = fields
        self.organism = organism
        self.id = self._generate_id()
        self.label = "pathway"
        self.properties = self._generate_properties()

    def _generate_id(self):
        """
        Generate a random pathway ID (KEGG-style).
        """
        organism_code = "hsa" if self.organism == "9606" else "org"
        pathway_number = f"{random.randint(1, 99999):05d}"
        return f"{organism_code}{pathway_number}"

    def _generate_properties(self):
        properties = {}

        # Generate pathway name
        if PathwayAdapterPathwayField.NAME in self.fields:
            pathway_names = [
                "Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)",
                "Pentose phosphate pathway", "Fatty acid biosynthesis",
                "Fatty acid degradation", "Steroid biosynthesis",
                "Purine metabolism", "Pyrimidine metabolism",
                "Amino acid biosynthesis", "Amino acid degradation",
                "Cell cycle", "Apoptosis", "DNA replication", "DNA repair",
                "mRNA surveillance pathway", "Protein processing in endoplasmic reticulum",
                "MAPK signaling pathway", "PI3K-Akt signaling pathway",
                "p53 signaling pathway", "Wnt signaling pathway",
                "Notch signaling pathway", "Hedgehog signaling pathway",
                "TGF-beta signaling pathway", "JAK-STAT signaling pathway",
                "NF-kappa B signaling pathway", "mTOR signaling pathway",
                "Insulin signaling pathway", "VEGF signaling pathway",
                "ErbB signaling pathway", "Calcium signaling pathway"
            ]
            properties["name"] = random.choice(pathway_names)

        # Generate description
        if PathwayAdapterPathwayField.DESCRIPTION in self.fields:
            pathway_name = properties.get("name", "pathway")
            descriptions = [
                f"A metabolic pathway involved in {pathway_name.lower()}",
                f"A signaling pathway that regulates {pathway_name.lower()}",
                f"A cellular process pathway controlling {pathway_name.lower()}",
                f"An important biological pathway for {pathway_name.lower()}"
            ]
            properties["description"] = random.choice(descriptions)

        # Set organism
        if PathwayAdapterPathwayField.ORGANISM in self.fields:
            properties["organism"] = self.organism

        # Generate pathway type
        if PathwayAdapterPathwayField.PATHWAY_TYPE in self.fields:
            pathway_types = [
                "metabolic pathway", "signaling pathway", "regulatory pathway",
                "transport pathway", "cellular process", "genetic information processing",
                "environmental information processing", "organismal systems"
            ]
            properties["pathway_type"] = random.choice(pathway_types)

        # Generate database source
        if PathwayAdapterPathwayField.DATABASE_SOURCE in self.fields:
            properties["database_source"] = random.choice([
                "KEGG", "Reactome", "WikiPathways", "BioCyc", "PANTHER", "GO"
            ])

        # Generate external IDs
        if PathwayAdapterPathwayField.EXTERNAL_IDS in self.fields:
            external_ids = []
            if random.random() < 0.8:  # 80% have Reactome ID
                external_ids.append(f"R-HSA-{random.randint(100000, 999999)}")
            if random.random() < 0.6:  # 60% have GO ID
                external_ids.append(f"GO:{random.randint(1, 9999999):07d}")
            if random.random() < 0.4:  # 40% have WikiPathways ID
                external_ids.append(f"WP{random.randint(1, 9999)}")
            
            if external_ids:
                properties["external_ids"] = "|".join(external_ids)

        # Generate gene count
        if PathwayAdapterPathwayField.GENE_COUNT in self.fields:
            properties["gene_count"] = random.randint(5, 200)

        # Generate compound count
        if PathwayAdapterPathwayField.COMPOUND_COUNT in self.fields:
            properties["compound_count"] = random.randint(2, 50)

        return properties