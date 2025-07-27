from __future__ import annotations

import random
import string
from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class GOAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    BIOLOGICAL_PROCESS = auto()
    MOLECULAR_FUNCTION = auto()
    CELLULAR_COMPONENT = auto()


class GOAdapterTermField(Enum):
    """
    Define possible fields the adapter can provide for GO terms.
    """
    ID = "id"
    NAME = "name"
    DEFINITION = "definition"
    NAMESPACE = "namespace"
    SYNONYMS = "synonyms"
    IS_OBSOLETE = "is_obsolete"


class GOAdapterEdgeType(Enum):
    """
    Enum for the types of edges the adapter can provide.
    """
    PROTEIN_HAS_GO_ANNOTATION = "protein_has_go_annotation"
    GO_TERM_IS_A_GO_TERM = "go_term_is_a_go_term"
    GO_TERM_PART_OF_GO_TERM = "go_term_part_of_go_term"


class GOAdapterProteinGOEdgeField(Enum):
    """
    Define possible fields for protein-GO annotation edges.
    """
    EVIDENCE_CODE = "evidence_code"
    EVIDENCE_SOURCE = "evidence_source"
    REFERENCE = "reference"
    QUALIFIER = "qualifier"
    ANNOTATION_DATE = "annotation_date"


class GOAdapterGOHierarchyEdgeField(Enum):
    """
    Define possible fields for GO term hierarchy edges.
    """
    RELATIONSHIP_TYPE = "relationship_type"
    SOURCE = "source"


class GOAdapter:
    """
    Gene Ontology BioCypher adapter. Generates GO term nodes and 
    protein-GO annotation edges for creating a knowledge graph.

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
        logger.info("Generating GO term nodes.")

        self.nodes = []

        # Generate biological process terms
        if GOAdapterNodeType.BIOLOGICAL_PROCESS in self.node_types:
            num_bp = 15 if self.test_mode else 100
            [self.nodes.append(BiologicalProcess(fields=self.node_fields)) 
             for _ in range(num_bp)]

        # Generate molecular function terms
        if GOAdapterNodeType.MOLECULAR_FUNCTION in self.node_types:
            num_mf = 10 if self.test_mode else 80
            [self.nodes.append(MolecularFunction(fields=self.node_fields)) 
             for _ in range(num_mf)]

        # Generate cellular component terms
        if GOAdapterNodeType.CELLULAR_COMPONENT in self.node_types:
            num_cc = 8 if self.test_mode else 60
            [self.nodes.append(CellularComponent(fields=self.node_fields)) 
             for _ in range(num_cc)]

        for node in self.nodes:
            yield (node.get_id(), node.get_label(), node.get_properties())

    def get_edges(self, target_proteins=None):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.

        Args:
            target_proteins: List of protein nodes to create annotations with.
        """
        logger.info("Generating GO annotation edges.")

        if not self.nodes:
            raise ValueError("No nodes found. Please run get_nodes() first.")

        # Create protein-GO annotation edges
        if (GOAdapterEdgeType.PROTEIN_HAS_GO_ANNOTATION in self.edge_types 
            and target_proteins):
            
            for go_term in self.nodes:
                # Each GO term is associated with 3-8 proteins
                num_annotations = random.randint(3, 8)
                selected_proteins = random.sample(
                    target_proteins, 
                    min(num_annotations, len(target_proteins))
                )
                
                for protein in selected_proteins:
                    relationship_id = "".join(
                        random.choice(string.ascii_letters + string.digits)
                        for _ in range(10)
                    )
                    
                    properties = {}
                    for field in self.edge_fields:
                        if field == GOAdapterProteinGOEdgeField.EVIDENCE_CODE:
                            properties["evidence_code"] = random.choice([
                                "IEA", "IDA", "IMP", "IGI", "IPI", "ISS", "TAS", "NAS", 
                                "IC", "ND", "IEP", "EXP", "HDA", "HMP", "HGI", "HEP"
                            ])
                        elif field == GOAdapterProteinGOEdgeField.EVIDENCE_SOURCE:
                            properties["evidence_source"] = random.choice([
                                "UniProtKB", "MGI", "SGD", "FlyBase", "WormBase", 
                                "TAIR", "RGD", "ZFIN", "PomBase"
                            ])
                        elif field == GOAdapterProteinGOEdgeField.REFERENCE:
                            properties["reference"] = f"PMID:{random.randint(10000000, 35000000)}"
                        elif field == GOAdapterProteinGOEdgeField.QUALIFIER:
                            if random.random() < 0.2:  # 20% have qualifiers
                                properties["qualifier"] = random.choice([
                                    "NOT", "contributes_to", "colocalizes_with"
                                ])
                        elif field == GOAdapterProteinGOEdgeField.ANNOTATION_DATE:
                            year = random.randint(2000, 2024)
                            month = random.randint(1, 12)
                            day = random.randint(1, 28)
                            properties["annotation_date"] = f"{year:04d}-{month:02d}-{day:02d}"

                    yield (
                        relationship_id,
                        protein.get_id() if hasattr(protein, 'get_id') else protein,
                        go_term.get_id(),
                        GOAdapterEdgeType.PROTEIN_HAS_GO_ANNOTATION.value,
                        properties,
                    )

        # Create GO term hierarchy edges (is_a relationships)
        if GOAdapterEdgeType.GO_TERM_IS_A_GO_TERM in self.edge_types:
            # Group nodes by namespace for hierarchy
            bp_nodes = [n for n in self.nodes if isinstance(n, BiologicalProcess)]
            mf_nodes = [n for n in self.nodes if isinstance(n, MolecularFunction)]
            cc_nodes = [n for n in self.nodes if isinstance(n, CellularComponent)]
            
            for node_group in [bp_nodes, mf_nodes, cc_nodes]:
                for i, child_term in enumerate(node_group):
                    if i > 0 and random.random() < 0.4:  # 40% chance of having a parent
                        parent_term = random.choice(node_group[:i])
                        
                        relationship_id = "".join(
                            random.choice(string.ascii_letters + string.digits)
                            for _ in range(10)
                        )
                        
                        properties = {}
                        for field in self.edge_fields:
                            if field == GOAdapterGOHierarchyEdgeField.RELATIONSHIP_TYPE:
                                properties["relationship_type"] = "is_a"
                            elif field == GOAdapterGOHierarchyEdgeField.SOURCE:
                                properties["source"] = "Gene Ontology Consortium"

                        yield (
                            relationship_id,
                            child_term.get_id(),
                            parent_term.get_id(),
                            GOAdapterEdgeType.GO_TERM_IS_A_GO_TERM.value,
                            properties,
                        )

        # Create GO term part_of relationships
        if GOAdapterEdgeType.GO_TERM_PART_OF_GO_TERM in self.edge_types:
            # Some terms have part_of relationships
            for node in self.nodes:
                if random.random() < 0.15:  # 15% chance of part_of relationship
                    # Find a potential parent from same namespace
                    same_namespace_nodes = [
                        n for n in self.nodes 
                        if type(n) == type(node) and n != node
                    ]
                    
                    if same_namespace_nodes:
                        parent_term = random.choice(same_namespace_nodes)
                        
                        relationship_id = "".join(
                            random.choice(string.ascii_letters + string.digits)
                            for _ in range(10)
                        )
                        
                        properties = {}
                        for field in self.edge_fields:
                            if field == GOAdapterGOHierarchyEdgeField.RELATIONSHIP_TYPE:
                                properties["relationship_type"] = "part_of"
                            elif field == GOAdapterGOHierarchyEdgeField.SOURCE:
                                properties["source"] = "Gene Ontology Consortium"

                        yield (
                            relationship_id,
                            node.get_id(),
                            parent_term.get_id(),
                            GOAdapterEdgeType.GO_TERM_PART_OF_GO_TERM.value,
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
            self.node_types = [type for type in GOAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [field for field in GOAdapterTermField]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in GOAdapterEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [
                field for field in chain(
                    GOAdapterProteinGOEdgeField,
                    GOAdapterGOHierarchyEdgeField,
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


class BiologicalProcess(Node):
    """
    Generates instances of GO Biological Process terms.
    """

    def __init__(self, fields: Optional[list] = None):
        self.fields = fields
        self.id = self._generate_id()
        self.label = "biological_process"
        self.properties = self._generate_properties()

    def _generate_id(self):
        """
        Generate a random GO ID for biological process.
        """
        return f"GO:{random.randint(1, 9999999):07d}"

    def _generate_properties(self):
        properties = {}

        # Generate term name
        if GOAdapterTermField.NAME in self.fields:
            bp_names = [
                "cell cycle", "apoptotic process", "DNA repair", "transcription",
                "translation", "protein folding", "signal transduction",
                "metabolic process", "transport", "immune response",
                "development", "cell differentiation", "homeostasis",
                "response to stimulus", "cellular respiration", "photosynthesis",
                "cell division", "protein modification", "RNA processing",
                "DNA replication"
            ]
            properties["name"] = random.choice(bp_names)

        # Generate definition
        if GOAdapterTermField.DEFINITION in self.fields:
            properties["definition"] = f"A biological process involving {properties.get('name', 'cellular activity')}"

        # Set namespace
        if GOAdapterTermField.NAMESPACE in self.fields:
            properties["namespace"] = "biological_process"

        # Generate synonyms
        if GOAdapterTermField.SYNONYMS in self.fields:
            if random.random() < 0.3:  # 30% have synonyms
                base_name = properties.get("name", "process")
                synonyms = [f"{base_name} pathway", f"regulation of {base_name}"]
                properties["synonyms"] = "|".join(synonyms)

        # Set obsolete status
        if GOAdapterTermField.IS_OBSOLETE in self.fields:
            properties["is_obsolete"] = random.random() < 0.05  # 5% are obsolete

        return properties


class MolecularFunction(Node):
    """
    Generates instances of GO Molecular Function terms.
    """

    def __init__(self, fields: Optional[list] = None):
        self.fields = fields
        self.id = self._generate_id()
        self.label = "molecular_function"
        self.properties = self._generate_properties()

    def _generate_id(self):
        """
        Generate a random GO ID for molecular function.
        """
        return f"GO:{random.randint(1, 9999999):07d}"

    def _generate_properties(self):
        properties = {}

        # Generate term name
        if GOAdapterTermField.NAME in self.fields:
            mf_names = [
                "protein kinase activity", "DNA binding", "RNA binding",
                "enzyme activity", "transcription factor activity",
                "transporter activity", "receptor activity", "catalytic activity",
                "binding", "structural molecule activity", "antioxidant activity",
                "molecular transducer activity", "protein phosphatase activity",
                "hydrolase activity", "transferase activity", "oxidoreductase activity",
                "ligase activity", "lyase activity", "isomerase activity",
                "motor activity"
            ]
            properties["name"] = random.choice(mf_names)

        # Generate definition
        if GOAdapterTermField.DEFINITION in self.fields:
            properties["definition"] = f"Molecular function related to {properties.get('name', 'protein activity')}"

        # Set namespace
        if GOAdapterTermField.NAMESPACE in self.fields:
            properties["namespace"] = "molecular_function"

        # Generate synonyms
        if GOAdapterTermField.SYNONYMS in self.fields:
            if random.random() < 0.3:  # 30% have synonyms
                base_name = properties.get("name", "activity")
                synonyms = [base_name.replace("activity", "function")]
                properties["synonyms"] = "|".join(synonyms)

        # Set obsolete status
        if GOAdapterTermField.IS_OBSOLETE in self.fields:
            properties["is_obsolete"] = random.random() < 0.05  # 5% are obsolete

        return properties


class CellularComponent(Node):
    """
    Generates instances of GO Cellular Component terms.
    """

    def __init__(self, fields: Optional[list] = None):
        self.fields = fields
        self.id = self._generate_id()
        self.label = "cellular_component"
        self.properties = self._generate_properties()

    def _generate_id(self):
        """
        Generate a random GO ID for cellular component.
        """
        return f"GO:{random.randint(1, 9999999):07d}"

    def _generate_properties(self):
        properties = {}

        # Generate term name
        if GOAdapterTermField.NAME in self.fields:
            cc_names = [
                "nucleus", "cytoplasm", "mitochondrion", "endoplasmic reticulum",
                "Golgi apparatus", "ribosome", "lysosome", "peroxisome",
                "cell membrane", "nuclear membrane", "cytoskeleton",
                "chromosome", "nucleolus", "vacuole", "chloroplast",
                "cell wall", "centriole", "flagellum", "cilium", "vesicle"
            ]
            properties["name"] = random.choice(cc_names)

        # Generate definition
        if GOAdapterTermField.DEFINITION in self.fields:
            properties["definition"] = f"A cellular component referring to {properties.get('name', 'cellular structure')}"

        # Set namespace
        if GOAdapterTermField.NAMESPACE in self.fields:
            properties["namespace"] = "cellular_component"

        # Generate synonyms
        if GOAdapterTermField.SYNONYMS in self.fields:
            if random.random() < 0.2:  # 20% have synonyms
                base_name = properties.get("name", "component")
                synonyms = [f"{base_name} complex"]
                properties["synonyms"] = "|".join(synonyms)

        # Set obsolete status
        if GOAdapterTermField.IS_OBSOLETE in self.fields:
            properties["is_obsolete"] = random.random() < 0.05  # 5% are obsolete

        return properties