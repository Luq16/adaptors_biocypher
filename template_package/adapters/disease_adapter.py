from __future__ import annotations

import random
import string
from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class DiseaseAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    DISEASE = auto()


class DiseaseAdapterDiseaseField(Enum):
    """
    Define possible fields the adapter can provide for diseases.
    """
    ID = "id"
    NAME = "name"
    DESCRIPTION = "description"
    SYNONYMS = "synonyms"
    DOID = "doid"
    MESH_ID = "mesh_id"
    OMIM_ID = "omim_id"
    ICD10_CODE = "icd10_code"
    DISEASE_TYPE = "disease_type"
    AFFECTED_SYSTEM = "affected_system"


class DiseaseAdapterEdgeType(Enum):
    """
    Enum for the types of edges the adapter can provide.
    """
    GENE_ASSOCIATED_WITH_DISEASE = "gene_associated_with_disease"
    DISEASE_IS_SUBTYPE_OF_DISEASE = "disease_is_subtype_of_disease"


class DiseaseAdapterGeneDiseaseEdgeField(Enum):
    """
    Define possible fields for gene-disease association edges.
    """
    ASSOCIATION_TYPE = "association_type"
    CONFIDENCE_SCORE = "confidence_score"
    EVIDENCE_TYPE = "evidence_type"
    SOURCE = "source"
    PUBMED_IDS = "pubmed_ids"


class DiseaseAdapterDiseaseHierarchyEdgeField(Enum):
    """
    Define possible fields for disease hierarchy edges.
    """
    RELATIONSHIP_TYPE = "relationship_type"
    SOURCE = "source"


class DiseaseAdapter:
    """
    Disease BioCypher adapter. Generates disease nodes and disease-gene 
    association edges for creating a knowledge graph.

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
        logger.info("Generating disease nodes.")

        self.nodes = []

        if DiseaseAdapterNodeType.DISEASE in self.node_types:
            # In a real implementation, this would fetch from MONDO, DO, OMIM, etc.
            # For now, generating example diseases
            num_diseases = 25 if self.test_mode else 200
            [self.nodes.append(Disease(fields=self.node_fields)) 
             for _ in range(num_diseases)]

        for node in self.nodes:
            yield (node.get_id(), node.get_label(), node.get_properties())

    def get_edges(self, target_genes=None):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.

        Args:
            target_genes: List of gene/protein nodes to create associations with.
        """
        logger.info("Generating disease association edges.")

        if not self.nodes:
            raise ValueError("No nodes found. Please run get_nodes() first.")

        # Create gene-disease association edges
        if (DiseaseAdapterEdgeType.GENE_ASSOCIATED_WITH_DISEASE in self.edge_types 
            and target_genes):
            
            for disease in self.nodes:
                # Each disease is associated with 2-5 genes
                num_associations = random.randint(2, 5)
                selected_genes = random.sample(
                    target_genes, 
                    min(num_associations, len(target_genes))
                )
                
                for gene in selected_genes:
                    relationship_id = "".join(
                        random.choice(string.ascii_letters + string.digits)
                        for _ in range(10)
                    )
                    
                    properties = {}
                    for field in self.edge_fields:
                        if field == DiseaseAdapterGeneDiseaseEdgeField.ASSOCIATION_TYPE:
                            properties["association_type"] = random.choice([
                                "causative", "risk_factor", "protective", "modifier", 
                                "biomarker", "therapeutic_target"
                            ])
                        elif field == DiseaseAdapterGeneDiseaseEdgeField.CONFIDENCE_SCORE:
                            properties["confidence_score"] = round(random.uniform(0.1, 1.0), 3)
                        elif field == DiseaseAdapterGeneDiseaseEdgeField.EVIDENCE_TYPE:
                            properties["evidence_type"] = random.choice([
                                "genetic_association", "functional_study", "expression_study",
                                "literature_mining", "clinical_study", "animal_model"
                            ])
                        elif field == DiseaseAdapterGeneDiseaseEdgeField.SOURCE:
                            properties["source"] = random.choice([
                                "DisGeNET", "OMIM", "ClinVar", "GWAS Catalog", "UniProt", "PubMed"
                            ])
                        elif field == DiseaseAdapterGeneDiseaseEdgeField.PUBMED_IDS:
                            # Generate 1-3 fake PubMed IDs
                            num_pubmed = random.randint(1, 3)
                            pubmed_ids = [str(random.randint(10000000, 35000000)) for _ in range(num_pubmed)]
                            properties["pubmed_ids"] = "|".join(pubmed_ids)

                    yield (
                        relationship_id,
                        gene.get_id() if hasattr(gene, 'get_id') else gene,
                        disease.get_id(),
                        DiseaseAdapterEdgeType.GENE_ASSOCIATED_WITH_DISEASE.value,
                        properties,
                    )

        # Create disease hierarchy edges
        if DiseaseAdapterEdgeType.DISEASE_IS_SUBTYPE_OF_DISEASE in self.edge_types:
            # Create some parent-child relationships between diseases
            for i, disease in enumerate(self.nodes):
                if i > 0 and random.random() < 0.3:  # 30% chance of having a parent
                    parent_disease = random.choice(self.nodes[:i])  # Choose from earlier diseases
                    
                    relationship_id = "".join(
                        random.choice(string.ascii_letters + string.digits)
                        for _ in range(10)
                    )
                    
                    properties = {}
                    for field in self.edge_fields:
                        if field == DiseaseAdapterDiseaseHierarchyEdgeField.RELATIONSHIP_TYPE:
                            properties["relationship_type"] = random.choice([
                                "is_a", "part_of", "subtype_of", "manifestation_of"
                            ])
                        elif field == DiseaseAdapterDiseaseHierarchyEdgeField.SOURCE:
                            properties["source"] = random.choice([
                                "MONDO", "Disease Ontology", "OMIM", "MeSH", "ICD-10"
                            ])

                    yield (
                        relationship_id,
                        disease.get_id(),
                        parent_disease.get_id(),
                        DiseaseAdapterEdgeType.DISEASE_IS_SUBTYPE_OF_DISEASE.value,
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
            self.node_types = [type for type in DiseaseAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [field for field in DiseaseAdapterDiseaseField]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in DiseaseAdapterEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [
                field for field in chain(
                    DiseaseAdapterGeneDiseaseEdgeField,
                    DiseaseAdapterDiseaseHierarchyEdgeField,
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


class Disease(Node):
    """
    Generates instances of diseases.
    """

    def __init__(self, fields: Optional[list] = None):
        self.fields = fields
        self.id = self._generate_id()
        self.label = "disease"
        self.properties = self._generate_properties()

    def _generate_id(self):
        """
        Generate a random disease ID (DOID-style).
        """
        return f"DOID:{random.randint(1, 99999):05d}"

    def _generate_properties(self):
        properties = {}

        # Generate disease name
        if DiseaseAdapterDiseaseField.NAME in self.fields:
            disease_names = [
                "Alzheimer's disease", "Parkinson's disease", "Type 2 diabetes", 
                "Hypertension", "Coronary artery disease", "Breast cancer", 
                "Lung cancer", "Asthma", "Rheumatoid arthritis", "Depression",
                "Schizophrenia", "Epilepsy", "Multiple sclerosis", "Crohn's disease",
                "Ulcerative colitis", "Osteoporosis", "Glaucoma", "Migraine",
                "Bipolar disorder", "Autism spectrum disorder", "Huntington's disease",
                "Cystic fibrosis", "Sickle cell anemia", "Hemophilia", "Muscular dystrophy"
            ]
            properties["name"] = random.choice(disease_names)

        # Generate description
        if DiseaseAdapterDiseaseField.DESCRIPTION in self.fields:
            descriptions = [
                "Progressive neurodegenerative disorder", "Chronic metabolic disorder",
                "Inflammatory autoimmune condition", "Malignant neoplasm",
                "Cardiovascular disorder", "Respiratory condition",
                "Psychiatric disorder", "Genetic disorder", "Infectious disease",
                "Endocrine disorder"
            ]
            properties["description"] = random.choice(descriptions)

        # Generate synonyms
        if DiseaseAdapterDiseaseField.SYNONYMS in self.fields:
            # Generate 1-3 synonyms
            base_name = properties.get("name", "Disease")
            synonyms = []
            if random.random() < 0.5:
                synonyms.append(f"{base_name} syndrome")
            if random.random() < 0.3:
                synonyms.append(f"Familial {base_name.lower()}")
            if random.random() < 0.2:
                synonyms.append(f"{base_name} type 1")
            
            if synonyms:
                properties["synonyms"] = "|".join(synonyms)

        # Generate DOID
        if DiseaseAdapterDiseaseField.DOID in self.fields:
            properties["doid"] = self.id

        # Generate MeSH ID
        if DiseaseAdapterDiseaseField.MESH_ID in self.fields:
            properties["mesh_id"] = f"D{random.randint(1000, 999999):06d}"

        # Generate OMIM ID
        if DiseaseAdapterDiseaseField.OMIM_ID in self.fields:
            if random.random() < 0.4:  # Not all diseases have OMIM IDs
                properties["omim_id"] = f"{random.randint(100000, 699999)}"

        # Generate ICD-10 code
        if DiseaseAdapterDiseaseField.ICD10_CODE in self.fields:
            letters = random.choices("ABCDEFGHIJKLMNOPQRSTUVWXYZ", k=1)[0]
            numbers = f"{random.randint(0, 99):02d}"
            decimal = f".{random.randint(0, 9)}" if random.random() < 0.5 else ""
            properties["icd10_code"] = f"{letters}{numbers}{decimal}"

        # Generate disease type
        if DiseaseAdapterDiseaseField.DISEASE_TYPE in self.fields:
            disease_types = [
                "genetic disorder", "autoimmune disease", "infectious disease",
                "neoplasm", "metabolic disorder", "neurological disorder",
                "cardiovascular disease", "respiratory disease", "psychiatric disorder",
                "rare disease"
            ]
            properties["disease_type"] = random.choice(disease_types)

        # Generate affected system
        if DiseaseAdapterDiseaseField.AFFECTED_SYSTEM in self.fields:
            systems = [
                "nervous system", "cardiovascular system", "respiratory system",
                "digestive system", "immune system", "endocrine system",
                "musculoskeletal system", "integumentary system", "urogenital system",
                "hematopoietic system"
            ]
            properties["affected_system"] = random.choice(systems)

        return properties