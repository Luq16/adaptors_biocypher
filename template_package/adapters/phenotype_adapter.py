from __future__ import annotations

import random
import string
from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class PhenotypeAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    PHENOTYPE = auto()


class PhenotypeAdapterPhenotypeField(Enum):
    """
    Define possible fields the adapter can provide for phenotypes.
    """
    ID = "id"
    NAME = "name"
    DESCRIPTION = "description"
    SYNONYMS = "synonyms"
    HPO_ID = "hpo_id"
    MAMMALIAN_PHENOTYPE_ID = "mp_id"
    ORGANISM = "organism"
    PHENOTYPE_TYPE = "phenotype_type"
    SEVERITY = "severity"
    ONSET = "onset"


class PhenotypeAdapterEdgeType(Enum):
    """
    Enum for the types of edges the adapter can provide.
    """
    GENE_HAS_PHENOTYPE = "gene_has_phenotype"
    DISEASE_HAS_PHENOTYPE = "disease_has_phenotype"
    PHENOTYPE_IS_SUBTYPE_OF_PHENOTYPE = "phenotype_is_subtype_of_phenotype"


class PhenotypeAdapterGenePhenotypeEdgeField(Enum):
    """
    Define possible fields for gene-phenotype association edges.
    """
    ASSOCIATION_TYPE = "association_type"
    EVIDENCE_CODE = "evidence_code"
    PENETRANCE = "penetrance"
    EXPRESSIVITY = "expressivity"
    SOURCE = "source"
    REFERENCE = "reference"


class PhenotypeAdapterDiseasePhenotypeEdgeField(Enum):
    """
    Define possible fields for disease-phenotype association edges.
    """
    FREQUENCY = "frequency"
    ONSET = "onset"
    SEVERITY = "severity"
    EVIDENCE_CODE = "evidence_code"
    SOURCE = "source"


class PhenotypeAdapterPhenotypeHierarchyEdgeField(Enum):
    """
    Define possible fields for phenotype hierarchy edges.
    """
    RELATIONSHIP_TYPE = "relationship_type"
    SOURCE = "source"


class PhenotypeAdapter:
    """
    Phenotype BioCypher adapter. Generates phenotype nodes and gene/disease-phenotype 
    association edges for creating a knowledge graph.

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
        logger.info("Generating phenotype nodes.")

        self.nodes = []

        if PhenotypeAdapterNodeType.PHENOTYPE in self.node_types:
            # In a real implementation, this would fetch from HPO, MP, etc.
            # For now, generating example phenotypes
            num_phenotypes = 30 if self.test_mode else 300
            [self.nodes.append(Phenotype(fields=self.node_fields, organism=self.organism)) 
             for _ in range(num_phenotypes)]

        for node in self.nodes:
            yield (node.get_id(), node.get_label(), node.get_properties())

    def get_edges(self, target_genes=None, target_diseases=None):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.

        Args:
            target_genes: List of gene/protein nodes to create associations with.
            target_diseases: List of disease nodes to create associations with.
        """
        logger.info("Generating phenotype association edges.")

        if not self.nodes:
            raise ValueError("No nodes found. Please run get_nodes() first.")

        # Create gene-phenotype association edges
        if (PhenotypeAdapterEdgeType.GENE_HAS_PHENOTYPE in self.edge_types 
            and target_genes):
            
            for phenotype in self.nodes:
                # Each phenotype is associated with 1-5 genes
                num_genes = random.randint(1, 5)
                selected_genes = random.sample(
                    target_genes, 
                    min(num_genes, len(target_genes))
                )
                
                for gene in selected_genes:
                    relationship_id = "".join(
                        random.choice(string.ascii_letters + string.digits)
                        for _ in range(10)
                    )
                    
                    properties = {}
                    for field in self.edge_fields:
                        if field == PhenotypeAdapterGenePhenotypeEdgeField.ASSOCIATION_TYPE:
                            properties["association_type"] = random.choice([
                                "causative", "contributory", "modifying", "protective",
                                "risk_factor", "biomarker"
                            ])
                        elif field == PhenotypeAdapterGenePhenotypeEdgeField.EVIDENCE_CODE:
                            properties["evidence_code"] = random.choice([
                                "IMP", "IGI", "IPI", "ISS", "IEA", "TAS", "NAS", "IC"
                            ])
                        elif field == PhenotypeAdapterGenePhenotypeEdgeField.PENETRANCE:
                            properties["penetrance"] = random.choice([
                                "complete", "incomplete", "age-dependent", "sex-dependent"
                            ])
                        elif field == PhenotypeAdapterGenePhenotypeEdgeField.EXPRESSIVITY:
                            properties["expressivity"] = random.choice([
                                "variable", "consistent", "mild", "moderate", "severe"
                            ])
                        elif field == PhenotypeAdapterGenePhenotypeEdgeField.SOURCE:
                            properties["source"] = random.choice([
                                "HPO", "MGI", "OMIM", "Orphanet", "ClinVar", "GWAS Catalog"
                            ])
                        elif field == PhenotypeAdapterGenePhenotypeEdgeField.REFERENCE:
                            properties["reference"] = f"PMID:{random.randint(10000000, 35000000)}"

                    yield (
                        relationship_id,
                        gene.get_id() if hasattr(gene, 'get_id') else gene,
                        phenotype.get_id(),
                        PhenotypeAdapterEdgeType.GENE_HAS_PHENOTYPE.value,
                        properties,
                    )

        # Create disease-phenotype association edges
        if (PhenotypeAdapterEdgeType.DISEASE_HAS_PHENOTYPE in self.edge_types 
            and target_diseases):
            
            for phenotype in self.nodes:
                # Each phenotype is associated with 1-3 diseases
                num_diseases = random.randint(1, 3)
                selected_diseases = random.sample(
                    target_diseases, 
                    min(num_diseases, len(target_diseases))
                )
                
                for disease in selected_diseases:
                    relationship_id = "".join(
                        random.choice(string.ascii_letters + string.digits)
                        for _ in range(10)
                    )
                    
                    properties = {}
                    for field in self.edge_fields:
                        if field == PhenotypeAdapterDiseasePhenotypeEdgeField.FREQUENCY:
                            properties["frequency"] = random.choice([
                                "very rare", "rare", "occasional", "frequent", 
                                "very frequent", "obligate"
                            ])
                        elif field == PhenotypeAdapterDiseasePhenotypeEdgeField.ONSET:
                            properties["onset"] = random.choice([
                                "congenital", "neonatal", "infantile", "childhood",
                                "juvenile", "adult", "late_adult", "variable"
                            ])
                        elif field == PhenotypeAdapterDiseasePhenotypeEdgeField.SEVERITY:
                            properties["severity"] = random.choice([
                                "mild", "moderate", "severe", "profound", "variable"
                            ])
                        elif field == PhenotypeAdapterDiseasePhenotypeEdgeField.EVIDENCE_CODE:
                            properties["evidence_code"] = random.choice([
                                "IEA", "TAS", "NAS", "IDA", "IC", "IMP"
                            ])
                        elif field == PhenotypeAdapterDiseasePhenotypeEdgeField.SOURCE:
                            properties["source"] = random.choice([
                                "HPO", "OMIM", "Orphanet", "DECIPHER", "ClinVar"
                            ])

                    yield (
                        relationship_id,
                        disease.get_id() if hasattr(disease, 'get_id') else disease,
                        phenotype.get_id(),
                        PhenotypeAdapterEdgeType.DISEASE_HAS_PHENOTYPE.value,
                        properties,
                    )

        # Create phenotype hierarchy edges
        if PhenotypeAdapterEdgeType.PHENOTYPE_IS_SUBTYPE_OF_PHENOTYPE in self.edge_types:
            # Create some parent-child relationships between phenotypes
            for i, child_phenotype in enumerate(self.nodes):
                if i > 0 and random.random() < 0.4:  # 40% chance of having a parent
                    parent_phenotype = random.choice(self.nodes[:i])  # Choose from earlier phenotypes
                    
                    relationship_id = "".join(
                        random.choice(string.ascii_letters + string.digits)
                        for _ in range(10)
                    )
                    
                    properties = {}
                    for field in self.edge_fields:
                        if field == PhenotypeAdapterPhenotypeHierarchyEdgeField.RELATIONSHIP_TYPE:
                            properties["relationship_type"] = random.choice([
                                "is_a", "part_of", "subtype_of"
                            ])
                        elif field == PhenotypeAdapterPhenotypeHierarchyEdgeField.SOURCE:
                            properties["source"] = random.choice([
                                "HPO", "MP", "Mammalian Phenotype Ontology"
                            ])

                    yield (
                        relationship_id,
                        child_phenotype.get_id(),
                        parent_phenotype.get_id(),
                        PhenotypeAdapterEdgeType.PHENOTYPE_IS_SUBTYPE_OF_PHENOTYPE.value,
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
            self.node_types = [type for type in PhenotypeAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [field for field in PhenotypeAdapterPhenotypeField]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in PhenotypeAdapterEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [
                field for field in chain(
                    PhenotypeAdapterGenePhenotypeEdgeField,
                    PhenotypeAdapterDiseasePhenotypeEdgeField,
                    PhenotypeAdapterPhenotypeHierarchyEdgeField,
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


class Phenotype(Node):
    """
    Generates instances of phenotypes.
    """

    def __init__(self, fields: Optional[list] = None, organism: str = "9606"):
        self.fields = fields
        self.organism = organism
        self.id = self._generate_id()
        self.label = "phenotype"
        self.properties = self._generate_properties()

    def _generate_id(self):
        """
        Generate a random phenotype ID (HPO-style for human, MP-style for mouse).
        """
        if self.organism == "9606":  # Human
            return f"HP:{random.randint(1, 9999999):07d}"
        elif self.organism == "10090":  # Mouse
            return f"MP:{random.randint(1, 9999999):07d}"
        else:
            return f"PHENO:{random.randint(1, 9999999):07d}"

    def _generate_properties(self):
        properties = {}

        # Generate phenotype name based on organism
        if PhenotypeAdapterPhenotypeField.NAME in self.fields:
            if self.organism == "9606":  # Human phenotypes
                phenotype_names = [
                    "Intellectual disability", "Seizure", "Growth delay", "Microcephaly",
                    "Macrocephaly", "Autism", "Ataxia", "Spasticity", "Muscle weakness",
                    "Hearing loss", "Visual impairment", "Cardiac arrhythmia",
                    "Cardiomyopathy", "Renal dysfunction", "Hepatomegaly",
                    "Short stature", "Tall stature", "Obesity", "Diabetes",
                    "Hypertension", "Anemia", "Thrombocytopenia", "Immunodeficiency",
                    "Skin abnormality", "Hair abnormality", "Dental abnormality",
                    "Skeletal abnormality", "Joint contracture", "Scoliosis"
                ]
            else:  # General/mouse phenotypes
                phenotype_names = [
                    "Embryonic lethality", "Postnatal lethality", "Growth retardation",
                    "Behavioral abnormality", "Neurological abnormality", 
                    "Cardiovascular abnormality", "Respiratory abnormality",
                    "Digestive system abnormality", "Reproductive system abnormality",
                    "Immune system abnormality", "Metabolic abnormality",
                    "Coat color abnormality", "Hair follicle abnormality",
                    "Bone abnormality", "Muscle abnormality", "Eye abnormality",
                    "Ear abnormality", "Kidney abnormality", "Liver abnormality"
                ]
            
            properties["name"] = random.choice(phenotype_names)

        # Generate description
        if PhenotypeAdapterPhenotypeField.DESCRIPTION in self.fields:
            phenotype_name = properties.get("name", "phenotype")
            properties["description"] = f"An observable characteristic or trait related to {phenotype_name.lower()}"

        # Generate synonyms
        if PhenotypeAdapterPhenotypeField.SYNONYMS in self.fields:
            if random.random() < 0.3:  # 30% have synonyms
                base_name = properties.get("name", "phenotype")
                synonyms = []
                if "abnormality" in base_name.lower():
                    synonyms.append(base_name.replace("abnormality", "anomaly"))
                if "disability" in base_name.lower():
                    synonyms.append(base_name.replace("disability", "impairment"))
                
                if synonyms:
                    properties["synonyms"] = "|".join(synonyms)

        # Set HPO ID for human phenotypes
        if PhenotypeAdapterPhenotypeField.HPO_ID in self.fields:
            if self.organism == "9606":
                properties["hpo_id"] = self.id

        # Set MP ID for mouse phenotypes
        if PhenotypeAdapterPhenotypeField.MAMMALIAN_PHENOTYPE_ID in self.fields:
            if self.organism == "10090":
                properties["mp_id"] = self.id

        # Set organism
        if PhenotypeAdapterPhenotypeField.ORGANISM in self.fields:
            properties["organism"] = self.organism

        # Generate phenotype type
        if PhenotypeAdapterPhenotypeField.PHENOTYPE_TYPE in self.fields:
            phenotype_types = [
                "morphological", "physiological", "behavioral", "developmental",
                "biochemical", "cellular", "molecular", "growth"
            ]
            properties["phenotype_type"] = random.choice(phenotype_types)

        # Generate severity
        if PhenotypeAdapterPhenotypeField.SEVERITY in self.fields:
            if random.random() < 0.6:  # 60% have severity info
                properties["severity"] = random.choice([
                    "mild", "moderate", "severe", "profound", "borderline"
                ])

        # Generate onset
        if PhenotypeAdapterPhenotypeField.ONSET in self.fields:
            if random.random() < 0.5:  # 50% have onset info
                if self.organism == "9606":  # Human
                    properties["onset"] = random.choice([
                        "congenital", "neonatal", "infantile", "childhood",
                        "juvenile", "adult", "late_adult"
                    ])
                else:  # Mouse/other
                    properties["onset"] = random.choice([
                        "embryonic", "postnatal", "weaning", "young_adult", "adult"
                    ])

        return properties