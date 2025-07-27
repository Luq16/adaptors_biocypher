from __future__ import annotations

import random
import string
from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class CompoundAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    COMPOUND = auto()


class CompoundAdapterCompoundField(Enum):
    """
    Define possible fields the adapter can provide for compounds.
    """
    ID = "id"
    NAME = "name"
    SMILES = "smiles"
    INCHI = "inchi"
    INCHI_KEY = "inchi_key"
    FORMULA = "formula"
    MOLECULAR_WEIGHT = "molecular_weight"
    DRUG_TYPE = "drug_type"
    DESCRIPTION = "description"


class CompoundAdapterEdgeType(Enum):
    """
    Enum for the types of edges the adapter can provide.
    """
    COMPOUND_TARGETS_PROTEIN = "compound_targets_protein"
    COMPOUND_TREATS_DISEASE = "compound_treats_disease"


class CompoundAdapterCompoundProteinEdgeField(Enum):
    """
    Define possible fields for compound-protein edges.
    """
    INTERACTION_TYPE = "interaction_type"
    ACTIVITY_TYPE = "activity_type"
    ACTIVITY_VALUE = "activity_value"
    CONFIDENCE_SCORE = "confidence_score"
    SOURCE = "source"


class CompoundAdapterCompoundDiseaseEdgeField(Enum):
    """
    Define possible fields for compound-disease edges.
    """
    INDICATION_TYPE = "indication_type"
    APPROVAL_STATUS = "approval_status"
    PHASE = "phase"
    SOURCE = "source"


class CompoundAdapter:
    """
    Compound BioCypher adapter. Generates compound nodes and compound-target 
    interaction edges for creating a knowledge graph.

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
        logger.info("Generating compound nodes.")

        self.nodes = []

        if CompoundAdapterNodeType.COMPOUND in self.node_types:
            # In a real implementation, this would fetch from ChEMBL, DrugBank, etc.
            # For now, generating example compounds
            num_compounds = 30 if self.test_mode else 500
            [self.nodes.append(Compound(fields=self.node_fields)) 
             for _ in range(num_compounds)]

        for node in self.nodes:
            yield (node.get_id(), node.get_label(), node.get_properties())

    def get_edges(self, target_proteins=None, target_diseases=None):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.

        Args:
            target_proteins: List of protein nodes to create edges with.
            target_diseases: List of disease nodes to create edges with.
        """
        logger.info("Generating compound interaction edges.")

        if not self.nodes:
            raise ValueError("No nodes found. Please run get_nodes() first.")

        # Create compound-protein edges
        if (CompoundAdapterEdgeType.COMPOUND_TARGETS_PROTEIN in self.edge_types 
            and target_proteins):
            
            for compound in self.nodes:
                # Each compound targets 1-3 proteins
                num_targets = random.randint(1, 3)
                selected_proteins = random.sample(
                    target_proteins, 
                    min(num_targets, len(target_proteins))
                )
                
                for protein in selected_proteins:
                    relationship_id = "".join(
                        random.choice(string.ascii_letters + string.digits)
                        for _ in range(10)
                    )
                    
                    properties = {}
                    for field in self.edge_fields:
                        if field == CompoundAdapterCompoundProteinEdgeField.INTERACTION_TYPE:
                            properties["interaction_type"] = random.choice([
                                "inhibitor", "activator", "antagonist", "agonist", "modulator"
                            ])
                        elif field == CompoundAdapterCompoundProteinEdgeField.ACTIVITY_TYPE:
                            properties["activity_type"] = random.choice([
                                "IC50", "EC50", "Ki", "Kd", "ED50"
                            ])
                        elif field == CompoundAdapterCompoundProteinEdgeField.ACTIVITY_VALUE:
                            properties["activity_value"] = round(random.uniform(0.1, 1000.0), 3)
                        elif field == CompoundAdapterCompoundProteinEdgeField.CONFIDENCE_SCORE:
                            properties["confidence_score"] = round(random.uniform(0.1, 1.0), 3)
                        elif field == CompoundAdapterCompoundProteinEdgeField.SOURCE:
                            properties["source"] = random.choice([
                                "ChEMBL", "DrugBank", "BindingDB", "PubChem"
                            ])

                    yield (
                        relationship_id,
                        compound.get_id(),
                        protein.get_id() if hasattr(protein, 'get_id') else protein,
                        CompoundAdapterEdgeType.COMPOUND_TARGETS_PROTEIN.value,
                        properties,
                    )

        # Create compound-disease edges
        if (CompoundAdapterEdgeType.COMPOUND_TREATS_DISEASE in self.edge_types 
            and target_diseases):
            
            for compound in self.nodes:
                # Each compound treats 1-2 diseases
                num_indications = random.randint(1, 2)
                selected_diseases = random.sample(
                    target_diseases, 
                    min(num_indications, len(target_diseases))
                )
                
                for disease in selected_diseases:
                    relationship_id = "".join(
                        random.choice(string.ascii_letters + string.digits)
                        for _ in range(10)
                    )
                    
                    properties = {}
                    for field in self.edge_fields:
                        if field == CompoundAdapterCompoundDiseaseEdgeField.INDICATION_TYPE:
                            properties["indication_type"] = random.choice([
                                "treatment", "prevention", "symptom_relief", "cure"
                            ])
                        elif field == CompoundAdapterCompoundDiseaseEdgeField.APPROVAL_STATUS:
                            properties["approval_status"] = random.choice([
                                "approved", "investigational", "experimental", "withdrawn"
                            ])
                        elif field == CompoundAdapterCompoundDiseaseEdgeField.PHASE:
                            properties["phase"] = random.choice([
                                "Phase I", "Phase II", "Phase III", "Phase IV", "Approved"
                            ])
                        elif field == CompoundAdapterCompoundDiseaseEdgeField.SOURCE:
                            properties["source"] = random.choice([
                                "DrugBank", "ClinicalTrials.gov", "FDA", "EMA"
                            ])

                    yield (
                        relationship_id,
                        compound.get_id(),
                        disease.get_id() if hasattr(disease, 'get_id') else disease,
                        CompoundAdapterEdgeType.COMPOUND_TREATS_DISEASE.value,
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
            self.node_types = [type for type in CompoundAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [field for field in CompoundAdapterCompoundField]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in CompoundAdapterEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [
                field for field in chain(
                    CompoundAdapterCompoundProteinEdgeField,
                    CompoundAdapterCompoundDiseaseEdgeField,
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


class Compound(Node):
    """
    Generates instances of chemical compounds.
    """

    def __init__(self, fields: Optional[list] = None):
        self.fields = fields
        self.id = self._generate_id()
        self.label = "compound"
        self.properties = self._generate_properties()

    def _generate_id(self):
        """
        Generate a random compound ID (ChEMBL-style).
        """
        return f"CHEMBL{random.randint(1000, 999999)}"

    def _generate_properties(self):
        properties = {}

        # Generate compound name
        if CompoundAdapterCompoundField.NAME in self.fields:
            compound_names = [
                "Aspirin", "Ibuprofen", "Paracetamol", "Metformin", "Atorvastatin",
                "Lisinopril", "Amlodipine", "Simvastatin", "Omeprazole", "Levothyroxine",
                "Azithromycin", "Amoxicillin", "Hydrochlorothiazide", "Metoprolol",
                "Losartan", "Albuterol", "Furosemide", "Prednisone", "Tramadol", "Sertraline"
            ]
            properties["name"] = random.choice(compound_names)

        # Generate SMILES
        if CompoundAdapterCompoundField.SMILES in self.fields:
            # Simplified SMILES examples
            smiles_examples = [
                "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
                "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  # Ibuprofen
                "CC(=O)NC1=CC=C(C=C1)O",  # Paracetamol
                "CN(C)C(=N)N=C(N)N",  # Metformin
                "COC1=CC=C(C=C1)CCN",  # Example
            ]
            properties["smiles"] = random.choice(smiles_examples)

        # Generate InChI
        if CompoundAdapterCompoundField.INCHI in self.fields:
            properties["inchi"] = f"InChI=1S/C{random.randint(6,20)}H{random.randint(4,30)}N{random.randint(0,4)}O{random.randint(1,6)}/..."

        # Generate InChI Key
        if CompoundAdapterCompoundField.INCHI_KEY in self.fields:
            # Generate a random InChI Key format
            key1 = "".join(random.choices(string.ascii_uppercase, k=14))
            key2 = "".join(random.choices(string.ascii_uppercase, k=10))
            key3 = "".join(random.choices(string.ascii_uppercase + string.digits + "-", k=1))
            properties["inchi_key"] = f"{key1}-{key2}-{key3}"

        # Generate molecular formula
        if CompoundAdapterCompoundField.FORMULA in self.fields:
            c_count = random.randint(6, 25)
            h_count = random.randint(4, 40)
            n_count = random.randint(0, 4)
            o_count = random.randint(1, 8)
            
            formula = f"C{c_count}H{h_count}"
            if n_count > 0:
                formula += f"N{n_count}" if n_count > 1 else "N"
            if o_count > 0:
                formula += f"O{o_count}" if o_count > 1 else "O"
            
            properties["formula"] = formula

        # Generate molecular weight
        if CompoundAdapterCompoundField.MOLECULAR_WEIGHT in self.fields:
            properties["molecular_weight"] = round(random.uniform(100.0, 800.0), 2)

        # Generate drug type
        if CompoundAdapterCompoundField.DRUG_TYPE in self.fields:
            drug_types = [
                "small molecule", "biologics", "peptide", "antibody", 
                "vaccine", "gene therapy", "antisense oligonucleotide"
            ]
            properties["drug_type"] = random.choice(drug_types)

        # Generate description
        if CompoundAdapterCompoundField.DESCRIPTION in self.fields:
            descriptions = [
                "Anti-inflammatory agent", "Analgesic compound", "Antihypertensive drug",
                "Antibiotic agent", "Antidiabetic medication", "Cholesterol-lowering drug",
                "Bronchodilator", "Proton pump inhibitor", "Antidepressant", "Diuretic"
            ]
            properties["description"] = random.choice(descriptions)

        return properties