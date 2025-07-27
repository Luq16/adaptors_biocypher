from __future__ import annotations

import random
import string
from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class DrugAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    DRUG = auto()


class DrugAdapterDrugField(Enum):
    """
    Define possible fields the adapter can provide for drugs.
    """
    NAME = "name"
    SMILES = "smiles"
    INCHI = "inchi"
    INCHI_KEY = "inchi_key"
    MOLECULAR_WEIGHT = "molecular_weight"
    CAS_NUMBER = "cas_number"
    DRUGBANK_ID = "drugbank_id"
    CHEMBL_ID = "chembl_id"
    PUBCHEM_CID = "pubchem_cid"
    ATC_CODES = "atc_codes"
    DRUG_GROUPS = "drug_groups"
    INDICATION = "indication"
    MECHANISM_OF_ACTION = "mechanism_of_action"
    PHARMACOLOGY = "pharmacology"
    TOXICITY = "toxicity"
    CATEGORIES = "categories"


class DrugAdapterEdgeType(Enum):
    """
    Enum for the types of edges the adapter can provide.
    """
    DRUG_TARGETS_PROTEIN = "drug_targets_protein"
    DRUG_TREATS_DISEASE = "drug_treats_disease"
    DRUG_DRUG_INTERACTION = "drug_drug_interaction"
    DRUG_HAS_SIDE_EFFECT = "drug_has_side_effect"
    DRUG_PARTICIPATES_IN_PATHWAY = "drug_participates_in_pathway"


class DrugAdapterDrugTargetEdgeField(Enum):
    """
    Define possible fields for drug-target edges.
    """
    INTERACTION_TYPE = "interaction_type"
    BINDING_AFFINITY = "binding_affinity"
    IC50 = "ic50"
    KD = "kd"
    EC50 = "ec50"
    ACTIVITY = "activity"
    SOURCE_DATABASE = "source_database"
    PUBMED_IDS = "pubmed_ids"
    ASSAY_TYPE = "assay_type"
    ORGANISM = "organism"


class DrugAdapterDrugDrugEdgeField(Enum):
    """
    Define possible fields for drug-drug interaction edges.
    """
    INTERACTION_TYPE = "interaction_type"
    SEVERITY = "severity"
    MECHANISM = "mechanism"
    CLINICAL_SIGNIFICANCE = "clinical_significance"
    SOURCE_DATABASE = "source_database"
    PUBMED_IDS = "pubmed_ids"


class DrugAdapter:
    """
    Drug BioCypher adapter. Generates drug nodes and drug-related interactions 
    for creating comprehensive pharmaceutical knowledge graphs.

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
        logger.info("Generating drug nodes.")

        self.nodes = []

        if DrugAdapterNodeType.DRUG in self.node_types:
            # Generate drugs
            num_drugs = 50 if self.test_mode else 500
            [self.nodes.append(DrugNode(fields=self.node_fields)) 
             for _ in range(num_drugs)]

        for node in self.nodes:
            yield (node.get_id(), node.get_label(), node.get_properties())

    def get_edges(self, target_proteins=None, target_diseases=None, target_side_effects=None, target_pathways=None):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.

        Args:
            target_proteins: List of protein nodes for drug-target interactions.
            target_diseases: List of disease nodes for drug-disease relationships.
            target_side_effects: List of side effect nodes for adverse reactions.
            target_pathways: List of pathway nodes for drug pathway participation.
        """
        logger.info("Generating drug interaction edges.")

        if not self.nodes:
            list(self.get_nodes())  # Generate nodes first

        # Drug-protein target interactions
        if (DrugAdapterEdgeType.DRUG_TARGETS_PROTEIN in self.edge_types 
            and target_proteins):
            
            for drug in self.nodes:
                # Each drug targets 1-5 proteins
                if random.random() < 0.7:  # 70% of drugs have targets
                    num_targets = random.randint(1, 5)
                    selected_proteins = random.sample(target_proteins, min(num_targets, len(target_proteins)))
                    
                    for protein in selected_proteins:
                        protein_id = protein[0] if isinstance(protein, tuple) else str(protein)
                        
                        relationship_id = "".join(
                            random.choice(string.ascii_letters + string.digits)
                            for _ in range(12)
                        )
                        
                        properties = {}
                        for field in self.edge_fields:
                            if field == DrugAdapterDrugTargetEdgeField.INTERACTION_TYPE:
                                properties["interaction_type"] = random.choice([
                                    "inhibitor", "agonist", "antagonist", "modulator",
                                    "activator", "inducer", "blocker", "substrate",
                                    "cofactor", "allosteric modulator"
                                ])
                            elif field == DrugAdapterDrugTargetEdgeField.BINDING_AFFINITY:
                                properties["binding_affinity"] = f"{random.uniform(0.1, 1000):.2f} nM"
                            elif field == DrugAdapterDrugTargetEdgeField.IC50:
                                properties["ic50"] = f"{random.uniform(0.01, 100):.3f} μM"
                            elif field == DrugAdapterDrugTargetEdgeField.KD:
                                properties["kd"] = f"{random.uniform(0.1, 500):.2f} nM"
                            elif field == DrugAdapterDrugTargetEdgeField.EC50:
                                properties["ec50"] = f"{random.uniform(0.01, 50):.3f} μM"
                            elif field == DrugAdapterDrugTargetEdgeField.ACTIVITY:
                                properties["activity"] = random.choice([
                                    "active", "inactive", "inconclusive", "moderate", "high", "low"
                                ])
                            elif field == DrugAdapterDrugTargetEdgeField.SOURCE_DATABASE:
                                properties["source_database"] = random.choice([
                                    "ChEMBL", "DrugBank", "STITCH", "BindingDB", 
                                    "PubChem", "DGIDB", "Pharos", "CTD"
                                ])
                            elif field == DrugAdapterDrugTargetEdgeField.PUBMED_IDS:
                                num_pubmed = random.randint(1, 4)
                                pubmed_ids = [str(random.randint(15000000, 35000000)) for _ in range(num_pubmed)]
                                properties["pubmed_ids"] = "|".join(pubmed_ids)
                            elif field == DrugAdapterDrugTargetEdgeField.ASSAY_TYPE:
                                properties["assay_type"] = random.choice([
                                    "biochemical", "cell-based", "functional", "binding",
                                    "phenotypic", "reporter gene", "enzymatic"
                                ])
                            elif field == DrugAdapterDrugTargetEdgeField.ORGANISM:
                                properties["organism"] = random.choice([
                                    "Homo sapiens", "Mus musculus", "Rattus norvegicus"
                                ])

                        yield (
                            relationship_id,
                            drug.get_id(),
                            protein_id,
                            DrugAdapterEdgeType.DRUG_TARGETS_PROTEIN.value,
                            properties,
                        )

        # Drug-disease treatment relationships
        if (DrugAdapterEdgeType.DRUG_TREATS_DISEASE in self.edge_types 
            and target_diseases):
            
            for drug in self.nodes:
                # Each drug treats 1-3 diseases
                if random.random() < 0.6:  # 60% of drugs have disease indications
                    num_diseases = random.randint(1, 3)
                    selected_diseases = random.sample(target_diseases, min(num_diseases, len(target_diseases)))
                    
                    for disease in selected_diseases:
                        disease_id = disease[0] if isinstance(disease, tuple) else str(disease)
                        
                        relationship_id = "".join(
                            random.choice(string.ascii_letters + string.digits)
                            for _ in range(12)
                        )
                        
                        properties = {
                            "indication_type": random.choice([
                                "approved", "investigational", "off-label", "experimental"
                            ]),
                            "therapeutic_class": random.choice([
                                "antineoplastic", "anti-inflammatory", "analgesic",
                                "antibiotic", "antiviral", "antidepressant", "antihypertensive"
                            ]),
                            "approval_status": random.choice([
                                "FDA approved", "EMA approved", "investigational", "withdrawn"
                            ])
                        }

                        yield (
                            relationship_id,
                            drug.get_id(),
                            disease_id,
                            DrugAdapterEdgeType.DRUG_TREATS_DISEASE.value,
                            properties,
                        )

        # Drug-drug interactions
        if DrugAdapterEdgeType.DRUG_DRUG_INTERACTION in self.edge_types:
            interactions_created = set()
            
            for i, drug_a in enumerate(self.nodes):
                for drug_b in self.nodes[i+1:]:  # Avoid duplicates
                    if random.random() < 0.05:  # 5% chance of interaction
                        
                        # Create unique interaction identifier
                        pair = tuple(sorted([drug_a.get_id(), drug_b.get_id()]))
                        
                        if pair not in interactions_created:
                            interactions_created.add(pair)
                            
                            relationship_id = "".join(
                                random.choice(string.ascii_letters + string.digits)
                                for _ in range(12)
                            )
                            
                            properties = {}
                            for field in self.edge_fields:
                                if field == DrugAdapterDrugDrugEdgeField.INTERACTION_TYPE:
                                    properties["interaction_type"] = random.choice([
                                        "pharmacokinetic", "pharmacodynamic", "synergistic",
                                        "antagonistic", "additive", "competitive"
                                    ])
                                elif field == DrugAdapterDrugDrugEdgeField.SEVERITY:
                                    properties["severity"] = random.choice([
                                        "minor", "moderate", "major", "contraindicated"
                                    ])
                                elif field == DrugAdapterDrugDrugEdgeField.MECHANISM:
                                    properties["mechanism"] = random.choice([
                                        "CYP450 inhibition", "CYP450 induction", "protein binding displacement",
                                        "renal elimination competition", "receptor competition",
                                        "metabolic pathway interference"
                                    ])
                                elif field == DrugAdapterDrugDrugEdgeField.CLINICAL_SIGNIFICANCE:
                                    properties["clinical_significance"] = random.choice([
                                        "monitor therapy", "consider therapy modification",
                                        "avoid combination", "dose adjustment required"
                                    ])
                                elif field == DrugAdapterDrugDrugEdgeField.SOURCE_DATABASE:
                                    properties["source_database"] = random.choice([
                                        "DrugBank", "DDinter", "SIDER", "Lexicomp", "Micromedex"
                                    ])

                            yield (
                                relationship_id,
                                drug_a.get_id(),
                                drug_b.get_id(),
                                DrugAdapterEdgeType.DRUG_DRUG_INTERACTION.value,
                                properties,
                            )

        # Drug-side effect relationships
        if (DrugAdapterEdgeType.DRUG_HAS_SIDE_EFFECT in self.edge_types 
            and target_side_effects):
            
            for drug in self.nodes:
                # Each drug has 1-6 side effects
                if random.random() < 0.8:  # 80% of drugs have side effects
                    num_side_effects = random.randint(1, 6)
                    selected_side_effects = random.sample(target_side_effects, min(num_side_effects, len(target_side_effects)))
                    
                    for side_effect in selected_side_effects:
                        side_effect_id = side_effect[0] if isinstance(side_effect, tuple) else str(side_effect)
                        
                        relationship_id = "".join(
                            random.choice(string.ascii_letters + string.digits)
                            for _ in range(12)
                        )
                        
                        properties = {
                            "frequency": random.choice([
                                "very common", "common", "uncommon", "rare", "very rare"
                            ]),
                            "severity": random.choice(["mild", "moderate", "severe"]),
                            "source": random.choice(["SIDER", "FDA AERS", "clinical trials"])
                        }

                        yield (
                            relationship_id,
                            drug.get_id(),
                            side_effect_id,
                            DrugAdapterEdgeType.DRUG_HAS_SIDE_EFFECT.value,
                            properties,
                        )

        # Drug-pathway participation
        if (DrugAdapterEdgeType.DRUG_PARTICIPATES_IN_PATHWAY in self.edge_types 
            and target_pathways):
            
            for drug in self.nodes:
                # Each drug participates in 1-3 pathways
                if random.random() < 0.5:  # 50% of drugs have pathway associations
                    num_pathways = random.randint(1, 3)
                    selected_pathways = random.sample(target_pathways, min(num_pathways, len(target_pathways)))
                    
                    for pathway in selected_pathways:
                        pathway_id = pathway[0] if isinstance(pathway, tuple) else str(pathway)
                        
                        relationship_id = "".join(
                            random.choice(string.ascii_letters + string.digits)
                            for _ in range(12)
                        )
                        
                        properties = {
                            "role": random.choice([
                                "inhibitor", "activator", "modulator", "substrate", "cofactor"
                            ]),
                            "pathway_type": random.choice([
                                "metabolic", "signaling", "transport", "degradation"
                            ])
                        }

                        yield (
                            relationship_id,
                            drug.get_id(),
                            pathway_id,
                            DrugAdapterEdgeType.DRUG_PARTICIPATES_IN_PATHWAY.value,
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
            self.node_types = [type for type in DrugAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [field for field in DrugAdapterDrugField]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in DrugAdapterEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in DrugAdapterDrugTargetEdgeField]


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


class DrugNode(Node):
    """
    Generates instances of drug nodes.
    """

    def __init__(self, fields: Optional[list] = None):
        self.fields = fields
        self.id = self._generate_id()
        self.label = "drug"
        self.properties = self._generate_properties()

    def _generate_id(self):
        """
        Generate a DrugBank-style identifier.
        """
        # DrugBank format: DB + 5 digits
        return f"DB{random.randint(10000, 99999):05d}"

    def _generate_properties(self):
        properties = {}

        # Generate drug name
        if DrugAdapterDrugField.NAME in self.fields:
            drug_names = [
                "Aspirin", "Ibuprofen", "Acetaminophen", "Metformin", "Atorvastatin",
                "Lisinopril", "Amlodipine", "Metoprolol", "Omeprazole", "Simvastatin",
                "Levothyroxine", "Azithromycin", "Furosemide", "Hydrochlorothiazide", "Losartan",
                "Gabapentin", "Tramadol", "Prednisone", "Amoxicillin", "Clopidogrel",
                "Warfarin", "Insulin", "Morphine", "Digoxin", "Fluoxetine",
                "Sertraline", "Citalopram", "Escitalopram", "Venlafaxine", "Duloxetine",
                "Donepezil", "Memantine", "Levodopa", "Carbidopa", "Quetiapine",
                "Risperidone", "Olanzapine", "Haloperidol", "Lithium", "Valproate",
                "Phenytoin", "Carbamazepine", "Lamotrigine", "Topiramate", "Levetiracetam"
            ]
            properties["name"] = random.choice(drug_names)

        # Generate SMILES
        if DrugAdapterDrugField.SMILES in self.fields:
            # Simplified SMILES patterns for common functional groups
            smiles_patterns = [
                "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin-like
                "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  # Ibuprofen-like
                "CC(=O)NC1=CC=C(O)C=C1",  # Acetaminophen-like
                "CN(C)C(=N)NC(=N)N",  # Metformin-like
                "CCC1=CC=C(C=C1)N(CC)CC",  # Generic aromatic
                "CC1=CC(=CC=C1)C(=O)O",  # Benzoic acid derivative
                "CCOC1=CC=C(C=C1)C(=O)O",  # Ethyl ester derivative
                "CC(C)(C)NCC(O)C1=CC(O)=CC=C1"  # Beta-blocker like
            ]
            properties["smiles"] = random.choice(smiles_patterns)

        # Generate molecular weight
        if DrugAdapterDrugField.MOLECULAR_WEIGHT in self.fields:
            properties["molecular_weight"] = round(random.uniform(150.0, 800.0), 2)

        # Generate DrugBank ID
        if DrugAdapterDrugField.DRUGBANK_ID in self.fields:
            properties["drugbank_id"] = self.id

        # Generate ChEMBL ID
        if DrugAdapterDrugField.CHEMBL_ID in self.fields:
            properties["chembl_id"] = f"CHEMBL{random.randint(100000, 999999)}"

        # Generate PubChem CID
        if DrugAdapterDrugField.PUBCHEM_CID in self.fields:
            properties["pubchem_cid"] = f"CID:{random.randint(1000, 999999)}"

        # Generate ATC codes
        if DrugAdapterDrugField.ATC_CODES in self.fields:
            atc_codes = [
                "A02BC01", "A10BA02", "B01AC06", "C01AA05", "C03AA03",
                "C09AA02", "J01MA12", "N02AA01", "N05AH04", "R03BB01"
            ]
            properties["atc_codes"] = random.choice(atc_codes)

        # Generate drug groups
        if DrugAdapterDrugField.DRUG_GROUPS in self.fields:
            properties["drug_groups"] = random.choice([
                "approved", "investigational", "experimental", "illicit", "withdrawn"
            ])

        # Generate indication
        if DrugAdapterDrugField.INDICATION in self.fields:
            indications = [
                "Pain relief", "Hypertension", "Diabetes", "Depression", "Anxiety",
                "Cancer", "Infection", "Inflammation", "Cardiovascular disease", "Epilepsy"
            ]
            properties["indication"] = random.choice(indications)

        # Generate mechanism of action
        if DrugAdapterDrugField.MECHANISM_OF_ACTION in self.fields:
            mechanisms = [
                "COX-1/COX-2 inhibitor", "Calcium channel blocker", "ACE inhibitor",
                "Beta-adrenergic receptor antagonist", "Proton pump inhibitor",
                "SSRI", "SNRI", "Dopamine receptor antagonist", "GABA receptor modulator"
            ]
            properties["mechanism_of_action"] = random.choice(mechanisms)

        # Generate InChI Key
        if DrugAdapterDrugField.INCHI_KEY in self.fields:
            # Generate realistic-looking InChI Key
            part1 = "".join(random.choices("ABCDEFGHIJKLMNOPQRSTUVWXYZ", k=14))
            part2 = "".join(random.choices("ABCDEFGHIJKLMNOPQRSTUVWXYZ", k=10))
            part3 = "".join(random.choices("ABCDEFGHIJKLMNOP", k=1))
            properties["inchi_key"] = f"{part1}-{part2}-{part3}"

        # Generate categories
        if DrugAdapterDrugField.CATEGORIES in self.fields:
            categories = [
                "Analgesics", "Anti-inflammatory Agents", "Antihypertensive Agents",
                "Antidiabetic Agents", "Antidepressants", "Anticonvulsants",
                "Antimicrobials", "Cardiovascular Agents", "Central Nervous System Agents"
            ]
            properties["categories"] = random.choice(categories)

        return properties