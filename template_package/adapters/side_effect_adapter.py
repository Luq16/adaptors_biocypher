from __future__ import annotations

import random
import string
from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class SideEffectAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    SIDE_EFFECT = auto()


class SideEffectAdapterSideEffectField(Enum):
    """
    Define possible fields the adapter can provide for side effects.
    """
    NAME = "name"
    SYNONYMS = "synonyms"
    MEDDRA_ID = "meddra_id"
    SEVERITY = "severity" 
    FREQUENCY_CLASS = "frequency_class"


class SideEffectAdapterEdgeType(Enum):
    """
    Enum for the types of edges the adapter can provide.
    """
    DRUG_HAS_SIDE_EFFECT = "drug_has_side_effect"
    SIDE_EFFECT_IS_SUBTYPE_OF_SIDE_EFFECT = "side_effect_is_subtype_of_side_effect"


class SideEffectAdapterDrugSideEffectEdgeField(Enum):
    """
    Define possible fields for drug-side effect edges.
    """
    SOURCE = "source"
    FREQUENCY = "frequency"
    PROPORTIONAL_REPORTING_RATIO = "proportional_reporting_ratio"
    CONFIDENCE_SCORE = "confidence_score"
    PUBMED_IDS = "pubmed_ids"
    DETECTION_METHOD = "detection_method"
    SEVERITY_SCORE = "severity_score"


class SideEffectAdapter:
    """
    Side Effect BioCypher adapter. Generates side effect nodes and drug-side effect 
    interaction edges for creating a knowledge graph of adverse drug reactions.

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
        logger.info("Generating side effect nodes.")

        self.nodes = []

        if SideEffectAdapterNodeType.SIDE_EFFECT in self.node_types:
            # Generate side effects
            num_side_effects = 30 if self.test_mode else 300
            [self.nodes.append(SideEffectNode(fields=self.node_fields)) 
             for _ in range(num_side_effects)]

        for node in self.nodes:
            yield (node.get_id(), node.get_label(), node.get_properties())

    def get_edges(self, target_drugs=None, target_side_effects=None):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.

        Args:
            target_drugs: List of drug nodes to create side effect relationships with.
            target_side_effects: List of side effect nodes for hierarchical relationships.
        """
        logger.info("Generating side effect edges.")

        if not self.nodes:
            list(self.get_nodes())  # Generate nodes first

        # Drug-side effect relationships
        if (SideEffectAdapterEdgeType.DRUG_HAS_SIDE_EFFECT in self.edge_types 
            and target_drugs):
            
            for drug in target_drugs:
                # Each drug has a 40% chance of having side effects
                if random.random() < 0.4:
                    # Select 1-4 random side effects for this drug
                    num_side_effects = random.randint(1, 4)
                    selected_side_effects = random.sample(self.nodes, min(num_side_effects, len(self.nodes)))
                    
                    for side_effect in selected_side_effects:
                        drug_id = drug[0] if isinstance(drug, tuple) else str(drug)
                        side_effect_id = side_effect.get_id()
                        
                        relationship_id = "".join(
                            random.choice(string.ascii_letters + string.digits)
                            for _ in range(12)
                        )
                        
                        properties = {}
                        for field in self.edge_fields:
                            if field == SideEffectAdapterDrugSideEffectEdgeField.SOURCE:
                                properties["source"] = random.choice([
                                    "SIDER", "OFFSIDES", "ADReCS", "FDA Adverse Event Reporting System",
                                    "EudraVigilance", "VigiBase", "Clinical Trials"
                                ])
                            elif field == SideEffectAdapterDrugSideEffectEdgeField.FREQUENCY:
                                properties["frequency"] = random.choice([
                                    "very common", "common", "uncommon", "rare", "very rare",
                                    ">10%", "1-10%", "0.1-1%", "0.01-0.1%", "<0.01%"
                                ])
                            elif field == SideEffectAdapterDrugSideEffectEdgeField.PROPORTIONAL_REPORTING_RATIO:
                                properties["proportional_reporting_ratio"] = round(random.uniform(0.5, 5.0), 3)
                            elif field == SideEffectAdapterDrugSideEffectEdgeField.CONFIDENCE_SCORE:
                                properties["confidence_score"] = round(random.uniform(0.1, 1.0), 3)
                            elif field == SideEffectAdapterDrugSideEffectEdgeField.PUBMED_IDS:
                                num_pubmed = random.randint(1, 3)
                                pubmed_ids = [str(random.randint(15000000, 35000000)) for _ in range(num_pubmed)]
                                properties["pubmed_ids"] = "|".join(pubmed_ids)
                            elif field == SideEffectAdapterDrugSideEffectEdgeField.DETECTION_METHOD:
                                properties["detection_method"] = random.choice([
                                    "spontaneous reporting", "clinical trial", "post-marketing surveillance",
                                    "electronic health records", "systematic review", "meta-analysis",
                                    "pharmacovigilance database mining", "case-control study"
                                ])
                            elif field == SideEffectAdapterDrugSideEffectEdgeField.SEVERITY_SCORE:
                                properties["severity_score"] = random.randint(1, 5)

                        yield (
                            relationship_id,
                            drug_id,
                            side_effect_id,
                            SideEffectAdapterEdgeType.DRUG_HAS_SIDE_EFFECT.value,
                            properties,
                        )

        # Side effect hierarchical relationships  
        if (SideEffectAdapterEdgeType.SIDE_EFFECT_IS_SUBTYPE_OF_SIDE_EFFECT in self.edge_types):
            # Create some hierarchical relationships between side effects
            for i, side_effect in enumerate(self.nodes):
                # 20% chance of having a parent side effect
                if random.random() < 0.2 and i > 0:
                    parent_side_effect = random.choice(self.nodes[:i])  # Choose from earlier nodes
                    
                    relationship_id = "".join(
                        random.choice(string.ascii_letters + string.digits)
                        for _ in range(12)
                    )
                    
                    yield (
                        relationship_id,
                        side_effect.get_id(),
                        parent_side_effect.get_id(),
                        SideEffectAdapterEdgeType.SIDE_EFFECT_IS_SUBTYPE_OF_SIDE_EFFECT.value,
                        {"relationship_type": "is_a"},
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
            self.node_types = [type for type in SideEffectAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [field for field in SideEffectAdapterSideEffectField]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in SideEffectAdapterEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in SideEffectAdapterDrugSideEffectEdgeField]


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


class SideEffectNode(Node):
    """
    Generates instances of side effect nodes.
    """

    def __init__(self, fields: Optional[list] = None):
        self.fields = fields
        self.id = self._generate_id()
        self.label = "side_effect"
        self.properties = self._generate_properties()

    def _generate_id(self):
        """
        Generate MedDRA-style identifier.
        """
        # MedDRA format: numeric ID (8 digits)
        return f"MEDDRA:{random.randint(10000000, 99999999)}"

    def _generate_properties(self):
        properties = {}

        # Generate side effect name
        if SideEffectAdapterSideEffectField.NAME in self.fields:
            side_effect_names = [
                "Nausea", "Headache", "Dizziness", "Fatigue", "Diarrhea",
                "Constipation", "Vomiting", "Dry mouth", "Drowsiness", "Insomnia",
                "Anxiety", "Depression", "Rash", "Itching", "Swelling",
                "Weight gain", "Weight loss", "Blurred vision", "Tinnitus", "Hair loss",
                "Muscle pain", "Joint pain", "Back pain", "Chest pain", "Abdominal pain",
                "Shortness of breath", "Cough", "Fever", "Chills", "Sweating",
                "Tremor", "Confusion", "Memory problems", "Concentration difficulties",
                "Loss of appetite", "Increased appetite", "Heartburn", "Bloating",
                "Hypertension", "Hypotension", "Tachycardia", "Bradycardia",
                "Arrhythmia", "Palpitations", "Edema", "Bruising", "Bleeding",
                "Anemia", "Thrombocytopenia", "Neutropenia", "Liver dysfunction"
            ]
            properties["name"] = random.choice(side_effect_names)

        # Generate synonyms
        if SideEffectAdapterSideEffectField.SYNONYMS in self.fields:
            if random.random() < 0.3:  # 30% have synonyms
                synonyms = [
                    "adverse effect", "adverse reaction", "side reaction", 
                    "unwanted effect", "unintended effect", "drug reaction"
                ]
                properties["synonyms"] = random.choice(synonyms)

        # Generate MedDRA ID
        if SideEffectAdapterSideEffectField.MEDDRA_ID in self.fields:
            properties["meddra_id"] = self.id.split(":")[1]

        # Generate severity
        if SideEffectAdapterSideEffectField.SEVERITY in self.fields:
            properties["severity"] = random.choice([
                "mild", "moderate", "severe", "life-threatening"
            ])

        # Generate frequency class
        if SideEffectAdapterSideEffectField.FREQUENCY_CLASS in self.fields:
            properties["frequency_class"] = random.choice([
                "very common", "common", "uncommon", "rare", "very rare"
            ])

        return properties