from __future__ import annotations

import random
import string
from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class ECAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    EC_NUMBER = auto()


class ECAdapterECField(Enum):
    """
    Define possible fields the adapter can provide for EC numbers.
    """
    NAME = "name"
    DESCRIPTION = "description"
    REACTION = "reaction"
    COFACTORS = "cofactors"
    COMMENTS = "comments"
    PROSITE = "prosite"
    PATHWAY = "pathway"
    CLASS_NAME = "class_name"
    SUBCLASS_NAME = "subclass_name"
    SUBSUBCLASS_NAME = "subsubclass_name"


class ECAdapterEdgeType(Enum):
    """
    Enum for the types of edges the adapter can provide.
    """
    PROTEIN_HAS_EC_ANNOTATION = "protein_has_ec_annotation"
    EC_IS_A_EC = "ec_is_a_ec"
    EC_PART_OF_EC = "ec_part_of_ec"


class ECAdapterProteinECEdgeField(Enum):
    """
    Define possible fields for protein-EC annotation edges.
    """
    SOURCE = "source"
    EVIDENCE_CODE = "evidence_code"
    CONFIDENCE_SCORE = "confidence_score"
    PUBMED_IDS = "pubmed_ids"
    ASSIGNED_BY = "assigned_by"
    ANNOTATION_DATE = "annotation_date"
    QUALIFIERS = "qualifiers"


class ECAdapter:
    """
    EC (Enzyme Commission) BioCypher adapter. Generates EC number nodes and 
    protein-EC annotation edges for enzyme classification in knowledge graphs.

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
        logger.info("Generating EC number nodes.")

        self.nodes = []

        if ECAdapterNodeType.EC_NUMBER in self.node_types:
            # Generate EC numbers
            num_ec_numbers = 40 if self.test_mode else 400
            [self.nodes.append(ECNumberNode(fields=self.node_fields)) 
             for _ in range(num_ec_numbers)]

        for node in self.nodes:
            yield (node.get_id(), node.get_label(), node.get_properties())

    def get_edges(self, target_proteins=None):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.

        Args:
            target_proteins: List of protein nodes to create EC annotations for.
        """
        logger.info("Generating EC annotation edges.")

        if not self.nodes:
            list(self.get_nodes())  # Generate nodes first

        # Protein-EC annotation relationships
        if (ECAdapterEdgeType.PROTEIN_HAS_EC_ANNOTATION in self.edge_types 
            and target_proteins):
            
            for protein in target_proteins:
                # Each protein has a 30% chance of having EC annotations
                if random.random() < 0.3:
                    # Select 1-3 random EC numbers for this protein
                    num_ecs = random.randint(1, 3)
                    selected_ecs = random.sample(self.nodes, min(num_ecs, len(self.nodes)))
                    
                    for ec_number in selected_ecs:
                        protein_id = protein[0] if isinstance(protein, tuple) else str(protein)
                        
                        relationship_id = "".join(
                            random.choice(string.ascii_letters + string.digits)
                            for _ in range(12)
                        )
                        
                        properties = {}
                        for field in self.edge_fields:
                            if field == ECAdapterProteinECEdgeField.SOURCE:
                                properties["source"] = random.choice([
                                    "UniProt", "ExPASy", "KEGG", "BRENDA", "MetaCyc",
                                    "Reactome", "IntEnz", "ENZYME"
                                ])
                            elif field == ECAdapterProteinECEdgeField.EVIDENCE_CODE:
                                properties["evidence_code"] = random.choice([
                                    "IDA", "IMP", "IGI", "IPI", "IEP", "ISS", "ISA", 
                                    "ISO", "ISM", "IGC", "IBA", "IBD", "IKR", "IRD"
                                ])
                            elif field == ECAdapterProteinECEdgeField.CONFIDENCE_SCORE:
                                properties["confidence_score"] = round(random.uniform(0.5, 1.0), 3)
                            elif field == ECAdapterProteinECEdgeField.PUBMED_IDS:
                                num_pubmed = random.randint(1, 3)
                                pubmed_ids = [str(random.randint(10000000, 35000000)) for _ in range(num_pubmed)]
                                properties["pubmed_ids"] = "|".join(pubmed_ids)
                            elif field == ECAdapterProteinECEdgeField.ASSIGNED_BY:
                                properties["assigned_by"] = random.choice([
                                    "UniProt", "ExPASy", "KEGG", "InterPro", "Pfam"
                                ])
                            elif field == ECAdapterProteinECEdgeField.ANNOTATION_DATE:
                                year = random.randint(2010, 2024)
                                month = random.randint(1, 12)
                                day = random.randint(1, 28)
                                properties["annotation_date"] = f"{year:04d}-{month:02d}-{day:02d}"
                            elif field == ECAdapterProteinECEdgeField.QUALIFIERS:
                                if random.random() < 0.2:  # 20% chance
                                    properties["qualifiers"] = random.choice([
                                        "enables", "contributes_to", "colocalizes_with"
                                    ])

                        yield (
                            relationship_id,
                            protein_id,
                            ec_number.get_id(),
                            ECAdapterEdgeType.PROTEIN_HAS_EC_ANNOTATION.value,
                            properties,
                        )

        # EC hierarchical relationships (is_a)
        if ECAdapterEdgeType.EC_IS_A_EC in self.edge_types:
            for i, ec_child in enumerate(self.nodes):
                # 15% chance of having a parent EC
                if random.random() < 0.15 and i > 0:
                    # Find a suitable parent (more general EC number)
                    child_parts = ec_child.get_id().replace("EC:", "").split(".")
                    
                    # Try to create hierarchical relationship
                    if len(child_parts) == 4:  # Full EC number
                        potential_parents = []
                        for parent_ec in self.nodes[:i]:
                            parent_parts = parent_ec.get_id().replace("EC:", "").split(".")
                            # Parent should be more general (fewer specific parts)
                            if (len(parent_parts) < 4 and 
                                all(child_parts[j] == parent_parts[j] for j in range(len(parent_parts)))):
                                potential_parents.append(parent_ec)
                        
                        if potential_parents:
                            parent_ec = random.choice(potential_parents)
                            
                            relationship_id = "".join(
                                random.choice(string.ascii_letters + string.digits)
                                for _ in range(12)
                            )
                            
                            yield (
                                relationship_id,
                                ec_child.get_id(),
                                parent_ec.get_id(),
                                ECAdapterEdgeType.EC_IS_A_EC.value,
                                {"relationship_type": "is_a"},
                            )

        # EC part_of relationships
        if ECAdapterEdgeType.EC_PART_OF_EC in self.edge_types:
            for i, ec_part in enumerate(self.nodes):
                # 10% chance of being part of another EC
                if random.random() < 0.1 and i > 0:
                    parent_ec = random.choice(self.nodes[:i])
                    
                    relationship_id = "".join(
                        random.choice(string.ascii_letters + string.digits)
                        for _ in range(12)
                    )
                    
                    yield (
                        relationship_id,
                        ec_part.get_id(),
                        parent_ec.get_id(),
                        ECAdapterEdgeType.EC_PART_OF_EC.value,
                        {"relationship_type": "part_of"},
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
            self.node_types = [type for type in ECAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [field for field in ECAdapterECField]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in ECAdapterEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in ECAdapterProteinECEdgeField]


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


class ECNumberNode(Node):
    """
    Generates instances of EC number nodes.
    """

    def __init__(self, fields: Optional[list] = None):
        self.fields = fields
        self.id = self._generate_id()
        self.label = "ec_number"
        self.properties = self._generate_properties()

    def _generate_id(self):
        """
        Generate an EC number identifier.
        """
        # EC format: EC:1.2.3.4 (class.subclass.subsubclass.serial)
        ec_class = random.randint(1, 7)  # 7 main enzyme classes
        subclass = random.randint(1, 99)
        subsubclass = random.randint(1, 99)
        serial = random.randint(1, 999)
        
        return f"EC:{ec_class}.{subclass}.{subsubclass}.{serial}"

    def _generate_properties(self):
        properties = {}

        # Extract EC parts for classification
        ec_parts = self.id.replace("EC:", "").split(".")
        ec_class = int(ec_parts[0])

        # Generate enzyme name
        if ECAdapterECField.NAME in self.fields:
            # Base names on EC class
            class_enzymes = {
                1: ["oxidase", "dehydrogenase", "reductase", "peroxidase", "catalase"],
                2: ["transferase", "transaminase", "kinase", "acetyltransferase", "methyltransferase"],
                3: ["hydrolase", "peptidase", "protease", "lipase", "esterase"],
                4: ["lyase", "decarboxylase", "synthase", "cyclase", "dehydratase"],
                5: ["isomerase", "mutase", "racemase", "epimerase", "tautomerase"],
                6: ["ligase", "synthetase", "carboxylase", "synthase", "coenzyme A ligase"],
                7: ["translocase", "transporter", "channel", "pump", "carrier"]
            }
            
            base_name = random.choice(class_enzymes.get(ec_class, ["enzyme"]))
            substrate = random.choice([
                "glucose", "ATP", "NADH", "protein", "DNA", "RNA", "lipid",
                "amino acid", "fatty acid", "steroid", "nucleotide", "phosphate"
            ])
            properties["name"] = f"{substrate} {base_name}"

        # Generate description
        if ECAdapterECField.DESCRIPTION in self.fields:
            class_descriptions = {
                1: "Catalyzes oxidation-reduction reactions",
                2: "Catalyzes the transfer of a group from one compound to another",
                3: "Catalyzes the hydrolysis of various bonds",
                4: "Catalyzes the addition or removal of groups to form double bonds",
                5: "Catalyzes isomerization reactions",
                6: "Catalyzes the formation of bonds with ATP cleavage",
                7: "Catalyzes the movement of ions or molecules across membranes"
            }
            properties["description"] = class_descriptions.get(ec_class, "Enzyme activity")

        # Generate reaction
        if ECAdapterECField.REACTION in self.fields:
            reactions = [
                "A + NAD+ = B + NADH + H+",
                "ATP + A = ADP + A-phosphate",
                "A-B + H2O = A + B",
                "A = B + C",
                "A = B (isomerization)",
                "A + B + ATP = A-B + ADP + Pi"
            ]
            properties["reaction"] = random.choice(reactions)

        # Generate cofactors
        if ECAdapterECField.COFACTORS in self.fields:
            cofactors = [
                "Mg2+", "Zn2+", "Fe2+", "Ca2+", "Mn2+", "Cu2+", "NAD+", "NADP+",
                "FAD", "FMN", "Coenzyme A", "Biotin", "Thiamine pyrophosphate",
                "Pyridoxal phosphate", "Folate", "Cobalamin"
            ]
            num_cofactors = random.randint(1, 3)
            properties["cofactors"] = "|".join(random.sample(cofactors, num_cofactors))

        # Generate class names
        if ECAdapterECField.CLASS_NAME in self.fields:
            class_names = {
                1: "Oxidoreductases", 2: "Transferases", 3: "Hydrolases",
                4: "Lyases", 5: "Isomerases", 6: "Ligases", 7: "Translocases"
            }
            properties["class_name"] = class_names.get(ec_class, "Unknown")

        # Generate subclass name
        if ECAdapterECField.SUBCLASS_NAME in self.fields:
            subclass_examples = [
                "Acting on CH-OH group of donors", "Acting on amino acids",
                "Acting on peptide bonds", "Carbon-carbon lyases",
                "Intramolecular transferases", "Forming carbon-oxygen bonds"
            ]
            properties["subclass_name"] = random.choice(subclass_examples)

        # Generate comments
        if ECAdapterECField.COMMENTS in self.fields:
            if random.random() < 0.3:  # 30% chance
                comments = [
                    "Highly conserved across species",
                    "Essential for cellular metabolism",
                    "Requires metal cofactor for activity",
                    "Subject to allosteric regulation",
                    "Found in multiple cellular compartments"
                ]
                properties["comments"] = random.choice(comments)

        # Generate pathway
        if ECAdapterECField.PATHWAY in self.fields:
            pathways = [
                "Glycolysis", "Citric acid cycle", "Pentose phosphate pathway",
                "Fatty acid synthesis", "Amino acid metabolism", "Purine metabolism",
                "Pyrimidine metabolism", "Glycogen metabolism", "Cholesterol biosynthesis"
            ]
            properties["pathway"] = random.choice(pathways)

        return properties