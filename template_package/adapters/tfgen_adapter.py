from __future__ import annotations

import random
import string
from enum import Enum, auto
from itertools import chain
from typing import Optional
from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class TFGenAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """
    TRANSCRIPTION_FACTOR = auto()
    REGULATORY_REGION = auto()


class TFGenAdapterTranscriptionFactorField(Enum):
    """
    Define possible fields the adapter can provide for transcription factors.
    """
    NAME = "name"
    GENE_NAME = "gene_name"
    UNIPROT_ID = "uniprot_id"
    ORGANISM = "organism"
    TF_CLASS = "tf_class"
    DNA_BINDING_DOMAIN = "dna_binding_domain"
    COFACTORS = "cofactors"
    CELLULAR_LOCALIZATION = "cellular_localization"


class TFGenAdapterRegulatoryRegionField(Enum):
    """
    Define possible fields for regulatory regions.
    """
    NAME = "name"
    TYPE = "type"
    CHROMOSOME = "chromosome"
    START_POSITION = "start_position"
    END_POSITION = "end_position"
    STRAND = "strand"
    SEQUENCE = "sequence"


class TFGenAdapterEdgeType(Enum):
    """
    Enum for the types of edges the adapter can provide.
    """
    TF_REGULATES_GENE = "tf_regulates_gene"
    TF_BINDS_REGULATORY_REGION = "tf_binds_regulatory_region"
    REGULATORY_REGION_REGULATES_GENE = "regulatory_region_regulates_gene"


class TFGenAdapterTFGeneEdgeField(Enum):
    """
    Define possible fields for TF-gene regulation edges.
    """
    REGULATION_TYPE = "regulation_type"
    EFFECT = "effect"
    CONFIDENCE_SCORE = "confidence_score"
    SOURCE_DATABASE = "source_database"
    PUBMED_IDS = "pubmed_ids"
    EVIDENCE_TYPE = "evidence_type"
    BINDING_SITE_COUNT = "binding_site_count"
    DISTANCE_TO_TSS = "distance_to_tss"
    CONTEXT = "context"
    MECHANISM = "mechanism"


class TFGenAdapter:
    """
    Transcription Factor-Gene Regulation BioCypher adapter. Generates transcription 
    factor nodes and regulatory interaction edges for gene regulation networks.

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
        logger.info("Generating transcription factor and regulatory region nodes.")

        self.nodes = []
        self.transcription_factors = []
        self.regulatory_regions = []

        # Generate transcription factors
        if TFGenAdapterNodeType.TRANSCRIPTION_FACTOR in self.node_types:
            num_tfs = 20 if self.test_mode else 200
            for _ in range(num_tfs):
                tf_node = TranscriptionFactorNode(fields=self.node_fields, organism=self.organism)
                self.transcription_factors.append(tf_node)
                self.nodes.append(tf_node)

        # Generate regulatory regions
        if TFGenAdapterNodeType.REGULATORY_REGION in self.node_types:
            num_regions = 30 if self.test_mode else 300
            for _ in range(num_regions):
                region_node = RegulatoryRegionNode(fields=self.node_fields)
                self.regulatory_regions.append(region_node)
                self.nodes.append(region_node)

        for node in self.nodes:
            yield (node.get_id(), node.get_label(), node.get_properties())

    def get_edges(self, target_genes=None):
        """
        Returns a generator of edge tuples for edge types specified in the
        adapter constructor.

        Args:
            target_genes: List of gene/protein nodes to create regulation relationships with.
        """
        logger.info("Generating transcription factor regulation edges.")

        if not self.nodes:
            list(self.get_nodes())  # Generate nodes first

        # TF-Gene regulation relationships
        if (TFGenAdapterEdgeType.TF_REGULATES_GENE in self.edge_types 
            and target_genes):
            
            for tf in self.transcription_factors:
                # Each TF regulates 2-8 genes
                if random.random() < 0.8:  # 80% of TFs have regulatory targets
                    num_targets = random.randint(2, 8)
                    selected_genes = random.sample(target_genes, min(num_targets, len(target_genes)))
                    
                    for gene in selected_genes:
                        gene_id = gene[0] if isinstance(gene, tuple) else str(gene)
                        
                        relationship_id = "".join(
                            random.choice(string.ascii_letters + string.digits)
                            for _ in range(12)
                        )
                        
                        properties = {}
                        for field in self.edge_fields:
                            if field == TFGenAdapterTFGeneEdgeField.REGULATION_TYPE:
                                properties["regulation_type"] = random.choice([
                                    "transcriptional regulation", "post-transcriptional regulation",
                                    "chromatin remodeling", "epigenetic regulation"
                                ])
                            elif field == TFGenAdapterTFGeneEdgeField.EFFECT:
                                properties["effect"] = random.choice([
                                    "activation", "repression", "dual", "context-dependent"
                                ])
                            elif field == TFGenAdapterTFGeneEdgeField.CONFIDENCE_SCORE:
                                properties["confidence_score"] = round(random.uniform(0.3, 1.0), 3)
                            elif field == TFGenAdapterTFGeneEdgeField.SOURCE_DATABASE:
                                properties["source_database"] = random.choice([
                                    "DoRothEA", "CollecTRI", "TRRUST", "RegNetwork",
                                    "JASPAR", "TRANSFAC", "ENCODE", "ChIP-Atlas"
                                ])
                            elif field == TFGenAdapterTFGeneEdgeField.PUBMED_IDS:
                                num_pubmed = random.randint(1, 4)
                                pubmed_ids = [str(random.randint(15000000, 35000000)) for _ in range(num_pubmed)]
                                properties["pubmed_ids"] = "|".join(pubmed_ids)
                            elif field == TFGenAdapterTFGeneEdgeField.EVIDENCE_TYPE:
                                properties["evidence_type"] = random.choice([
                                    "ChIP-seq", "EMSA", "reporter assay", "qPCR",
                                    "RNA-seq", "microarray", "luciferase assay",
                                    "DNase-seq", "ATAC-seq", "ChIP-chip"
                                ])
                            elif field == TFGenAdapterTFGeneEdgeField.BINDING_SITE_COUNT:
                                properties["binding_site_count"] = random.randint(1, 5)
                            elif field == TFGenAdapterTFGeneEdgeField.DISTANCE_TO_TSS:
                                # Distance to transcription start site
                                properties["distance_to_tss"] = random.randint(-50000, 10000)
                            elif field == TFGenAdapterTFGeneEdgeField.CONTEXT:
                                properties["context"] = random.choice([
                                    "cell cycle", "stress response", "development",
                                    "metabolism", "differentiation", "apoptosis",
                                    "immune response", "cancer"
                                ])
                            elif field == TFGenAdapterTFGeneEdgeField.MECHANISM:
                                properties["mechanism"] = random.choice([
                                    "direct binding", "indirect regulation", "chromatin looping",
                                    "enhancer-promoter interaction", "cooperative binding",
                                    "competitive binding", "allosteric regulation"
                                ])

                        yield (
                            relationship_id,
                            tf.get_id(),
                            gene_id,
                            TFGenAdapterEdgeType.TF_REGULATES_GENE.value,
                            properties,
                        )

        # TF-Regulatory Region binding relationships
        if (TFGenAdapterEdgeType.TF_BINDS_REGULATORY_REGION in self.edge_types
            and self.regulatory_regions):
            
            for tf in self.transcription_factors:
                # Each TF binds to 1-4 regulatory regions
                if random.random() < 0.6:  # 60% of TFs have binding sites
                    num_regions = random.randint(1, 4)
                    selected_regions = random.sample(self.regulatory_regions, min(num_regions, len(self.regulatory_regions)))
                    
                    for region in selected_regions:
                        relationship_id = "".join(
                            random.choice(string.ascii_letters + string.digits)
                            for _ in range(12)
                        )
                        
                        properties = {
                            "binding_affinity": f"{random.uniform(0.1, 100):.2f} nM",
                            "binding_strength": random.choice(["weak", "moderate", "strong"]),
                            "cooperative_binding": random.choice(["true", "false"]),
                            "evidence": random.choice(["ChIP-seq", "EMSA", "DNase footprinting"])
                        }

                        yield (
                            relationship_id,
                            tf.get_id(),
                            region.get_id(),
                            TFGenAdapterEdgeType.TF_BINDS_REGULATORY_REGION.value,
                            properties,
                        )

        # Regulatory Region-Gene relationships
        if (TFGenAdapterEdgeType.REGULATORY_REGION_REGULATES_GENE in self.edge_types
            and target_genes and self.regulatory_regions):
            
            for region in self.regulatory_regions:
                # Each regulatory region affects 1-3 genes
                if random.random() < 0.7:  # 70% of regions regulate genes
                    num_genes = random.randint(1, 3)
                    selected_genes = random.sample(target_genes, min(num_genes, len(target_genes)))
                    
                    for gene in selected_genes:
                        gene_id = gene[0] if isinstance(gene, tuple) else str(gene)
                        
                        relationship_id = "".join(
                            random.choice(string.ascii_letters + string.digits)
                            for _ in range(12)
                        )
                        
                        properties = {
                            "regulatory_effect": random.choice([
                                "enhancer", "silencer", "promoter", "insulator"
                            ]),
                            "distance": random.randint(100, 100000),
                            "interaction_type": random.choice([
                                "direct", "looping", "chromatin_mediated"
                            ])
                        }

                        yield (
                            relationship_id,
                            region.get_id(),
                            gene_id,
                            TFGenAdapterEdgeType.REGULATORY_REGION_REGULATES_GENE.value,
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
            self.node_types = [type for type in TFGenAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [field for field in TFGenAdapterTranscriptionFactorField]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in TFGenAdapterEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in TFGenAdapterTFGeneEdgeField]


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


class TranscriptionFactorNode(Node):
    """
    Generates instances of transcription factor nodes.
    """

    def __init__(self, fields: Optional[list] = None, organism: str = "9606"):
        self.fields = fields
        self.organism = organism
        self.id = self._generate_id()
        self.label = "transcription_factor"
        self.properties = self._generate_properties()

    def _generate_id(self):
        """
        Generate a transcription factor identifier.
        """
        # Use UniProt-style ID for TFs
        prefix = "".join(random.choices(string.ascii_uppercase, k=random.choice([1, 2])))
        suffix = "".join(random.choices(string.ascii_uppercase + string.digits, k=random.choice([4, 5])))
        return f"{prefix}{suffix}"

    def _generate_properties(self):
        properties = {}

        # Add organism information
        if TFGenAdapterTranscriptionFactorField.ORGANISM in self.fields:
            properties["organism"] = self.organism

        # Generate TF name
        if TFGenAdapterTranscriptionFactorField.NAME in self.fields:
            tf_names = [
                "Nuclear factor kappa B subunit 1", "Transcription factor AP-1",
                "CCAAT/enhancer-binding protein alpha", "Signal transducer and activator of transcription 3",
                "Nuclear factor of activated T-cells 1", "E2F transcription factor 1",
                "Forkhead box protein O1", "Hypoxia-inducible factor 1-alpha",
                "Serum response factor", "Specificity protein 1", "c-Myc proto-oncogene protein",
                "p53 tumor suppressor", "Retinoid X receptor alpha", "Peroxisome proliferator-activated receptor gamma",
                "GATA-binding factor 1", "Runt-related transcription factor 1",
                "T-box transcription factor TBX21", "Homeobox protein NANOG",
                "SRY-box transcription factor SOX2", "POU domain class 5 transcription factor 1"
            ]
            properties["name"] = random.choice(tf_names)

        # Generate gene name
        if TFGenAdapterTranscriptionFactorField.GENE_NAME in self.fields:
            gene_names = [
                "NFKB1", "JUN", "CEBPA", "STAT3", "NFATC1", "E2F1", "FOXO1", "HIF1A",
                "SRF", "SP1", "MYC", "TP53", "RXRA", "PPARG", "GATA1", "RUNX1",
                "TBX21", "NANOG", "SOX2", "POU5F1", "CTCF", "YY1", "ELK1", "ATF2"
            ]
            properties["gene_name"] = random.choice(gene_names)

        # Generate UniProt ID
        if TFGenAdapterTranscriptionFactorField.UNIPROT_ID in self.fields:
            properties["uniprot_id"] = self.id

        # Generate TF class
        if TFGenAdapterTranscriptionFactorField.TF_CLASS in self.fields:
            tf_classes = [
                "Basic helix-loop-helix", "Zinc finger", "Homeodomain", "Basic leucine zipper",
                "Helix-turn-helix", "High mobility group", "Nuclear receptor", "Immunoglobulin fold",
                "beta-Scaffold factors with minor groove contacts", "Other transcription factors"
            ]
            properties["tf_class"] = random.choice(tf_classes)

        # Generate DNA binding domain
        if TFGenAdapterTranscriptionFactorField.DNA_BINDING_DOMAIN in self.fields:
            dbd_types = [
                "C2H2 zinc finger", "Basic helix-loop-helix", "Homeodomain",
                "Basic leucine zipper", "Nuclear receptor", "High mobility group box",
                "Helix-turn-helix", "Immunoglobulin fold", "beta-Sheet", "Other"
            ]
            properties["dna_binding_domain"] = random.choice(dbd_types)

        # Generate cofactors
        if TFGenAdapterTranscriptionFactorField.COFACTORS in self.fields:
            cofactors = [
                "p300", "CBP", "Mediator complex", "TFIID", "RNA polymerase II",
                "Chromatin remodeling complexes", "Histone modifying enzymes"
            ]
            if random.random() < 0.6:  # 60% chance
                properties["cofactors"] = random.choice(cofactors)

        # Generate cellular localization
        if TFGenAdapterTranscriptionFactorField.CELLULAR_LOCALIZATION in self.fields:
            properties["cellular_localization"] = random.choice([
                "nucleus", "cytoplasm", "nucleus and cytoplasm", "membrane-associated"
            ])

        return properties


class RegulatoryRegionNode(Node):
    """
    Generates instances of regulatory region nodes.
    """

    def __init__(self, fields: Optional[list] = None):
        self.fields = fields
        self.id = self._generate_id()
        self.label = "regulatory_region"
        self.properties = self._generate_properties()

    def _generate_id(self):
        """
        Generate a regulatory region identifier.
        """
        # Format: chr:start-end
        chromosome = random.randint(1, 22) if random.random() < 0.95 else random.choice(["X", "Y"])
        start = random.randint(1000000, 200000000)
        end = start + random.randint(200, 2000)
        return f"chr{chromosome}:{start}-{end}"

    def _generate_properties(self):
        properties = {}

        # Parse coordinates from ID
        parts = self.id.split(":")
        chromosome = parts[0]
        coords = parts[1].split("-")
        start_pos = int(coords[0])
        end_pos = int(coords[1])

        # Generate name
        if TFGenAdapterRegulatoryRegionField.NAME in self.fields:
            region_types = ["enhancer", "promoter", "silencer", "insulator"]
            region_type = random.choice(region_types)
            properties["name"] = f"{region_type}_{chromosome}_{start_pos}"

        # Generate type
        if TFGenAdapterRegulatoryRegionField.TYPE in self.fields:
            properties["type"] = random.choice([
                "enhancer", "promoter", "silencer", "insulator", "operator", "TFIID_binding_site"
            ])

        # Generate chromosome
        if TFGenAdapterRegulatoryRegionField.CHROMOSOME in self.fields:
            properties["chromosome"] = chromosome

        # Generate positions
        if TFGenAdapterRegulatoryRegionField.START_POSITION in self.fields:
            properties["start_position"] = start_pos
        
        if TFGenAdapterRegulatoryRegionField.END_POSITION in self.fields:
            properties["end_position"] = end_pos

        # Generate strand
        if TFGenAdapterRegulatoryRegionField.STRAND in self.fields:
            properties["strand"] = random.choice(["+", "-"])

        # Generate sequence
        if TFGenAdapterRegulatoryRegionField.SEQUENCE in self.fields:
            sequence_length = end_pos - start_pos
            nucleotides = ["A", "T", "G", "C"]
            properties["sequence"] = "".join(random.choices(nucleotides, k=min(sequence_length, 50)))

        return properties