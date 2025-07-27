"""
BioCypher template package adapters.

This package contains various adapters for integrating biological data sources
into BioCypher knowledge graphs. Each adapter follows the BioCypher adapter
pattern with standardized node and edge generation methods.

Available Adapters:
- ExampleAdapter: Basic example adapter for demonstration
- UniProtAdapter: Protein data from UniProt
- CompoundAdapter: Chemical compound and drug data
- DiseaseAdapter: Disease and disorder information
- GOAdapter: Gene Ontology terms and annotations
- PPIAdapter: Protein-protein interactions
- PathwayAdapter: Biological pathways and processes
- PhenotypeAdapter: Phenotypic data and associations
- OrthologyAdapter: Orthologous relationships between proteins

Usage:
    from template_package.adapters import UniProtAdapter
    
    adapter = UniProtAdapter(test_mode=True)
    nodes = list(adapter.get_nodes())
    edges = list(adapter.get_edges())
"""

from .example_adapter import ExampleAdapter
from .uniprot_adapter import UniProtAdapter
from .compound_adapter import CompoundAdapter
from .disease_adapter import DiseaseAdapter
from .go_adapter import GOAdapter
from .ppi_adapter import PPIAdapter
from .pathway_adapter import PathwayAdapter
from .phenotype_adapter import PhenotypeAdapter
from .orthology_adapter import OrthologyAdapter

__all__ = [
    "ExampleAdapter",
    "UniProtAdapter", 
    "CompoundAdapter",
    "DiseaseAdapter",
    "GOAdapter",
    "PPIAdapter",
    "PathwayAdapter",
    "PhenotypeAdapter",
    "OrthologyAdapter",
]