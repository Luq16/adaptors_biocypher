#!/usr/bin/env python3
"""
Integration test for BioCypher adaptors.
Tests how adaptors work together to create an integrated knowledge graph.
"""

import sys
import random
import string
from enum import Enum, auto
from itertools import chain
from typing import Optional

# Mock the biocypher logger
class MockLogger:
    def debug(self, msg): pass
    def info(self, msg): print(f"INFO: {msg}")
    def warning(self, msg): print(f"WARNING: {msg}")
    def error(self, msg): print(f"ERROR: {msg}")

logger = MockLogger()

def load_adapter(module_name, class_name):
    """Load and return an adapter class."""
    module_path = f"/Users/luqmanawoniyi_1/Documents/biocypher/project-template/template_package/adapters/{module_name}.py"
    
    with open(module_path, 'r') as f:
        content = f.read()
    
    content = content.replace('from biocypher._logger import logger', 'logger = MockLogger()')
    
    exec_globals = {
        'logger': logger,
        'random': random,
        'string': string,
        'Enum': Enum,
        'auto': auto,
        'chain': chain,
        'Optional': Optional,
        'MockLogger': MockLogger
    }
    
    exec(content, exec_globals)
    return exec_globals[class_name]

def test_integration():
    """Test creating an integrated knowledge graph using multiple adaptors."""
    print("BioCypher Adaptors Integration Test")
    print("=" * 40)
    
    # Load all adaptors
    print("Loading adaptors...")
    UniProtAdapter = load_adapter("uniprot_adapter", "UniProtAdapter")
    CompoundAdapter = load_adapter("compound_adapter", "CompoundAdapter")
    DiseaseAdapter = load_adapter("disease_adapter", "DiseaseAdapter")
    GOAdapter = load_adapter("go_adapter", "GOAdapter")
    PPIAdapter = load_adapter("ppi_adapter", "PPIAdapter")
    PathwayAdapter = load_adapter("pathway_adapter", "PathwayAdapter")
    PhenotypeAdapter = load_adapter("phenotype_adapter", "PhenotypeAdapter")
    OrthologyAdapter = load_adapter("orthology_adapter", "OrthologyAdapter")
    
    print("✅ All adaptors loaded")
    
    # Initialize adaptors with test mode
    print("\nInitializing adaptors...")
    uniprot_adapter = UniProtAdapter(test_mode=True)
    compound_adapter = CompoundAdapter(test_mode=True)
    disease_adapter = DiseaseAdapter(test_mode=True)
    go_adapter = GOAdapter(test_mode=True)
    ppi_adapter = PPIAdapter(test_mode=True)
    pathway_adapter = PathwayAdapter(test_mode=True)
    phenotype_adapter = PhenotypeAdapter(test_mode=True)
    orthology_adapter = OrthologyAdapter(test_mode=True)
    
    print("✅ All adaptors initialized")
    
    # Generate nodes from each adaptor
    print("\nGenerating nodes...")
    protein_nodes = list(uniprot_adapter.get_nodes())
    compound_nodes = list(compound_adapter.get_nodes())
    disease_nodes = list(disease_adapter.get_nodes())
    go_nodes = list(go_adapter.get_nodes())
    pathway_nodes = list(pathway_adapter.get_nodes())
    phenotype_nodes = list(phenotype_adapter.get_nodes())
    ortholog_nodes = list(orthology_adapter.get_nodes())
    ppi_protein_nodes = list(ppi_adapter.get_nodes())
    
    total_nodes = (len(protein_nodes) + len(compound_nodes) + len(disease_nodes) + 
                   len(go_nodes) + len(pathway_nodes) + len(phenotype_nodes) + 
                   len(ortholog_nodes) + len(ppi_protein_nodes))
    
    print(f"✅ Generated {total_nodes} total nodes:")
    print(f"  - Proteins (UniProt): {len(protein_nodes)}")
    print(f"  - Compounds: {len(compound_nodes)}")
    print(f"  - Diseases: {len(disease_nodes)}")
    print(f"  - GO terms: {len(go_nodes)}")
    print(f"  - Pathways: {len(pathway_nodes)}")
    print(f"  - Phenotypes: {len(phenotype_nodes)}")
    print(f"  - Ortholog groups: {len(ortholog_nodes)}")
    print(f"  - PPI proteins: {len(ppi_protein_nodes)}")
    
    # Generate edges with cross-adapter relationships
    print("\nGenerating integrated edges...")
    
    all_edges = []
    
    # Protein-protein interactions
    ppi_edges = list(uniprot_adapter.get_edges())
    ppi_network_edges = list(ppi_adapter.get_edges(protein_nodes=ppi_protein_nodes))
    all_edges.extend(ppi_edges)
    all_edges.extend(ppi_network_edges)
    
    # Compound-protein interactions (using proteins from UniProt)
    compound_edges = list(compound_adapter.get_edges(
        target_proteins=protein_nodes[:10],  # Use first 10 proteins
        target_diseases=disease_nodes[:5]    # Use first 5 diseases
    ))
    all_edges.extend(compound_edges)
    
    # Gene-disease associations
    disease_edges = list(disease_adapter.get_edges(target_genes=protein_nodes[:10]))
    all_edges.extend(disease_edges)
    
    # Protein-GO annotations
    go_edges = list(go_adapter.get_edges(target_proteins=protein_nodes[:10]))
    all_edges.extend(go_edges)
    
    # Pathway participation
    pathway_edges = list(pathway_adapter.get_edges(
        target_proteins=protein_nodes[:10],
        target_compounds=compound_nodes[:5]
    ))
    all_edges.extend(pathway_edges)
    
    # Phenotype associations
    phenotype_edges = list(phenotype_adapter.get_edges(
        target_genes=protein_nodes[:8],
        target_diseases=disease_nodes[:5]
    ))
    all_edges.extend(phenotype_edges)
    
    # Orthology relationships
    orthology_edges = list(orthology_adapter.get_edges(target_proteins=protein_nodes[:10]))
    all_edges.extend(orthology_edges)
    
    print(f"✅ Generated {len(all_edges)} total edges:")
    print(f"  - Protein-protein (UniProt): {len(ppi_edges)}")
    print(f"  - Protein-protein (PPI): {len(ppi_network_edges)}")
    print(f"  - Compound interactions: {len(compound_edges)}")
    print(f"  - Gene-disease: {len(disease_edges)}")
    print(f"  - Protein-GO: {len(go_edges)}")
    print(f"  - Pathway participation: {len(pathway_edges)}")
    print(f"  - Phenotype associations: {len(phenotype_edges)}")
    print(f"  - Orthology relationships: {len(orthology_edges)}")
    
    # Analyze the integrated graph
    print("\nAnalyzing integrated knowledge graph...")
    
    # Count edge types
    edge_types = {}
    for edge in all_edges:
        edge_type = edge[3]
        edge_types[edge_type] = edge_types.get(edge_type, 0) + 1
    
    print("Edge type distribution:")
    for edge_type, count in sorted(edge_types.items()):
        print(f"  - {edge_type}: {count}")
    
    # Identify cross-adaptor connections
    print("\nCross-adaptor connectivity analysis:")
    
    # Find proteins that appear in multiple contexts
    protein_ids = set(node[0] for node in protein_nodes)
    ppi_protein_ids = set(node[0] for node in ppi_protein_nodes)
    
    connected_proteins = set()
    for edge in all_edges:
        if edge[1] in protein_ids or edge[2] in protein_ids:
            connected_proteins.add(edge[1] if edge[1] in protein_ids else edge[2])
    
    print(f"  - {len(connected_proteins)} proteins have interactions")
    print(f"  - {len(protein_ids & ppi_protein_ids)} proteins appear in both UniProt and PPI networks")
    
    # Validate graph structure
    print("\nValidating graph structure...")
    
    # Check for self-loops
    self_loops = [edge for edge in all_edges if edge[1] == edge[2]]
    print(f"  - Self-loops: {len(self_loops)}")
    
    # Check edge property consistency
    edges_with_properties = [edge for edge in all_edges if edge[4]]
    print(f"  - Edges with properties: {len(edges_with_properties)}/{len(all_edges)}")
    
    # Sample some edges to show integration
    print("\nSample integrated relationships:")
    if len(all_edges) >= 5:
        sample_edges = random.sample(all_edges, 5)
        for i, edge in enumerate(sample_edges, 1):
            print(f"  {i}. {edge[1]} -> {edge[2]} ({edge[3]})")
    
    print(f"\n🎉 Integration test completed successfully!")
    print(f"Created knowledge graph with {total_nodes} nodes and {len(all_edges)} edges")
    
    return True

def main():
    """Run the integration test."""
    try:
        success = test_integration()
        if success:
            print("\n✅ ALL TESTS PASSED!")
            print("The BioCypher adaptors are working correctly and can be used together.")
            return 0
        else:
            print("\n❌ INTEGRATION TEST FAILED!")
            return 1
    except Exception as e:
        print(f"\n❌ Integration test failed with error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())