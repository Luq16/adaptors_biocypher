#!/usr/bin/env python3
"""
Test edge generation functionality for BioCypher adaptors.
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

class MockNode:
    """Mock node class for testing edges."""
    def __init__(self, node_id, label, properties=None):
        self.id = node_id
        self.label = label
        self.properties = properties or {}
    
    def get_id(self):
        return self.id
    
    def get_label(self):
        return self.label
    
    def get_properties(self):
        return self.properties

def create_mock_nodes(node_type, count=10):
    """Create mock nodes for testing."""
    nodes = []
    for i in range(count):
        if node_type == "protein":
            node_id = f"P{i:05d}"
        elif node_type == "compound":
            node_id = f"CHEMBL{i+100000}"
        elif node_type == "disease": 
            node_id = f"DOID:{i+10000:05d}"
        else:
            node_id = f"{node_type.upper()}{i:05d}"
        
        nodes.append(MockNode(node_id, node_type))
    
    return nodes

def test_adapter_edges(module_name, adapter_class_name):
    """Test edge generation for a specific adapter."""
    print(f"\nTesting {module_name} edges...")
    
    try:
        # Read and execute the module
        module_path = f"/Users/luqmanawoniyi_1/Documents/biocypher/project-template/template_package/adapters/{module_name}.py"
        
        with open(module_path, 'r') as f:
            content = f.read()
        
        # Replace the biocypher import
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
        
        # Get the adapter class
        adapter_class = exec_globals[adapter_class_name]
        adapter = adapter_class(test_mode=True)
        
        # Generate nodes first
        nodes = list(adapter.get_nodes())
        print(f"  Generated {len(nodes)} nodes")
        
        # Test edge generation based on adapter type
        edges = []
        
        if module_name == "compound_adapter":
            # Test compound-protein and compound-disease edges
            mock_proteins = create_mock_nodes("protein", 5)
            mock_diseases = create_mock_nodes("disease", 3)
            edges = list(adapter.get_edges(target_proteins=mock_proteins, target_diseases=mock_diseases))
            
        elif module_name == "disease_adapter":
            # Test gene-disease edges
            mock_genes = create_mock_nodes("protein", 5)
            edges = list(adapter.get_edges(target_genes=mock_genes))
            
        elif module_name == "go_adapter":
            # Test protein-GO edges
            mock_proteins = create_mock_nodes("protein", 5)
            edges = list(adapter.get_edges(target_proteins=mock_proteins))
            
        elif module_name == "ppi_adapter":
            # Test protein-protein interactions
            edges = list(adapter.get_edges(protein_nodes=nodes))
            
        elif module_name == "pathway_adapter":
            # Test pathway participation edges
            mock_proteins = create_mock_nodes("protein", 5)
            mock_compounds = create_mock_nodes("compound", 3)
            edges = list(adapter.get_edges(target_proteins=mock_proteins, target_compounds=mock_compounds))
            
        elif module_name == "phenotype_adapter":
            # Test phenotype association edges
            mock_genes = create_mock_nodes("protein", 5)
            mock_diseases = create_mock_nodes("disease", 3)
            edges = list(adapter.get_edges(target_genes=mock_genes, target_diseases=mock_diseases))
            
        elif module_name == "orthology_adapter":
            # Test orthology edges
            mock_proteins = create_mock_nodes("protein", 5)
            edges = list(adapter.get_edges(target_proteins=mock_proteins))
            
        elif module_name == "uniprot_adapter":
            # Test protein-protein edges
            edges = list(adapter.get_edges())
        
        print(f"  Generated {len(edges)} edges")
        
        if edges:
            # Show sample edges
            sample_edge = edges[0]
            print(f"  Sample edge: {sample_edge[1]} -> {sample_edge[2]} ({sample_edge[3]})")
            print(f"  Edge properties: {list(sample_edge[4].keys())}")
            
            # Validate edge structure
            assert len(sample_edge) == 5, f"Edge should have 5 elements, got {len(sample_edge)}"
            assert isinstance(sample_edge[4], dict), "Edge properties should be a dictionary"
            
        return True
        
    except Exception as e:
        print(f"  ❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Test edge generation for all adaptors."""
    print("BioCypher Adaptors Edge Generation Test")
    print("=" * 45)
    
    adaptors_to_test = [
        ("uniprot_adapter", "UniProtAdapter"),
        ("compound_adapter", "CompoundAdapter"),
        ("disease_adapter", "DiseaseAdapter"),
        ("go_adapter", "GOAdapter"),
        ("ppi_adapter", "PPIAdapter"),
        ("pathway_adapter", "PathwayAdapter"),
        ("phenotype_adapter", "PhenotypeAdapter"),
        ("orthology_adapter", "OrthologyAdapter"),
    ]
    
    results = {}
    
    for module_name, class_name in adaptors_to_test:
        results[module_name] = test_adapter_edges(module_name, class_name)
    
    # Summary
    print(f"\n{'='*45}")
    print("EDGE GENERATION TEST SUMMARY")
    print(f"{'='*45}")
    
    passed = sum(results.values())
    total = len(results)
    
    for adaptor, success in results.items():
        status = "✅ PASS" if success else "❌ FAIL"
        print(f"{adaptor:20} {status}")
    
    print(f"\nOverall: {passed}/{total} adaptors passed edge tests")
    
    if passed == total:
        print("🎉 All adaptors generate edges correctly!")
        return 0
    else:
        print("⚠️  Some adaptors have edge generation issues.")
        return 1

if __name__ == "__main__":
    sys.exit(main())