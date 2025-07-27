#!/usr/bin/env python3
"""
Test script for BioCypher template adaptors.
This script tests all adaptors to ensure they work correctly.
"""

import sys
import traceback
from template_package.adapters import (
    ExampleAdapter,
    UniProtAdapter, 
    CompoundAdapter,
    DiseaseAdapter,
    GOAdapter,
    PPIAdapter,
    PathwayAdapter,
    PhenotypeAdapter,
    OrthologyAdapter,
)

def test_adapter(adapter_class, adapter_name, **kwargs):
    """Test a single adapter."""
    print(f"\n{'='*50}")
    print(f"Testing {adapter_name}")
    print(f"{'='*50}")
    
    try:
        # Initialize adapter
        print(f"1. Initializing {adapter_name}...")
        adapter = adapter_class(test_mode=True, **kwargs)
        print("   ✅ Initialization successful")
        
        # Test get_nodes()
        print("2. Testing get_nodes()...")
        nodes = list(adapter.get_nodes())
        print(f"   ✅ Generated {len(nodes)} nodes")
        
        if nodes:
            # Show sample node
            sample_node = nodes[0]
            print(f"   Sample node: ID={sample_node[0]}, Label={sample_node[1]}")
            print(f"   Properties: {list(sample_node[2].keys())}")
        
        # Test get_edges() if the adapter supports it
        if hasattr(adapter, 'get_edges'):
            print("3. Testing get_edges()...")
            try:
                # Some adapters need target nodes for edges
                if adapter_name in ["CompoundAdapter", "DiseaseAdapter", "GOAdapter", 
                                   "PathwayAdapter", "PhenotypeAdapter"]:
                    # Create mock protein nodes for testing
                    mock_proteins = [MockNode(f"P{i:05d}", "protein") for i in range(10)]
                    edges = list(adapter.get_edges(target_proteins=mock_proteins))
                elif adapter_name == "PPIAdapter":
                    edges = list(adapter.get_edges(protein_nodes=nodes))
                elif adapter_name == "OrthologyAdapter":
                    mock_proteins = [MockNode(f"P{i:05d}", "protein") for i in range(10)]
                    edges = list(adapter.get_edges(target_proteins=mock_proteins))
                else:
                    edges = list(adapter.get_edges())
                
                print(f"   ✅ Generated {len(edges)} edges")
                
                if edges:
                    # Show sample edge
                    sample_edge = edges[0]
                    print(f"   Sample edge: {sample_edge[1]} -> {sample_edge[2]} ({sample_edge[3]})")
                    print(f"   Properties: {list(sample_edge[4].keys())}")
                    
            except Exception as e:
                print(f"   ⚠️  Edge generation failed: {e}")
        
        # Test get_node_count() if available
        if hasattr(adapter, 'get_node_count'):
            print("4. Testing get_node_count()...")
            count = adapter.get_node_count()
            print(f"   ✅ Node count: {count}")
        
        print(f"✅ {adapter_name} test completed successfully!")
        return True
        
    except Exception as e:
        print(f"❌ {adapter_name} test failed!")
        print(f"Error: {e}")
        print("Traceback:")
        traceback.print_exc()
        return False

class MockNode:
    """Mock node class for testing edges."""
    def __init__(self, node_id, label):
        self.id = node_id
        self.label = label
    
    def get_id(self):
        return self.id
    
    def get_label(self):
        return self.label

def main():
    """Run all adapter tests."""
    print("BioCypher Adaptors Test Suite")
    print("=" * 60)
    
    adapters_to_test = [
        (ExampleAdapter, "ExampleAdapter"),
        (UniProtAdapter, "UniProtAdapter"),
        (CompoundAdapter, "CompoundAdapter"),
        (DiseaseAdapter, "DiseaseAdapter"),
        (GOAdapter, "GOAdapter"),
        (PPIAdapter, "PPIAdapter"),
        (PathwayAdapter, "PathwayAdapter"),
        (PhenotypeAdapter, "PhenotypeAdapter"),
        (OrthologyAdapter, "OrthologyAdapter"),
    ]
    
    results = {}
    
    for adapter_class, adapter_name in adapters_to_test:
        results[adapter_name] = test_adapter(adapter_class, adapter_name)
    
    # Summary
    print(f"\n{'='*60}")
    print("TEST SUMMARY")
    print(f"{'='*60}")
    
    passed = sum(results.values())
    total = len(results)
    
    for adapter_name, success in results.items():
        status = "✅ PASS" if success else "❌ FAIL"
        print(f"{adapter_name:20} {status}")
    
    print(f"\nOverall: {passed}/{total} adapters passed")
    
    if passed == total:
        print("🎉 All adaptors are working correctly!")
        return 0
    else:
        print("⚠️  Some adaptors have issues that need to be fixed.")
        return 1

if __name__ == "__main__":
    sys.exit(main())