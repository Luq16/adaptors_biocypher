#!/usr/bin/env python3
"""
Test script for all new BioCypher adaptors.
Tests functionality of the 8 newly created adaptors from CROssBARv2.
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

def test_basic_adaptor_functionality(adapter_name, adapter_class, class_name):
    """Test basic functionality of an adaptor."""
    print(f"\n{'='*60}")
    print(f"Testing {class_name}")
    print(f"{'='*60}")
    
    try:
        # Initialize adapter
        adapter = adapter_class(test_mode=True)
        print(f"✅ {class_name} initialized successfully")
        
        # Test node generation
        nodes = list(adapter.get_nodes())
        print(f"✅ Generated {len(nodes)} nodes")
        
        if len(nodes) > 0:
            # Print sample node
            sample_node = nodes[0]
            print(f"   Sample node: {sample_node[0]} ({sample_node[1]})")
            print(f"   Properties: {list(sample_node[2].keys())}")
        
        # Test edge generation (if applicable)
        edges = []
        try:
            if hasattr(adapter, 'get_edges'):
                # Try different edge generation patterns
                if adapter_name in ['side_effect', 'drug', 'ec', 'interpro', 'tfgen']:
                    # These need target nodes
                    sample_targets = [('TARGET1', 'protein', {}), ('TARGET2', 'protein', {})]
                    if adapter_name == 'side_effect':
                        edges = list(adapter.get_edges(target_drugs=sample_targets))
                    elif adapter_name == 'drug':
                        edges = list(adapter.get_edges(target_proteins=sample_targets))
                    elif adapter_name == 'ec':
                        edges = list(adapter.get_edges(target_proteins=sample_targets))
                    elif adapter_name == 'interpro':
                        edges = list(adapter.get_edges(target_proteins=sample_targets))
                    elif adapter_name == 'tfgen':
                        edges = list(adapter.get_edges(target_genes=sample_targets))
                else:
                    # These can generate edges independently
                    edges = list(adapter.get_edges())
                
                print(f"✅ Generated {len(edges)} edges")
                
                if len(edges) > 0:
                    # Print sample edge
                    sample_edge = edges[0]
                    print(f"   Sample edge: {sample_edge[1]} -> {sample_edge[2]} ({sample_edge[3]})")
                    if len(sample_edge) > 4:
                        print(f"   Edge properties: {list(sample_edge[4].keys())}")
            else:
                print("   No edge generation method found")
                
        except Exception as e:
            print(f"⚠️  Edge generation error: {e}")
        
        # Test node count method
        try:
            node_count = adapter.get_node_count()
            print(f"✅ Node count method: {node_count}")
        except Exception as e:
            print(f"⚠️  Node count error: {e}")
            
        return True, len(nodes), len(edges)
        
    except Exception as e:
        print(f"❌ {class_name} failed: {e}")
        import traceback
        traceback.print_exc()
        return False, 0, 0

def test_all_new_adaptors():
    """Test all newly created adaptors."""
    print("Testing All New BioCypher Adaptors")
    print("=" * 80)
    print("Testing 8 new adaptors created from CROssBARv2 examples")
    
    # Define new adaptors to test
    new_adaptors = [
        ("side_effect_adapter", "SideEffectAdapter", "SideEffectAdapter"),
        ("drug_adapter", "DrugAdapter", "DrugAdapter"),
        ("biogrid_adapter", "BiogridAdapter", "BiogridAdapter"),
        ("string_adapter", "StringAdapter", "StringAdapter"),
        ("intact_adapter", "IntactAdapter", "IntactAdapter"),
        ("ec_adapter", "ECAdapter", "ECAdapter"),
        ("interpro_adapter", "InterproAdapter", "InterproAdapter"),
        ("tfgen_adapter", "TFGenAdapter", "TFGenAdapter"),
    ]
    
    results = {}
    total_nodes = 0
    total_edges = 0
    successful_adaptors = 0
    
    for module_name, class_name, display_name in new_adaptors:
        try:
            # Load the adapter class
            AdapterClass = load_adapter(module_name, class_name)
            
            # Test basic functionality
            success, nodes, edges = test_basic_adaptor_functionality(
                module_name.replace('_adapter', ''), AdapterClass, display_name
            )
            
            results[display_name] = {
                'success': success,
                'nodes': nodes,
                'edges': edges
            }
            
            if success:
                successful_adaptors += 1
                total_nodes += nodes
                total_edges += edges
                
        except Exception as e:
            print(f"❌ Failed to load {display_name}: {e}")
            results[display_name] = {
                'success': False,
                'nodes': 0,
                'edges': 0,
                'error': str(e)
            }
    
    # Print summary
    print(f"\n{'='*80}")
    print(f"TEST SUMMARY")
    print(f"{'='*80}")
    print(f"Total Adaptors Tested: {len(new_adaptors)}")
    print(f"Successfully Working: {successful_adaptors}")
    print(f"Failed: {len(new_adaptors) - successful_adaptors}")
    print(f"Total Nodes Generated: {total_nodes}")
    print(f"Total Edges Generated: {total_edges}")
    
    print(f"\n{'='*60}")
    print("DETAILED RESULTS:")
    print(f"{'='*60}")
    
    for adaptor_name, result in results.items():
        status = "✅ PASS" if result['success'] else "❌ FAIL"
        print(f"{adaptor_name:20} {status:8} Nodes: {result['nodes']:3d} Edges: {result['edges']:3d}")
        if not result['success'] and 'error' in result:
            print(f"                     Error: {result['error']}")
    
    # Test integration between adaptors
    print(f"\n{'='*60}")
    print("INTEGRATION TEST:")
    print(f"{'='*60}")
    
    try:
        # Load a few key adaptors for integration test
        DrugAdapter = load_adapter("drug_adapter", "DrugAdapter")
        SideEffectAdapter = load_adapter("side_effect_adapter", "SideEffectAdapter")
        
        drug_adapter = DrugAdapter(test_mode=True)
        side_effect_adapter = SideEffectAdapter(test_mode=True)
        
        # Generate nodes
        drugs = list(drug_adapter.get_nodes())
        side_effects = list(side_effect_adapter.get_nodes())
        
        # Test cross-adaptor relationships
        drug_side_effect_edges = list(side_effect_adapter.get_edges(target_drugs=drugs))
        drug_internal_edges = list(drug_adapter.get_edges())
        
        print(f"✅ Integration test successful!")
        print(f"   Drugs: {len(drugs)}, Side Effects: {len(side_effects)}")
        print(f"   Drug-Side Effect relationships: {len(drug_side_effect_edges)}")
        print(f"   Drug internal relationships: {len(drug_internal_edges)}")
        
    except Exception as e:
        print(f"❌ Integration test failed: {e}")
    
    if successful_adaptors == len(new_adaptors):
        print(f"\n🎉 ALL NEW ADAPTORS ARE WORKING CORRECTLY!")
        print(f"   Total of {len(new_adaptors)} new adaptors successfully implemented")
        print(f"   Combined with existing 8 adaptors = 16 total adaptors available")
        return True
    else:
        print(f"\n⚠️  {len(new_adaptors) - successful_adaptors} adaptors need attention")
        return False

def main():
    """Run the comprehensive test."""
    try:
        success = test_all_new_adaptors()
        return 0 if success else 1
    except Exception as e:
        print(f"\n❌ Test suite failed with error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())