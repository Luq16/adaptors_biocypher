#!/usr/bin/env python3
"""
Simple test for BioCypher template adaptors without full BioCypher dependency.
This tests the basic structure and logic of our adaptors.
"""

import sys
import os
import random
import string
from enum import Enum, auto
from itertools import chain
from typing import Optional

# Mock the biocypher logger to avoid import issues
class MockLogger:
    def debug(self, msg): print(f"DEBUG: {msg}")
    def info(self, msg): print(f"INFO: {msg}")
    def warning(self, msg): print(f"WARNING: {msg}")
    def error(self, msg): print(f"ERROR: {msg}")

logger = MockLogger()

# Add the template_package to path
sys.path.insert(0, '/Users/luqmanawoniyi_1/Documents/biocypher/project-template')

def test_basic_adaptor_structure():
    """Test that we can import and instantiate our adaptors."""
    print("Testing Basic Adaptor Structure")
    print("=" * 40)
    
    # Test importing enum classes from each adaptor
    adaptor_modules = [
        'uniprot_adapter',
        'compound_adapter', 
        'disease_adapter',
        'go_adapter',
        'ppi_adapter',
        'pathway_adapter',
        'phenotype_adapter',
        'orthology_adapter'
    ]
    
    results = {}
    
    for module_name in adaptor_modules:
        try:
            print(f"\nTesting {module_name}...")
            
            # Import the module file directly
            module_path = f"/Users/luqmanawoniyi_1/Documents/biocypher/project-template/template_package/adapters/{module_name}.py"
            
            with open(module_path, 'r') as f:
                content = f.read()
            
            # Replace the biocypher import with our mock
            content = content.replace('from biocypher._logger import logger', 'logger = MockLogger()')
            
            # Create a temporary module
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
            
            # Find the main adapter class
            adapter_class = None
            for name, obj in exec_globals.items():
                if name.endswith('Adapter') and name != 'ExampleAdapter':
                    if isinstance(obj, type) and hasattr(obj, 'get_nodes'):
                        adapter_class = obj
                        break
            
            if adapter_class:
                print(f"  ✅ Found adapter class: {adapter_class.__name__}")
                
                # Test instantiation
                adapter = adapter_class(test_mode=True)
                print(f"  ✅ Successfully instantiated {adapter_class.__name__}")
                
                # Test get_nodes
                nodes = list(adapter.get_nodes())
                print(f"  ✅ Generated {len(nodes)} nodes")
                
                if nodes:
                    sample_node = nodes[0]
                    print(f"  Sample: ID={sample_node[0]}, Label={sample_node[1]}, Props={len(sample_node[2])}")
                
                results[module_name] = True
            else:
                print(f"  ❌ No adapter class found in {module_name}")
                results[module_name] = False
                
        except Exception as e:
            print(f"  ❌ Error testing {module_name}: {e}")
            results[module_name] = False
    
    return results

def test_example_adapter():
    """Test the original example adapter works."""
    print("\nTesting Example Adapter")
    print("=" * 30)
    
    try:
        # Read and execute the example adapter
        with open('/Users/luqmanawoniyi_1/Documents/biocypher/project-template/template_package/adapters/example_adapter.py', 'r') as f:
            content = f.read()
        
        # Replace the biocypher import with our mock
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
        
        # Test the ExampleAdapter
        ExampleAdapter = exec_globals['ExampleAdapter']
        adapter = ExampleAdapter()
        
        print("✅ ExampleAdapter instantiated")
        
        # Test nodes
        nodes = list(adapter.get_nodes())
        print(f"✅ Generated {len(nodes)} nodes")
        
        if nodes:
            # Test edges
            edges = list(adapter.get_edges())
            print(f"✅ Generated {len(edges)} edges")
            
            if edges:
                sample_edge = edges[0]
                print(f"Sample edge: {sample_edge[1]} -> {sample_edge[2]} ({sample_edge[3]})")
        
        return True
        
    except Exception as e:
        print(f"❌ ExampleAdapter test failed: {e}")
        return False

def main():
    """Run all tests."""
    print("BioCypher Adaptors Simple Test Suite")
    print("=" * 50)
    
    # Test example adapter first
    example_result = test_example_adapter()
    
    # Test our new adaptors
    adaptor_results = test_basic_adaptor_structure()
    
    # Summary
    print(f"\n{'='*50}")
    print("TEST SUMMARY")
    print(f"{'='*50}")
    
    print(f"ExampleAdapter: {'✅ PASS' if example_result else '❌ FAIL'}")
    
    for adaptor, result in adaptor_results.items():
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{adaptor:20} {status}")
    
    total_passed = sum(adaptor_results.values()) + (1 if example_result else 0)
    total_tests = len(adaptor_results) + 1
    
    print(f"\nOverall: {total_passed}/{total_tests} adaptors passed")
    
    if total_passed == total_tests:
        print("🎉 All adaptors have correct basic structure!")
        return 0
    else:
        print("⚠️  Some adaptors have structural issues.")
        return 1

if __name__ == "__main__":
    sys.exit(main())