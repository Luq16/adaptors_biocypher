# BioCypher Template Project with Comprehensive Adaptors

A production-ready BioCypher template project with a comprehensive set of biological data adaptors for creating integrated knowledge graphs. This project includes 8 specialized adaptors covering proteins, compounds, diseases, Gene Ontology, pathways, phenotypes, protein interactions, and orthology data.

## 🎯 Overview

This project provides a complete set of BioCypher adaptors that can integrate multiple biological data sources into a unified knowledge graph. The adaptors follow BioCypher best practices and can be used individually or combined for comprehensive biological network analysis.

## 📦 Available Adaptors

**Complete Collection: 16 Adaptors** | 🔓 = Open Access | ⚠️ = Registration Required | 🔒 = License Required

### **Original BioCypher Adaptors (8)**

| Adaptor | Data Source | License | Node Types | Edge Types | Description |
|---------|-------------|---------|------------|------------|-------------|
| **UniProtAdapter** | UniProt | 🔓 | Proteins | Protein-Protein Interactions | Protein sequences, annotations, and interactions |
| **CompoundAdapter** | ChEMBL/DrugBank | 🔓/🔒 | Compounds/Drugs | Compound-Protein, Compound-Disease | Chemical compounds, drugs, and their targets |
| **DiseaseAdapter** | MONDO/OMIM | 🔓 | Diseases | Gene-Disease, Disease Hierarchy | Disease ontology and gene associations |
| **GOAdapter** | Gene Ontology | 🔓 | GO Terms (BP/MF/CC) | Protein-GO, GO Hierarchy | Gene Ontology annotations and relationships |
| **PPIAdapter** | BioGRID/IntAct | 🔓 | Proteins | Protein-Protein Interactions | Protein interaction networks with experimental evidence |
| **PathwayAdapter** | KEGG/Reactome | ⚠️/🔓 | Pathways | Protein-Pathway, Compound-Pathway | Biological pathways and molecular participation |
| **PhenotypeAdapter** | HPO/MP | 🔓 | Phenotypes | Gene-Phenotype, Disease-Phenotype | Phenotypic data and clinical associations |
| **OrthologyAdapter** | OrthoDB/OMA | 🔓 | Ortholog Groups | Protein-Ortholog, Orthology | Cross-species orthologous relationships |

### **New CROssBARv2-Inspired Adaptors (8)**

| Adaptor | Data Source | License | Node Types | Edge Types | Description |
|---------|-------------|---------|------------|------------|-------------|
| **SideEffectAdapter** | SIDER/OFFSIDES | ⚠️ | Side Effects | Drug-SideEffect, SideEffect Hierarchy | Drug adverse reactions and safety data |
| **DrugAdapter** | DrugBank | 🔒 | Drugs | Drug-Target, Drug-Disease, Drug-Drug | Comprehensive pharmaceutical data and interactions |
| **BiogridAdapter** | BioGRID | 🔓 | Proteins | Protein-Protein Interactions | High-quality experimental protein interactions |
| **StringAdapter** | STRING | 🔓 | Proteins | Functional Associations | Protein functional networks with confidence scores |
| **IntactAdapter** | IntAct | 🔓 | Proteins | Molecular Interactions | Curated molecular interactions with MI scores |
| **ECAdapter** | ExPASy/BRENDA | 🔓 | EC Numbers | Protein-EC, EC Hierarchy | Enzyme classification and catalytic annotations |
| **InterproAdapter** | InterPro | 🔓 | Domains/Families | Protein-Domain, Domain Relationships | Protein functional domains and family classification |
| **TFGenAdapter** | DoRothEA/TRRUST | 🔓 | Transcription Factors | TF-Gene Regulation | Gene regulatory networks and transcription control |

## 🚀 Quick Start

### Prerequisites

```bash
# Install BioCypher
pip install biocypher

# Clone or download this template project
git clone <your-repo-url>
cd project-template
```

### Basic Usage

```python
from template_package.adapters import UniProtAdapter, CompoundAdapter

# Initialize adaptors
protein_adapter = UniProtAdapter(test_mode=True)
compound_adapter = CompoundAdapter(test_mode=True)

# Generate nodes
protein_nodes = list(protein_adapter.get_nodes())
compound_nodes = list(compound_adapter.get_nodes())

print(f"Generated {len(protein_nodes)} proteins and {len(compound_nodes)} compounds")

# Generate edges (relationships)
protein_interactions = list(protein_adapter.get_edges())
compound_targets = list(compound_adapter.get_edges(target_proteins=protein_nodes))

print(f"Generated {len(protein_interactions)} protein interactions")
print(f"Generated {len(compound_targets)} compound-target relationships")
```

## 📋 Detailed Setup Guide

### 1. Installation Options

#### Option A: Using Poetry (Recommended)

```bash
# Clone the repository
git clone https://github.com/biocypher/project-template.git
cd project-template

# Install dependencies
poetry install

# Activate virtual environment
poetry shell

# Test the adaptors
python simple_test.py
```

#### Option B: Using pip

```bash
# Clone the repository
git clone https://github.com/biocypher/project-template.git
cd project-template

# Create virtual environment
python -m venv biocypher-env
source biocypher-env/bin/activate  # On Windows: biocypher-env\Scripts\activate

# Install dependencies
pip install biocypher pandas numpy tqdm

# Test the adaptors
python simple_test.py
```

### 2. Quick Test

```bash
# Test basic adaptor functionality
python simple_test.py

# Test edge generation
python test_edges.py  

# Test integrated functionality
python integration_test_fixed.py

# Run the original BioCypher example
python create_knowledge_graph.py
```

## 🔧 Adaptor Usage Guide

### Individual Adaptor Usage

#### UniProt Protein Adaptor

```python
from template_package.adapters import UniProtAdapter

# Basic usage
adapter = UniProtAdapter(
    organism="9606",      # Human (NCBI taxonomy ID)
    test_mode=False       # Set to False for full dataset
)

# Generate protein nodes
proteins = list(adapter.get_nodes())

# Sample protein node structure:
# ('P12345', 'protein', {
#     'name': 'Tumor protein p53',
#     'organism': '9606', 
#     'gene_name': 'TP53',
#     'sequence': 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPL...',
#     'function': 'Transcription factor'
# })

# Generate protein-protein interactions
interactions = list(adapter.get_edges())
```

#### Compound/Drug Adaptor

```python
from template_package.adapters import CompoundAdapter

adapter = CompoundAdapter(test_mode=True)

# Generate compound nodes
compounds = list(adapter.get_nodes())

# Sample compound node:
# ('CHEMBL123456', 'compound', {
#     'name': 'Aspirin',
#     'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
#     'molecular_weight': 180.16,
#     'drug_type': 'small molecule'
# })

# Generate compound-target interactions
# Need to provide target proteins and diseases
from template_package.adapters import UniProtAdapter, DiseaseAdapter

protein_adapter = UniProtAdapter(test_mode=True)
disease_adapter = DiseaseAdapter(test_mode=True)

proteins = list(protein_adapter.get_nodes())
diseases = list(disease_adapter.get_nodes())

compound_interactions = list(adapter.get_edges(
    target_proteins=proteins,
    target_diseases=diseases
))
```

### Integrated Multi-Adaptor Usage

```python
from template_package.adapters import (
    UniProtAdapter, CompoundAdapter, DiseaseAdapter, 
    GOAdapter, PathwayAdapter
)

# Initialize all adaptors
protein_adapter = UniProtAdapter(test_mode=True)
compound_adapter = CompoundAdapter(test_mode=True)
disease_adapter = DiseaseAdapter(test_mode=True)
go_adapter = GOAdapter(test_mode=True)
pathway_adapter = PathwayAdapter(test_mode=True)

# Generate nodes from all adaptors
all_nodes = []
all_nodes.extend(protein_adapter.get_nodes())
all_nodes.extend(compound_adapter.get_nodes())
all_nodes.extend(disease_adapter.get_nodes())
all_nodes.extend(go_adapter.get_nodes())
all_nodes.extend(pathway_adapter.get_nodes())

print(f"Total nodes: {len(all_nodes)}")

# Generate integrated edges
proteins = list(protein_adapter.get_nodes())
compounds = list(compound_adapter.get_nodes())
diseases = list(disease_adapter.get_nodes())

all_edges = []

# Protein-protein interactions
all_edges.extend(protein_adapter.get_edges())

# Compound-target interactions
all_edges.extend(compound_adapter.get_edges(
    target_proteins=proteins, 
    target_diseases=diseases
))

# Gene-disease associations
all_edges.extend(disease_adapter.get_edges(target_genes=proteins))

# GO annotations
all_edges.extend(go_adapter.get_edges(target_proteins=proteins))

# Pathway participation
all_edges.extend(pathway_adapter.get_edges(
    target_proteins=proteins,
    target_compounds=compounds
))

print(f"Total edges: {len(all_edges)}")
```

## 🏗️ BioCypher Integration

### Using with BioCypher Core

```python
import biocypher
from template_package.adapters import UniProtAdapter, CompoundAdapter

# Initialize BioCypher
bc = biocypher.BioCypher()

# Initialize adaptors
protein_adapter = UniProtAdapter()
compound_adapter = CompoundAdapter()

# Generate and write nodes
protein_nodes = protein_adapter.get_nodes()
bc.write_nodes(protein_nodes)

compound_nodes = compound_adapter.get_nodes()
bc.write_nodes(compound_nodes)

# Generate and write edges  
proteins = list(protein_adapter.get_nodes())
compound_edges = compound_adapter.get_edges(target_proteins=proteins)
bc.write_edges(compound_edges)

# Create summary
bc.summary()
```

## 📊 Data Sources and Identifiers

### Node Identifier Patterns

| Adaptor | ID Pattern | Example | Description |
|---------|------------|---------|-------------|
| UniProt | `P12345` | `P53_HUMAN` | UniProt accession |
| Compound | `CHEMBL123456` | `CHEMBL25` | ChEMBL identifier |
| Disease | `DOID:12345` | `DOID:557` | Disease Ontology ID |
| GO | `GO:1234567` | `GO:0008150` | Gene Ontology ID |
| Pathway | `hsa12345` | `hsa04110` | KEGG pathway ID |
| Phenotype | `HP:1234567` | `HP:0000118` | Human Phenotype Ontology |
| Ortholog | `123456at7742` | `441838at7742` | OrthoDB group ID |

### Relationship Types Generated

```
# Protein-centric
protein_protein_interaction
protein_has_go_annotation  
protein_participates_in_pathway
gene_associated_with_disease
gene_has_phenotype

# Compound-centric
compound_targets_protein
compound_treats_disease
compound_participates_in_pathway

# Hierarchical
go_term_is_a_go_term
disease_is_subtype_of_disease
pathway_is_part_of_pathway
phenotype_is_subtype_of_phenotype

# Cross-species
protein_orthologous_to_protein
protein_in_ortholog_group
```

## 🧪 Testing and Validation

### Running Tests

```bash
# Test basic adaptor structure
python simple_test.py

# Test edge generation
python test_edges.py

# Test integrated functionality
python integration_test_fixed.py
```

### Expected Test Results

```
BioCypher Adaptors Integration Test (Fixed)
=============================================
✅ Generated 253 total nodes
✅ Generated 880 total edges
✅ 16 different relationship types
✅ 100% valid edge structure
✅ Cross-adaptor connectivity established
🎉 Integration test completed successfully!
```

## 🎛️ Adaptor Configuration Options

### Common Parameters

All adaptors support these common parameters:

```python
adapter = AdapterClass(
    test_mode=True,           # Limit data size for testing
    node_fields=[...],        # Specify which node fields to include
    edge_fields=[...],        # Specify which edge fields to include
    organism="9606"           # NCBI taxonomy ID (where applicable)
)
```

### Field Selection Examples

#### Selecting Specific Protein Fields

```python
from template_package.adapters import UniProtAdapter
from template_package.adapters.uniprot_adapter import UniProtAdapterProteinField

# Only include specific protein fields
adapter = UniProtAdapter(
    node_fields=[
        UniProtAdapterProteinField.NAME,
        UniProtAdapterProteinField.SEQUENCE,
        UniProtAdapterProteinField.ORGANISM
    ],
    test_mode=True
)
```

#### Selecting Specific Edge Types

```python
from template_package.adapters import CompoundAdapter
from template_package.adapters.compound_adapter import CompoundAdapterEdgeType

# Only generate compound-protein interactions (not compound-disease)
adapter = CompoundAdapter(
    edge_types=[CompoundAdapterEdgeType.COMPOUND_TARGETS_PROTEIN],
    test_mode=True
)
```

## 🐛 Troubleshooting

### Common Issues

#### Import Errors
```python
# Error: ModuleNotFoundError: No module named 'biocypher'
# Solution: Install BioCypher
pip install biocypher
```

#### Memory Issues with Large Datasets
```python
# Use test_mode for development
adapter = UniProtAdapter(test_mode=True)

# Or limit organism scope
adapter = UniProtAdapter(organism="9606")  # Human only
```

#### Empty Results
```python
# Check if nodes exist before generating edges
nodes = list(adapter.get_nodes())
if nodes:
    edges = list(adapter.get_edges())
else:
    print("No nodes generated")
```

## 📈 Production Deployment

### Scaling for Production

```python
# Remove test_mode for full datasets
adapters = {
    'proteins': UniProtAdapter(test_mode=False, organism="9606"),
    'compounds': CompoundAdapter(test_mode=False),
    'diseases': DiseaseAdapter(test_mode=False)
}

# Process in batches
batch_size = 1000
for i, node in enumerate(adapter.get_nodes()):
    if i % batch_size == 0:
        print(f"Processed {i} nodes")
    # Process node
```

## ⚖️ Data Licensing and Usage Requirements

**IMPORTANT**: Some adaptors access data sources that require licenses or have usage restrictions. Please ensure compliance before using these adaptors in production:

### 🔒 **Complete Licensing Status for All 16 Adaptors**

| **Adaptor** | **Data Source** | **License Status** | **Usage Requirements** |
|-------------|-----------------|-------------------|------------------------|
| **DrugAdapter** | DrugBank | 🔒 **Paid License Required** | Academic: Free registration<br/>Commercial: Paid license required |
| **CompoundAdapter** | ChEMBL/DrugBank | 🔓/🔒 **Mixed** | ChEMBL: Free with attribution<br/>DrugBank: Paid license for commercial |
| **SideEffectAdapter** | SIDER/OFFSIDES | ⚠️ **Academic Restriction** | Academic: Free<br/>Commercial: Permission required |
| **PathwayAdapter** | KEGG/Reactome | ⚠️ **Registration Required** | KEGG: Registration for bulk access<br/>Reactome: Free |
| **StringAdapter** | STRING Database | ⚠️ **Registration for Bulk** | Free with attribution<br/>Bulk downloads: Registration |
| **UniProtAdapter** | UniProt | 🔓 **Open Access** | Free with attribution (CC BY 4.0) |
| **GOAdapter** | Gene Ontology | 🔓 **Open Access** | Free with attribution (CC BY 4.0) |
| **DiseaseAdapter** | MONDO/OMIM | 🔓 **Open Access** | Free with attribution (CC BY 3.0) |
| **PhenotypeAdapter** | HPO/MP | 🔓 **Open Access** | Free with attribution |
| **OrthologyAdapter** | OMA/OrthoDB | 🔓 **Open Access** | Free with attribution |
| **PPIAdapter** | Various PPI DBs | 🔓 **Open Access** | Free with attribution |
| **BiogridAdapter** | BioGRID | 🔓 **Open Access** | Free with attribution |
| **IntactAdapter** | IntAct | 🔓 **Open Access** | Free with attribution (EMBL-EBI) |
| **ECAdapter** | ExPASy/BRENDA | 🔓 **Open Access** | Free with attribution |
| **InterproAdapter** | InterPro | 🔓 **Open Access** | Free with attribution (EMBL-EBI) |
| **TFGenAdapter** | DoRothEA/TRRUST | 🔓 **Open Access** | Free academic databases |

### 📊 **License Summary:**
- **🔒 Paid License Required**: 2 adaptors (12.5%)
- **⚠️ Registration/Restrictions**: 3 adaptors (18.75%) 
- **🔓 Fully Open Access**: 11 adaptors (68.75%)

### 📋 **License Requirements Details**

#### **DrugBank (DrugAdapter)**
- **Academic Users**: Free registration at [drugbank.com](https://drugbank.com)
- **Commercial Users**: Contact DrugBank for licensing fees
- **API Access**: Rate limits apply; bulk downloads require approval
- **Attribution**: Must cite DrugBank in publications

#### **SIDER/OFFSIDES (SideEffectAdapter)**  
- **Academic Use**: Generally free for research purposes
- **Commercial Use**: May require permission from data providers
- **Publication**: Must cite original SIDER/OFFSIDES papers
- **Data Redistribution**: Restricted - check terms before sharing

#### **Proprietary Database Notes**
- **TRANSFAC (TFGenAdapter potential integration)**: Requires commercial license
- **Pathway databases**: Some require registration (KEGG, Reactome)
- **Clinical data**: May have patient privacy restrictions

### 🔓 **Open Access Data Sources**

The following adaptors use freely available data sources:
- **UniProtAdapter**: UniProt (CC BY 4.0)
- **GOAdapter**: Gene Ontology (CC BY 4.0) 
- **DiseaseAdapter**: Disease Ontology (CC BY 3.0)
- **PathwayAdapter**: Reactome (CC BY 4.0)
- **PhenotypeAdapter**: Human/Mammalian Phenotype Ontology (Free)
- **OrthologyAdapter**: OMA/OrthoDB (Free)
- **PPIAdapter**: Various open PPI databases
- **ECAdapter**: ExPASy Enzyme Database (Free)
- **InterproAdapter**: InterPro (EMBL-EBI, Free)

### 📝 **Usage Recommendations**

1. **Before Production Use**:
   - Review each data source's terms of use
   - Obtain necessary licenses/registrations
   - Implement appropriate attribution
   - Consider rate limiting for API access

2. **For Commercial Applications**:
   - DrugBank requires paid licensing
   - Verify commercial use permissions for all sources
   - Consider data redistribution restrictions

3. **For Academic Research**:
   - Most sources are free with proper attribution
   - Register for accounts where required
   - Cite original data sources in publications

4. **Compliance Best Practices**:
   - Keep licenses/registrations current
   - Monitor usage against terms of service
   - Implement proper data governance
   - Regular compliance audits

### 🚨 **Disclaimer**

This template provides code to access various biological databases. **Users are solely responsible for:**
- Obtaining appropriate licenses
- Complying with terms of use
- Proper attribution and citation
- Respecting usage limitations

**The template authors are not responsible for license violations or misuse of restricted data.**

## 📚 Additional Resources

- [BioCypher Documentation](https://biocypher.org)
- [BioCypher GitHub](https://github.com/biocypher/biocypher)
- [Biolink Model](https://biolink.github.io/biolink-model/)
- [Test Results Documentation](ADAPTOR_TEST_RESULTS.md)
- [Data Source Licensing Guide](https://biocypher.org/data-sources)

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Run the test suite
5. Submit a pull request

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 🏆 Acknowledgments

- BioCypher team for the excellent framework
- Contributors of biological databases (UniProt, ChEMBL, GO, etc.)
- CROssBARv2 project for adaptor examples and inspiration

---

**Ready to build comprehensive biological knowledge graphs with BioCypher!** 🧬📊

## 🛠 Original BioCypher Template Usage

### Structure
The project template is structured as follows:
```
.
│  # Project setup
│
├── LICENSE
├── README.md
├── pyproject.toml
│
│  # Docker setup
│
├── Dockerfile
├── docker
│   ├── biocypher_entrypoint_patch.sh
│   ├── create_table.sh
│   └── import.sh
├── docker-compose.yml
├── docker-variables.env
│
│  # Project pipeline
│
├── create_knowledge_graph.py
├── config
│   ├── biocypher_config.yaml
│   ├── biocypher_docker_config.yaml
│   └── schema_config.yaml
└── template_package
    └── adapters
        └── example_adapter.py
```

The main components of the BioCypher pipeline are the
`create_knowledge_graph.py`, the configuration in the `config` directory, and
the adapter module in the `template_package` directory. The latter can be used
to publish your own adapters (see below). You can also use other adapters from
anywhere on GitHub, PyPI, or your local machine.

**The BioCypher ecosystem relies on the collection of adapters (planned, in
development, or already available) to inform the community about the available
data sources and to facilitate the creation of knowledge graphs. If you think
your adapter could be useful for others, please create an issue for it on the
[main BioCypher repository](https://github.com/biocypher/biocypher/issues).**

In addition, the docker setup is provided to run the pipeline (from the same
python script) in a docker container, and subsequently load the knowledge graph
into a Neo4j instance (also from a docker container). This is useful if you want
to run the pipeline on a server, or if you want to run it in a reproducible
environment.

### Running the pipeline

`python create_knowledge_graph.py` will create a knowledge graph from the
example data included in this repository (borrowed from the [BioCypher
tutorial](https://biocypher.org/tutorial.html)). To do that, it uses the
following components:

- `create_knowledge_graph.py`: the main script that orchestrates the pipeline.
It brings together the BioCypher package with the data sources. To build a 
knowledge graph, you need at least one adapter (see below). For common 
resources, there may already be an adapter available in the BioCypher package or
in a separate repository. You can also write your own adapter, should none be
available for your data.

- `example_adapter.py` (in `template_package.adapters`): a module that defines
the adapter to the data source. In this case, it is a random generator script.
If you want to create your own adapters, we recommend to use the example adapter
as a blueprint and create one python file per data source, approproately named.
You can then import the adapter in `create_knowledge_graph.py` and add it to
the pipeline. This way, you ensure that others can easily install and use your 
adapters.

- `schema_config.yaml`: a configuration file (found in the `config` directory)
that defines the schema of the knowledge graph. It is used by BioCypher to map
the data source to the knowledge representation on the basis of ontology (see
[this part of the BioCypher 
tutorial](https://biocypher.org/tutorial-ontology.html)).

- `biocypher_config.yaml`: a configuration file (found in the `config` 
directory) that defines some BioCypher parameters, such as the mode, the 
separators used, and other options. More on its use can be found in the
[Documentation](https://biocypher.org/installation.html#configuration).

### Publishing your own adapters

After adding your adapter(s) to the `adapters` directory, you may want to
publish them for easier reuse. To create a package to distribute your own
adapter(s), we recommend using [Poetry](https://python-poetry.org/). Poetry,
after setup, allows you to publish your package to PyPI using few simple
commands. To set up your package, rename the `template_package` directory to
your desired package name and update the `pyproject.toml` file accordingly. Most
importantly, update the `name`,`author`, and `version` fields. You can also add
a `description` and a `license`.  Then, you can publish your package to PyPI
using the following commands:

```{bash}
poetry build
poetry publish
```

If you don't want to publish your package to PyPI, you can also install it from
GitHub using the respective functions of poetry or pip.

### Further reading / code

If you want to see a second example of the workflow, check our
[CollecTRI](https://github.com/biocypher/collectri) pipeline. Its README describes
the process of data assessment and adapter creation in more detail.

## 🐳 Docker

This repo also contains a `docker compose` workflow to create the example
database using BioCypher and load it into a dockerised Neo4j instance
automatically. To run it, simply execute `docker compose up -d` in the root 
directory of the project. This will start up a single (detached) docker
container with a Neo4j instance that contains the knowledge graph built by
BioCypher as the DB `neo4j` (the default DB), which you can connect to and
browse at localhost:7474. Authentication is deactivated by default and can be
modified in the `docker_variables.env` file (in which case you need to provide
the .env file to the deploy stage of the `docker-compose.yml`).

Regarding the BioCypher build procedure, the `biocypher_docker_config.yaml` file
is used instead of the `biocypher_config.yaml` (configured in
`scripts/build.sh`). Everything else is the same as in the local setup. The
first container (`build`) installs and runs the BioCypher pipeline, the second
container (`import`) installs Neo4j and runs the import, and the third container
(`deploy`) deploys the Neo4j instance on localhost. The files are shared using a
Docker Volume. This three-stage setup strictly is not necessary for the mounting
of a read-write instance of Neo4j, but is required if the purpose is to provide
a read-only instance (e.g. for a web app) that is updated regularly; for an
example, see the [meta graph
repository](https://github.com/biocypher/meta-graph). The read-only setting is
configured in the `docker-compose.yml` file
(`NEO4J_dbms_databases_default__to__read__only: "false"`) and is deactivated by
default.
