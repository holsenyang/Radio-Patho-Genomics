# Multi-omics Data Analysis Pipeline

A comprehensive bioinformatics pipeline for multi-omics data analysis, including co-expression network construction, functional enrichment, and pathway activity analysis.

## Project Overview

This pipeline performs integrated analysis of transcriptomic data across multiple datasets (Batch, TCGA, CPTAC) to identify biologically meaningful patterns and pathways relevant to the research focus. The analysis includes quality control, network analysis, module identification, functional annotation, and pathway activity assessment.

## Complete Analysis Pipeline

### 1. Data Preprocessing and Quality Control
- **Script**: `01_data_preprocessing.R`
- **Purpose**: Raw data loading, normalization, batch effect correction, and quality control
- **Input**: Raw expression matrices (TPM/FPKM counts)
- **Output**: Cleaned, normalized expression datasets ready for analysis

### 2. Co-expression Network Construction
- **Script**: `02_network_construction.R`
- **Purpose**: Weighted Gene Co-expression Network Analysis (WGCNA)
- **Key Steps**:
  - Soft threshold power selection
  - Network construction and module detection
  - Module-trait associations
  - Sample clustering and outlier detection
- **Output**: Gene modules, network topology, module-trait relationships

### 3. Hub Gene Identification
- **Script**: `03_hub_gene_analysis.R`
- **Purpose**: Identify key driver genes within modules
- **Key Steps**:
  - Module membership and gene significance calculation
  - Hub gene identification based on connectivity
  - Cross-dataset hub gene validation
  - Visualization of MM-GS relationships
- **Output**: Hub gene lists, validation across datasets

### 4. Module Preservation Analysis
- **Script**: `04_module_preservation.R`
- **Purpose**: Assess module reproducibility across datasets
- **Key Steps**:
  - Z-summary statistics calculation
  - Module preservation evaluation
  - Cross-dataset consistency assessment
- **Output**: Preservation statistics, reproducibility metrics

### 5. Functional Enrichment Analysis
- **Script**: `05_enrichment_analysis.R`
- **Purpose**: Biological interpretation of gene modules
- **Databases**: GO, KEGG, Reactome, WikiPathways, BioCarta, Hallmark
- **Key Steps**:
  - Pathway enrichment analysis
  - FDR correction and significance filtering
  - Cross-dataset enrichment reproducibility
  - Top pathway identification
- **Output**: Enriched pathways, significance tables, reproducibility analysis

### 6. Pathway Activity Analysis (GSVA)
- **Script**: `06_gsva_analysis.R`
- **Purpose**: Gene Set Variation Analysis for pathway activity
- **Key Steps**:
  - Pathway activity scoring using significant enriched pathways
  - Association with clinical features
  - Cross-dataset pathway activity consistency
  - Heatmap visualization
- **Output**: GSVA scores, pathway activity matrices, association results

## Dataset Integration

The pipeline integrates three main datasets:
- **Batch Data**: Internal cohort samples (n=36)
- **TCGA Data**: The Cancer Genome Atlas samples (n=407)  
- **CPTAC Data**: Clinical Proteomic Tumor Analysis Consortium samples (n=variable)

## Quick Start

### Option 1: Run Complete Pipeline
```r
source("run_analysis.R")