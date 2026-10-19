# Bioinformatics Pipeline: Mutual Information and Gene Distance Across Chromosomes and Breast Cancer Subtypes

This repository contains Bash and R scripts designed to automate genomic data processing, calculate statistical dependencies 
(Mutual Information) alongside physical gene distances, and generate comparative figures across human chromosomes and breast cancer subtypes.

---

## What Does This Pipeline Do?

The pipeline systematically iterates through human chromosomes (1 to 22, and X) and various biological subtypes to execute a multi-step analytical workflow:

1. **Gene Filtering (Data Preparation):** 
   - Uses `awk` to parse annotation datasets (e.g., BioMart exports) to isolate `protein_coding` genes per chromosome.
2. **Expression Matrix Subsetting:** 
   - Loops through specific breast cancer subtypes and control groups to extract matching expression profiles into dedicated tab-separated tables.
3. **Statistical Analysis & Distance Calculation:** 
   - Computes **Mutual Information (MI)**—evaluating statistical dependencies between genes—and integrates physical genomic distances.
4. **Figure Generation (Active Execution):** 
   - Invokes automated R scripts to process the filtered datasets and output publication-ready figures comparing Mutual Information versus gene distances and ordinal metrics.

---

The master script in bash calls the programs in R to do the tasks. 
