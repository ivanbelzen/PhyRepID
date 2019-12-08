
## Introduction

The protein repeat evolution (PRE/PhyRepD) pipeline quantifies repeat evolution using comparative (phylo)genomics of orthologous groups (OGs). 

The pipeline consists of ~~three~~ four components: 

 1. **Repeat detection:** Detection of a broad spectrum of protein repeats consisting of either structural domains or short linear motif sequences.
 2. **Repeat optimization:** Improving the sensitivity and precision of detection by making OG-specific repeat models.
 3. **Tree reconciliation:** Inferring evolutionary events in the repeat region through phylogenetic comparison of repeat trees to gene trees.
 4. **Post-processing / duplication analysis:** Finally, to quantify repeat evolution and compare protein families, a duplication score was derived from the post-ancestral duplications and the number of proteins in the OG.

## Overview of the workflow

The workflow for components 1-3 is defined in two Snakemake files which link together Python and Bash scripts. Post-processing and analysis of the inferred repeat duplications, as well as comparisons with other datasets are done with various Python and R scripts. 

 - Data collection (Snakemake)
 - Snakemake: repeat detection
 - Scripts collection: repeat duplication analysis 


## Prerequisites

### Python
Version 2.7 with pip (9.0.1)
 *pipeline\_methods.py*
contains configuration like paths and environmental variables, as well as methods used by other Python scripts

**Data retrieval and processing**
SPARQLWrapper (1.8.1)
requests (2.18.4)
pandas (0.23.1)
numpy (1.13.1)

**Tree reconciliation** 
anaconda-client (1.2.2)
conda (4.3.27)
treefix (1.1.10)
ete3 (3.0.0b36)

**Analysis:**
matplotlib (2.2.2)
matplotlib-venn (0.11.5)
seaborn (0.9.0)
scipy (1.1.0)
scikit-learn (0.19.1)
goatools (0.8.4)

### Software
HMMer
TreeFix

