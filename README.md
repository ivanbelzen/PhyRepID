# PhyRepID pipeline

## Overview

The Phylogenomics Repeat Identification 
(PhyRepID) pipeline quantifies repeat evolution using comparative (phylo)genomics of orthologous groups (OGs). 

The pipeline consists of four components: 

 1. **Data collection:** Retrieve protein sequences and gene trees from ENSEMBL from orthologous groups based  on a human protein-coding gene and its relation to 13 other selected representative vertebrate species.
 2. **Repeat detection & optimization:** Detection of a broad spectrum of protein repeats consisting of either structural domains or motif sequences. Followed by improving the sensitivity and precision of detection by making OG-specific repeat models.
 4. **Phylogenetic tree reconciliation:** Inferring evolutionary events in the repeat region through phylogenetic comparison of repeat trees to gene trees.
 5. **Post-processing & duplication analysis:** Finally, to quantify repeat evolution and compare protein families, a PRD (protein repeat dupllication) score was derived from the post-ancestral duplications and the number of proteins in the OG.

The workflow is defined in Snakemake.
[Detailed documentation can be found in the docs](docs/pipeline_documentation.md)

## Dependencies

### Python
Version 2.7 with pip (9.0.1)

 *pipeline\_methods.py* contains configuration like paths and environmental variables, as well as functions used by other Python scripts

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

### Software

- HMMER3.1.b2 download and compile source code [http://hmmer.org/download.html  
](http://hmmer.org/download.html) 
BSD-3-Clause licensed.
 Finn, R.D., Clements, J. & Eddy, S.R., 2011. HMMER web server: interactive sequence similarity searching. Nucleic acids research, 39(Web Server issue), pp.W29–37.
 
- Pfam-A 31.0 models [ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/)
CC0 licensed.
Finn, R.D. et al., 2015. The Pfam protein families database: towards a more sustainablefuture. Nucleic acids research, 44(D1), pp.D279–D285.


-   MAFFT L-INS-I 
Katoh, K., 2005. MAFFT version 5: improvement in accuracy of multiple sequence alignment. Nucleic acids research, 33(2), pp.511–518    

  

 -   IQ-TREE 
Nguyen, L.-T. et al., 2014. IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies. Molecular biology and evolution, 32(1), pp.268–274.

UBoot ultrafast bootstrapping
Hoang, D.T. et al., 2018. UFBoot2: Improving the Ultrafast Bootstrap Approximation. Molecular biology and evolution, 35(2), pp.518–522.

ModelFinder 
Kalyaanamoorthy, S. et al., 2017. ModelFinder: fast model selection for accurate phylogenetic estimates. Nature methods, 14(6), pp.587–589.
  

-   Species tree manually constructed based on the ENSEMBL species tree (Zerbino et al. 2018) annotated with taxonomic identifiers for the nodes and divergence times as branch lengths with the help of TimeTree.org (Hedges et al. 2015)

Zerbino, D.R. et al., 2018. Ensembl 2018. Nucleic acids research, 46(D1), pp.D754–D761.

Hedges, S.B. et al., 2015. Tree of life reveals clock-like speciation and diversification. Molecular biology and evolution, 32(4), pp.835–845.  

-   TreeFix
Wu, Y.-C. et al., 2013. TreeFix: statistically informed gene tree error correction using species trees. Systematic biology, 62(1), pp.110–120.
  

-   ETE toolkit 
Huerta-Cepas, J., Serra, F. & Bork, P., 2016. ETE 3: Reconstruction, Analysis, and Visualization of Phylogenomic Data. Molecular biology and evolution, 33(6), pp.1635–1638.

### Input data:  

- Orthologous proteins were acquired from ENSEMBL Compara v91 (Zerbino et al. 2018) using a SPARQL query on the EBI endpoint (https://www.ebi.ac.uk/rdf/services/sparql) on 12-04-2018. Protein-coding genes homologous in homo sapiens (http://identifiers.org/taxonomy/9606) and 13 vertebrate species [[taxon id list]]
  
- The gene trees and protein sequences were retrieved via the ENSEMBL API ([http://rest.ensembl.org](http://rest.ensembl.org))
  
- The HGNC gene symbols were retrieved from ENSEMBL BioMart on 30-09-2018 from Human genes GRCh38.p12 (ENSEMBL release 93).

 
- A list of human genes under positive selection as measured by dN/dS since the divergence of vertebrates was acquired from Sebastien Moretti on 04-09-2018
Moretti, S. et al., 2014. Selectome update: quality control and computational improvements to a database of positive selection. Nucleic acids research, 42(Database issue), pp.D917–21.


- Z-scores for synonymous, non-synonymous and loss-of-function mutations were acquired for each human gene, as well as the probability of being intolerant for loss-of-function mutations from the ExAC database v0.3.1 
Karczewski, K.J. et al., 2016. The ExAC Browser: Displaying reference data information from over 60,000 exomes. Available at: http://dx.doi.org/10.1101/070581.


