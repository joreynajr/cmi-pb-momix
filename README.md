## Benchmark of multi-omics joint Dimensionality Reduction (jDR) approaches in CMI-PB
I have taken the original [momix-notebook repository](https://github.com/ComputationalSystemsBiology/momix-notebook)
and investigated 5 representative jDR approaches. As part of this analysis I take a look at: 
1. ability to identify factors associated with task data  
2. metagenes associated with biological annotations (Reactome, GO, Hallmarks)
 
The benchmarked methods are:
* [Integrative NMF (intNMF)](https://cran.r-project.org/web/packages/IntNMF/index.html) 
* [Joint and individual variation explained (JIVE)](https://cran.r-project.org/web/packages/r.jive/index.html) 
* [Multiple co-inertia analysis (MCIA)](https://bioconductor.org/packages/release/bioc/html/omicade4.html) 
* [Regularized Generalized Canonical Correlation Analysis (RGCCA)](https://cran.r-project.org/web/packages/RGCCA/index.html) 
* [matrix-tri-factorization (scikit-fusion)](https://github.com/marinkaz/scikit-fusion) 

Due to installation or running difficulties the following methods have been skipped (as of 2021.08.09):
* [iCluster](https://cran.r-project.org/web/packages/iCluster/index.html)
* [Multi-Omics Factor Analysis (MOFA)](https://github.com/bioFAM/MOFA)
* [Multi-Study Factor Analysis (MSFA)](https://github.com/rdevito/MSFA) 
* [tensorial Independent Component Analysis (tICA)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1455-8)

## Input data

The input data is derived from the CMI-PB datatables.

##  Cite momix
The preprint describing momix is available in BioRxiv
https://www.biorxiv.org/content/10.1101/2020.01.14.905760v1
