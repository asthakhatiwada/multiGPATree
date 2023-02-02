# multiGPATree
multiGPATree is a statistical approach to integrate GWAS summary statistics and functional annotation information for more than one trait within a unified framework. Specifically, by combining a multivariate decision tree algorithm with a hierarchical modeling framework, multiGPATree simultaneously implements association mapping for one or more trait and identifies key combinations of functional annotations related to one or more disease risk-associated SNPs. 

# Installation
To install the multiGPATree, use the following commands. Please ensure to download the 'mvpart' package prior to installing the multiGPATree package.

```{r}
#install.packages("devtools")
library(devtools)
install_github("cran/mvpart")
install_github("asthakhatiwada/multiGPATree")
```

# Usage
The 'multiGPATree_vignette.pdf' provides a framework for the step-by-step post-GWAS analysis using the 'multiGPATree' package. Help file generated using the following code also provides a good starting point, including some example command lines, to use the 'multiGPATree' package:

```{r}
library(multiGPATree)
package?multiGPATree
```

