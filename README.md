# multiGPATree
multiGPATree is a statistical approach to integrate GWAS summary statistics and functional annotation information for more than one trait within a unified framework. Specifically, by combining a multivariate decision tree algorithm with a hierarchical modeling framework, multiGPATree simultaneously implements association mapping for one or more trait and identifies key combinations of functional annotations related to one or more disease risk-associated SNPs. 

# Installation
To install multiGPATree, use the following commands. Please ensure to download the 'mvpart' package prior to installing the multiGPATree package.

```{r}
#install.packages("devtools")
library(devtools)
install_github("cran/mvpart")
install_github("asthakhatiwada/multiGPATree")
library(multiGPATree)
package?multiGPATree
```

# Demo
Please review our [package vignette](https://asthakhatiwada.github.io/multiGPATree/vignettes/multiGPATree_vignette_html.html) for demonstration and overview of the functions included in the multiGPATree package.


