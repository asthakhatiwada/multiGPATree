## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  prompt = TRUE,
  comment = ""
)

## ---- message=FALSE, warning=FALSE--------------------------------------------
# install.packages("devtools")
# library(devtools)
# devtools::install_github("cran/mvpart")
# library(mvpart)
# devtools::install_github("asthakhatiwada/multiGPATree")
library(multiGPATree)

## -----------------------------------------------------------------------------
data("simdata")
class(simdata)
names(simdata)
dim(simdata$gwasPval)
dim(simdata$annMat)
head(simdata$gwasPval)
head(simdata$annMat)

## ---- message=FALSE-----------------------------------------------------------
fit.mGPATree <- multiGPATree(gwasPval = simdata$gwasPval,
                            annMat = simdata$annMat,
                            initAlpha = 0.1,
                            cpTry = 0.005)
fit.mGPATree

## ---- eval=TRUE, include=TRUE, message=FALSE, warning=FALSE-------------------
fit.mGPATree.pruned <- prune(fit.mGPATree@fit$P1_P2, cp = 0.20)
fit.mGPATree.pruned

## ---- eval=TRUE, include=TRUE, message=FALSE, warning=FALSE-------------------
plot(fit.mGPATree@fit$P1_P2)

## ---- eval=TRUE, include=TRUE, message=FALSE, warning=FALSE-------------------
leaf(fit.mGPATree@fit$P1_P2)

## ---- eval=TRUE, include=TRUE, message=FALSE, warning=FALSE-------------------
assoc.mGPATree <- multiGPATree::assoc(fit.mGPATree@fit$P1_P2,
                                      FDR = 0.01, 
                                      fdrControl="global")
head(assoc.mGPATree)
table(assoc.mGPATree$P1_P2)
table(assoc.mGPATree$P1_P2, assoc.mGPATree$leaf)
table(assoc.mGPATree$P1)
table(assoc.mGPATree$P1, assoc.mGPATree$leaf)
table(assoc.mGPATree$P2)
table(assoc.mGPATree$P2, assoc.mGPATree$leaf)

table(assoc.mGPATree$P1, assoc.mGPATree$P2)
table(assoc.mGPATree$P1, assoc.mGPATree$P2, assoc.mGPATree$P1_P2)

