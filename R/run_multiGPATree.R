#' Implement the multiGPATree Method
#'
#' This function will implement Stage 1 and 2 of the multiGPATree. method.
#'
#' @author  Aastha Khatiwada
#'
#' @param gwasPval A matrix of M X D dimension, where M is the number of SNPs and D is the number of phenotypes. The columns contains the GWAS association p-values for 1, .., D phenotypes. P-values must be between 0 and 1.
#' @param annMat A matrix of binary annotations, where row and column correspond to SNPs and annotations, respectively.
#' @param initAlpha Initial value for alpha estimate. Default is 0.1.
#' @param cpTry Complexity parameter (cp) value to be used. cpTry can be between 0 and 1 or NULL. Default is 0.001. When cpTry is NULL, multiGPATree will select the optimal cp to be used.
#' @param ncore number of cores to use. we recommend changing this if you are analyzing more than one GWAS pair
#' @details The multiGPATree() function fits the multiGPATree model. It requires atleast two GWAS and expects users to provide GWAS p-value to gwasPval and binary annotation data to annMat.
#' It is assumed that number of rows of matrix in gwasPval and annMat are equal and correspond to the same SNP.
#'
#' The assoc() function implements association mapping.
#'
#' The plot() function takes in an object of class multiGPATree and will plot the functional annotation tree from the multiGPATree model.
#'
#' The leaf() function takes in an object of class multiGPATree and will provide information regarding the functional annotations that are enriched (1) or not enriched (0) for SNPs in any leaf of the multiGPATree model plot.
#'
#' The prune() function takes in an object of class multiGPATree and a cp parameter and will prune the multiGPATree model result. This function can be useful when the tree obtained from multiGPATree model is larger.
#' @examples
#' \dontrun{
#' library(multiGPATree)
#'
#' # load multiGPATree example data
#' data(simdata)
#'
#' # fitting the multiGPATree model
#' fit <- multiGPATree(simdata$gwasPval, simdata$annMat)
#'
#' # get functional annotation information
#' leaf(fit)
#'
#' # association mapping
#' assoc.mgpatree <- assoc(fit, FDR = 0.01, fdrControl = 'global')
#'
#' # pruning the multiGPATree model fit
#' pruned.fit <- prune(fit, cp = 0.005)
#'
#' # plotting the multiGPATree model results
#' plot(fit)
#' plot(pruned.fit)
#' }
#' @return Contructs a multiGPATree class object
#' @importFrom methods new
#' @export


run_multiGPATree <- function(gwasPval, annMat, initAlpha = 0.1, cpTry = 0.001, ncore=2){
  
  # testdata <- multiGPATree:::data_simulation(M = 1000,
  #                                            nGWAS = 2,
  #                                            nAnn = 5,
  #                                            percent_ones = 0.1,
  #                                            percent_overlap_ann = 0.5,
  #                                            percent_overlap_gwas = 1,
  #                                            trueAlphaVec = c(0.4,0.4))
  # gwasPval = testdata$gwasPval
  # annMat = testdata$annMat
  # initAlpha = 0.1
  # cpTry = 0.001
  
  # start here
  GPATREE_S1 <- mGPATree_s1(gwasPval = gwasPval,
                            annMat = annMat,
                            initAlpha = initAlpha)
  
  GPATREE_S2 <- mGPATree_s2(gwasPval = gwasPval,
                            annMat = annMat,
                            alphaStage1 = GPATREE_S1$alpha,
                            initPi = GPATREE_S1$pi,
                            cpTry = cpTry)
  
  return_list <- list()
  return_list$numIterConvergence <- c('stage1' = GPATREE_S1$numIterConvergence, 
                                      'stage2' =  GPATREE_S2$numIterConvergence)
  return_list$alpha <- GPATREE_S1$alpha
  return_list$fit <- GPATREE_S2$fit
  return_list$fitSelectVar <- GPATREE_S2$fitSelectVar
  return_list$Z <- as.matrix(GPATREE_S2$Z)
  return_list$Zmarg <- as.matrix(GPATREE_S2$Zmarg) 
  colnames(return_list$Zmarg) <- colnames(gwasPval)
  rownames(return_list$Zmarg) <- rownames(gwasPval)
  return_list$pi <- as.matrix(GPATREE_S2$pi)
  return_list$licVec <- GPATREE_S2$licVec
  return_list$lcVec <- GPATREE_S2$lcVec
  # return(return_list)
  
  new( Class = "multiGPATree",
       fit = return_list)
}
