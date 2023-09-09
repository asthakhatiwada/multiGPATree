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
#' 
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
#' @import foreach doParallel
#' @importFrom utils combn 
#' @export


multiGPATree <- function(gwasPval, annMat, initAlpha = 0.1, cpTry = 0.001, ncore = 2){

  # testdata <- multiGPATree:::data_simulation(M = 1000,
  #                                            nGWAS = 2,
  #                                            nAnn = 10,
  #                                            percent_ones = 0.1,
  #                                            percent_overlap_ann = 0.5,
  #                                            percent_overlap_gwas = 1,
  #                                            trueAlphaVec = c(0.4,0.4))
  # head(testdata$annMat)
  # tdat <- multiGPATree:::data_simulation(M = 1000,
  #                                        nGWAS = 2,
  #                                        nAnn = 10,
  #                                        percent_ones = 0.1,
  #                                        percent_overlap_ann = 0.5,
  #                                        percent_overlap_gwas = 1,
  #                                        trueAlphaVec = c(0.3,0.3))
  # # # 
  # # # head(testdata$gwasPval)
  # # # class(testdata$gwasPval)
  # testdata$gwasPval <- as.matrix(cbind(testdata$gwasPval, tdat$gwasPval[, 1]))
  # colnames(testdata$gwasPval) <- c('P1', 'P2', 'P3')
  # # head(testdata$gwasPval)
  # 
  # 
  # gwasPval = testdata$gwasPval
  # annMat = testdata$annMat
  # initAlpha = 0.1
  # cpTry = 0.001
  # ncore <- detectCores()
  # head(annMat)
  
  
  #start here

  if ( initAlpha<= 0 | initAlpha >= 1 ) {
    stop( "Inappropriate value for 'initAlpha' argument. It should be between zero and one." )
  }

  if ( !is.matrix(gwasPval) ) {
    gwasPval <- as.matrix(gwasPval, colClasses=c("numeric"))
  }
  

  if ( !is.null(annMat) ) {
    if ( !is.matrix(annMat) ) {
      annMat <- as.matrix(annMat)
    }
  }

  if ( !is.null(annMat) ) {
    if ( nrow(gwasPval) != nrow(annMat) ) {
      stop( "Number of SNPs are different between p-value and annotation matrices. They should coincide. Please check your p-value and annotation matrices.")
    }
  }
  

  if ( any( gwasPval < 0 | gwasPval > 1 ) ) {
    stop( "Some p-values are smaller than zero or larger than one. p-value should be ranged between zero and one. Please check your p-value matrix." )
  }

  if ( !is.null(annMat) ) {
    if ( any( annMat != 0 & annMat != 1 ) ) {
      stop( "Some elements in annotation matrix has values other than zero or one. All the elements in annotation matrix should be either zero (not annotated) or one (annotated). Please check your annotation matrix." )
    }
  }

  if (!is.null(cpTry)) {
    if ( cpTry< 0 | cpTry > 1 ) {
      stop("Inappropriate value for complexity parameter (cp). cp should be between 0 and 1")
    }
  }


  nGWAS <- ncol(gwasPval)                 # No. of GWAS
  # nGWAS
  Npairs <- nGWAS*(nGWAS-1)/2            # No. of pairs
  Npairs
  if ( !is.null(annMat) ) {
    nAnn <- ncol(annMat)
  }

  message( "Info: Number of GWAS data: ", nGWAS )
  if (nGWAS == 1){
    message( "Info: data contains only one GWAS. Please use the GPATree package instead.")
  }
  
  if (nGWAS > 1){
    message( "Info: Number of GWAS pairs: ", Npairs )
  }
  
  if ( !is.null(annMat) ) {
    message( "Info: Number of annotations: ", nAnn )
  }
  message( "Info: Number of SNPs: ", nrow(gwasPval))
  
  
  if (!is.null(cpTry)) {
    message( "Info: complexity parameter (cp) to be used: ", cpTry )
  }
  if (is.null(cpTry)) {
    cpTry <- 0.001
    message( "Info: complexity parameter (cp) will be determined by multiGPATree.")
  }
  
  
  
  
  # if(nGWAS == 2){
  #   
  #   return_list <- run_multiGPATree(gwasPval = gwasPval,
  #                                   annMat = annMat, 
  #                                   initAlpha = initAlpha, 
  #                                   cpTry = cpTry)
  # }
  
  # gwasreturn <- gwasPval
  if (nGWAS >= 2){
    pair_comb <- combn(nGWAS, 2)
    datacall <- NULL
    registerDoParallel(cores = ncore, cl = ncore)
    return_list <- foreach(datacall = 1:Npairs, .combine = 'append') %dopar% {
      # datacall <- 1
      message( "Info: currently running pair ", datacall)
      gwasdata <- gwasPval[, pair_comb[, datacall]]
      name_list <- paste0(colnames(gwasdata), collapse = '_')
      result <- run_multiGPATree(gwasPval = gwasdata,
                                 annMat = annMat, 
                                 initAlpha = initAlpha, 
                                 cpTry = cpTry,
                                 ncore = ncore)
      # multigpatree_return <- list(result)
      multigpatree_return <- list(result)
      names(multigpatree_return) <- name_list
      # class(multigpatree_return) <- "multiGPATree"
      return(multigpatree_return)
    }
  }
  
  new( Class = "multiGPATreeMethod",
       fit = return_list,
       gwasPval = gwasPval,
       annMat = annMat )

}

