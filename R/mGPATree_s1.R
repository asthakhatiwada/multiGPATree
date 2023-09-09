#' Implement the Stage 1 of the multiGPATree. Method
#'
#' This function will implement the stage 1 of the multiGPATree. method.
#'
#' @author  Aastha Khatiwada
#'
#' @param gwasPval A matrix of M X 2 dimension where M is the number of SNPs. The first column contains the SNP id and is labeled 'SNPid' and the second column contains the GWAS association p-values and is called P1. Values in P1 must be between 0 and 1.
#' @param annMat A matrix of binary annotations, where row and column correspond to SNPs and annotations, respectively.
#' @param initAlpha Initial value for alpha estimate. Default is 0.1.
#'
#' @import stats
#' @return This function returns a \code{List} including:
#' \itemize{
#' \item num_iter_stage1: number of iterations taken for Stage 1 of multiGPATree. to converge.
#' \item pi_stage1: predicted posterior probability of being a non-null SNP in Stage 1 of multiGPATree. Method.
#' \item alpha_stage1: estimated alpha of multiGPATree. Method.
#' \item beta_stage1: beta parameters from the linear model fitted at convergence of Stage 1 of multiGPATree. Method.
#' }
#'
#' @export

mGPATree_s1 <- function(gwasPval, annMat, initAlpha){

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
  message( "Info: Processing Stage 1 of multiGPATree...")
  M <- nrow(gwasPval) # number of SNPs
  nAnn <- ncol(annMat) # number of annotations
  nGWAS <- ncol(gwasPval) # number of GWAS
  
  # annMat <- as.data.frame(annMat)
  for (i in 1:ncol(annMat)) {  annMat[, i] <- as.factor(annMat[, i])   }
  
  # define pi matrix for multiple GWAS datas
  
  binaryList <- vector( "list", nGWAS )
  for ( k in 1:nGWAS ) {
    binaryList[[k]] <- c( 0, 1 )
  }
  
  binaryMat <- as.matrix(expand.grid( binaryList ))
  nComp <- nrow(binaryMat)
  combVec <- apply( binaryMat, 1, function(bm) paste( bm, collapse="" ) )
  piMat <- matrix(1/nComp, M, nComp)
  colnames(piMat) <- combVec
  
  
  # initialization of parameters
  # initial values for EM ####
  
  verbose <- TRUE
  iter_max <- 20000
  iter_cur <- 2
  
  # storing
  
  alphaMat <- matrix(NA, nrow = iter_max, ncol = nGWAS)
  colnames(alphaMat) <- paste('alpha', 1:nGWAS, sep = '')
  licList <- rep( NA, iter_max )
  lcList <- rep(NA, iter_max)
  
  # initial stores
  
  alphaMat[1, ] <- rep(initAlpha, nGWAS)
  licList[1] <- -10000000
  lcList[1] <- -10000000
  

  # EM start ####

  options(digits = 20)

  while( iter_cur < iter_max) {


    # E-step 
    
    # gwas pvalue distributions
    
    distList <- list()
    for (k in 1:nGWAS) {
      nullgrp <- rep(1, M)
      nonnullgrp <- alphaMat[iter_cur-1, k] * gwasPval[, k]^(alphaMat[iter_cur-1, k]-1)
      distList[[k]] <- cbind(nullgrp, nonnullgrp)
    }
    
    # create combination of the products of the distributions
    gMat <- matrix(NA, nrow = M, ncol = nrow(binaryMat))
    for (i in 1:nrow(binaryMat)) {
      # i = 1 
      gvec <- 1
      for (j in 1:nGWAS) {
        gvec <- gvec * distList[[j]][ , binaryMat[i, j] + 1]
      }
      gMat[, i] <- gvec
    }
    colnames(gMat) <- combVec
    
    g1Mat <- piMat * gMat
    
    ziMat <- g1Mat/rowSums(g1Mat)
    colnames(ziMat) <- combVec
    
    ziMat[ ziMat < 10e-10 ] <- 10e-10 # change this to 0.1/0.25
    ziMat[ ziMat > 1-10e-10 ] <- 1-10e-10
    # M-step
    
    # update piMat
    
    reg_tree <- lm(ziMat ~ annMat)
    piMat <- predict(reg_tree)
    piMat[ piMat < 10e-10 ] <- 10e-10 # change this to 0.1/0.25
    piMat[ piMat > 1 - 10e-10 ] <- 1 - 10e-10
    # update alpha
    alphaVec <- c()
    for (i in 1:nGWAS) {
      which_row <- which(binaryMat[, i] == 1)
      alphaVec[i] = max( min(- sum(rowSums(ziMat[, which_row]))/sum(rowSums(ziMat[, which_row])*log(gwasPval[,i])), 0.999), 0.001)
    }
    
    # Incomplete log-Likelihood
    
    lic <- sum(log(rowSums(g1Mat)))
    
    # Complete data log-likelihood
    
    lc <- sum(rowSums(ziMat * log(g1Mat)))
    
    # if ( verbose == TRUE ) {
    #   print(paste('iter_cur: ', iter_cur))
    #   # print(paste('cpTry: ', cpTry))
    #   print(paste('alpha', 1:nGWAS, ': ', alphaVec, sep = ''))
    #   print(paste('lic: ', lic))
    #   print(paste('lic in previous iteration: ', licList[iter_cur-1]))
    #   print(paste('lc: ', lc))
    #   print(paste('lc in previous iteration: ', lcList[iter_cur-1]))
    # }
    
    # store values
    alphaMat[iter_cur, ] <- alphaVec
    licList[iter_cur] <- lic
    lcList[iter_cur] <- lc
    
    # stopping rules
    
    if ( abs( licList[iter_cur] - licList[iter_cur-1] ) <= 1e-4 &
         sum(abs( alphaMat[iter_cur, ] - alphaMat[iter_cur-1, ] ) <= 1e-4) == nGWAS ){
      message('Info: Incomplete log-likelihood and alpha converged in Stage 1 of multiGPATree.')
      break
    }
    
    iter_cur <- iter_cur + 1

  }


  # list of variables to output
  message( "Info: Stage 1 of multiGPATree completed.")
  mEM_s1 <- list()
  mEM_s1$numIterConvergence <- iter_cur
  mEM_s1$alpha <- alphaMat[iter_cur, ]
  mEM_s1$Z <- ziMat
  mEM_s1$pi <- piMat
  mEM_s1$licVec <- licList[1:iter_cur]
  mEM_s1$lcVec <- lcList[1:iter_cur]
  mEM_s1$annMat <- annMat
  mEM_s1$gwasPval <- gwasPval
  
  
  return(mEM_s1)

}

# fit1 <- mGPATree_s1(simdata$gwasPval, simdata$annMat, initAlpha = 0.1)

