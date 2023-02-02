#' Generate simulation data for multiple phenotypes
#'
#' This function will simulate annotation data and p-values for multiple phenotypes to apply the multGPATree method.
#'
#' @author  Aastha Khatiwada
#'
#' @param M a numeric value that specifies the number of SNPs to be included in the simulation data.
#' @param nGWAS a numeric value that specifies the number of GWAS phenotypes to be included in the simulation.
#' @param nAnn a numeric value that specifies the number of functional annotations to be included in the simulation data.
#' @param percent_ones_ann a numeric value between 0 and 1 that specifies the percent of annotated SNPs for all functional annotations.
#' @param percent_overlap_ann a numeric value between 0 and 1 that specifies the percent of overlap between annotated SNPs for A1-A2 and A3-A4 and A5-A6.
#' @param percent_overlap_gwas a numeric value between 0 and 1 that specifies the percent of overlap P1 and P2.
#' @param trueAlphaVec a numeric vector between 0 and 1 to be used as the shape parameter for Beta distribution for GWAS phenotype p-values.
#'
#' @return A list containing the annotation dataframe, annMat, a dataframe containing GWAS association p-values, gwasPval.
#'
#'
#'


data_simulation <- function(M, nGWAS, nAnn, percent_ones_ann, percent_overlap_ann, percent_overlap_gwas, trueAlphaVec){
  
  # M = 1000
  # nGWAS = 2
  # nAnn = 10
  # percent_ones_ann = 0.20 # 10, 15, 20
  # percent_overlap_ann = 0.35 # >=35 percent
  # percent_overlap_gwas = 0.75 # >= 25, 50, 75
  # trueAlphaVec = c(0.2, 0.2)
  
  
  # start here
  
  freq_ones <- ceiling(percent_ones_ann * M)
  freq_overlap_ann <- ceiling(percent_overlap_ann * freq_ones)
  freq_overlap_gwas <- ceiling(percent_overlap_gwas * freq_overlap_ann)
  A1 <- c( rep(1, freq_ones), 
           rep(0, M-freq_ones) )
  A2 <- c( rep(0, freq_ones-freq_overlap_ann), 
           rep(1, freq_ones), 
           rep(0, M-(2*freq_ones-freq_overlap_ann)) )
  A3 <- c( rep(0, 2*freq_ones-freq_overlap_ann), 
           rep(1, freq_ones), 
           rep(0, M-(3*freq_ones-freq_overlap_ann)) )
  A4 <- c( rep(0, 3*freq_ones-2*freq_overlap_ann), 
           rep(1, freq_ones), 
           rep(0, M-(4*freq_ones-2*freq_overlap_ann)) )
  A5 <- c( rep(0, 4*freq_ones-2*freq_overlap_ann), 
           rep(1, freq_ones-freq_overlap_gwas), 
           rep(0, M - (5*freq_ones-2*freq_overlap_ann-freq_overlap_gwas) - freq_overlap_gwas ), 
           rep(1, freq_overlap_gwas ) )
  A6 <- c( rep(0, 5*freq_ones-3*freq_overlap_ann ), 
           rep(1, freq_ones-freq_overlap_gwas),
           rep(0, M-(6*freq_ones-3*freq_overlap_ann-freq_overlap_gwas)-freq_overlap_gwas ),
           rep(1, freq_overlap_gwas) )
  A <- cbind(A1, A2, A3, A4, A5, A6)
  # heatmap(as.matrix(A), Colv = NA, Rowv = NA, scale = 'column')
  
  if (nAnn >6){
    for (i in 7:nAnn) {
      prop_bin <- runif(M, 0.1, 0.3)
      new_ann <- rbinom(M, 1, prop_bin)
      A <- cbind(A, new_ann)
      colnames(A)[i] <- paste('A',i,sep = '')
    }
  }
  
  # generating latent Z
  
  binaryList <- vector( "list", nGWAS )
  for ( k in 1:nGWAS ) {
    binaryList[[k]] <- c( 0, 1 )
  }
  
  binaryMat <- as.matrix(expand.grid( binaryList ))
  nComp <- nrow(binaryMat)
  combVec <- apply( binaryMat, 1, function(bm) paste( bm, collapse="" ) )
  Zmat <- matrix(0, M, nComp)
  colnames(Zmat) <- combVec
  
  Z <- rep(0, M)
  Z[A1 == 1 & A2 == 1] <- 1
  Z[A3 == 1 & A4 == 1] <- 2
  Z[(M-freq_overlap_gwas+1):M] <- 3
  
  for (i in 1:ncol(Zmat)) {
    Zmat[Z==i-1, i] <- 1
  }
  

  alpha <- trueAlphaVec
  gwasPval <- matrix(runif(nGWAS*M, 0, 1), nrow = M, ncol = nGWAS)
  colnames(gwasPval) <- paste('P', 1:nGWAS, sep = '')
  gwasPval[ which(Z==1 | Z==3), 1 ] <- rbeta( length(which(Z==1 | Z==3)), alpha[1], 1) 
  gwasPval[ which(Z==2 | Z==3), 2 ] <- rbeta( length(which(Z==2 | Z==3)), alpha[2], 1) 
  rownames(A) <- rownames(gwasPval) <- paste('SNP_', 1:M, sep = '')

  Ztrue <- matrix(0, nrow = M, ncol = nGWAS)
  Ztrue[Z==1 | Z==3, 1] <- 1
  Ztrue[Z==2 | Z==3, 2] <- 1
  colnames(Ztrue) <- colnames(gwasPval)
  
  return(list(annMat = as.matrix(A),
              gwasPval = gwasPval #,
              # Ztrue = Ztrue,
              # zMat = Zmat
              ))

}



# Example
# alpha <- c(0.3, 0.4)
# simdata <- data_simulation(M = 1000,
#                       nGWAS = 2,
#                       nAnn = 10,
#                       percent_ones_ann = 0.15,
#                       percent_overlap_ann = 0.5,
#                       percent_overlap_gwas = 1,
#                       trueAlphaVec = alpha)
# save(simdata, file = "/Users/khatiwadaa/Desktop/multiGPATree/data/simdata.RData")
# checkheat <- as.matrix(cbind(test$zMat, test$annMat))
# head(checkheat)
# heatmap(checkheat)
# heatmap(checkheat, Colv = NA, Rowv = NA, scale = 'column')
# 
# par(mfrow=c(1,2))
# for (i in 1:length(alpha)) {
#   hist(test$gwasPval[, i], breaks = 100, xlab = paste('P', i, sep = ''),
#        main = paste('Beta(', alpha[i], ', 1)', sep = ''))
# }
