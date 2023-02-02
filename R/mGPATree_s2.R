#' Implement the multiGPATree Method for multiple phenotypes
#'
#' This function will implement the multiGPATree method for multiple phenotypes while leveraging pleiotropy.
#'
#' @author  Aastha Khatiwada
#' @import mvpart
#' @param gwasPval A matrix of M X D dimension, where M is the number of SNPs and D is the number of phenotypes. The columns contains the GWAS association p-values for the respective phenotypes. The pvalues must be between 0 and 1.
#' @param annMat A matrix of binary annotations, where row and column correspond to SNPs and annotations, respectively.
#' @param alphaStage1 Alpha estimated in stage 1
#' @param initPi final alpha at convergence of stage 1
#' @param cpTry Complexity parameter (cp) value to be used to build multivaraite CART model. cpTry can be between 0 and 1 or NULL. Default is 0.001. When cpTry is NULL, multiGPATree will select the optimal cp to be used.
#'
#' @return This function returns a \code{List} including:
#' \itemize{
#' \item numIterConvergence: number of iterations taken for multiGPATree to converge.
#' \item alphaVec: alpha estimated by multiGPATree.
#' \item fit: MVPART model selected by multiGPATree.
#' \item fitSelectVar: annotations included in MVPART tree.
#' \item ziMat: posterior probability of being a non-null SNP.
#' \item piMat: predicted posterior probability of of being a non-null SNP.
#' \item licVec: list of incomplete log likelihood for all iteration of multiGPATree Method.
#' \item lcVec: list of complete log likelihood for all iteration of multiGPATree Method.
#' }
#'


mGPATree_s2 <- function(gwasPval, annMat, alphaStage1, initPi, cpTry){
  
  
  # start here
  
  message( "Info: Processing Stage 2 of multiGPATree...")
  
  if ( !is.matrix(gwasPval) ) {
    gwasPval <- as.matrix(gwasPval)
  }
  
  if ( !is.null(annMat) ) {
    if ( !is.matrix(annMat) ) {
      annMat <- as.matrix(annMat)
    }
  }
  
  M <- nrow(gwasPval) # number of SNPs
  nAnn <- ncol(annMat) # number of annotations
  nGWAS <- ncol(gwasPval) # number of GWAS
  
  annMat <- as.data.frame(annMat)
  for (i in 1:ncol(annMat)) {  annMat[, i] <- as.factor(annMat[, i])   }
  
  # initialize piMat using piMat from Stage 1
  
  piMat <- initPi
  # piMat[ piMat < 0.01 ] <- 0.01
  # piMat[ piMat > 0.99 ] <- 0.99
  
  # hist(rowSums(piMat))
  
  # define binaryList for multiple GWAS datas
  
  binaryList <- vector( "list", nGWAS )
  for ( k in 1:nGWAS ) {
    binaryList[[k]] <- c( 0, 1 )
  }
  
  binaryMat <- as.matrix(expand.grid( binaryList ))
  nComp <- nrow(binaryMat)
  combVec <- apply( binaryMat, 1, function(bm) paste( bm, collapse="" ) )
  
  # cp try list if lc doesn't increase ####
  cp_log10_prune_list <- seq(-1, -5, by =-0.01) #suggested by Dr. Chung
  cp_log10_scale_list <- 10^cp_log10_prune_list
  
  # initialization of parameters
  # initial values for EM ####
  
  verbose <- TRUE
  iter_max <- 20000
  iter_cur <- 2
  
  # storing
  
  # alphaMat <- matrix(NA, nrow = iter_max, ncol = nGWAS)
  # colnames(alphaMat) <- paste('alpha', 1:nGWAS, sep = '')
  licList <- rep( NA, iter_max )
  lcList <- rep(NA, iter_max)
  
  # initial stores
  
  # alphaMat[1, ] <- rep(initAlpha, nGWAS)
  licList[1] <- -10000000
  lcList[1] <- -10000000
  
  # gwas pvalue distributions
  
  distList <- list()
  for (k in 1:nGWAS) {
    nullgrp <- rep(1, M)
    nonnullgrp <- alphaStage1[k] * gwasPval[, k]^(alphaStage1[k]-1)
    distList[[k]] <- cbind(nullgrp, nonnullgrp)
  }
  
  # create combination of the products of the distributions
  gMat <- matrix(NA, nrow = M, ncol = nrow(binaryMat))
  for (i in 1:nrow(binaryMat)) {
    gvec <- 1
    for (j in 1:nGWAS) {
      gvec <- gvec * distList[[j]][ , binaryMat[i, j] + 1]
    }
    gMat[, i] <- gvec
  }
  colnames(gMat) <- combVec
  
  
  while( iter_cur < iter_max ) {
    
    
    # E-step 
    
    
    g1Mat <- piMat * gMat
    ziMat <- g1Mat/rowSums(g1Mat)
    colnames(ziMat) <- combVec
    ziMat[ ziMat < 10e-10 ] <- 10e-10
    ziMat[ ziMat >= 1 - 10e-10 ] <- 1 - 10e-10
    
    # M-step
    
    # update piMat
    
    mvpart_res <- data.matrix(ziMat)
    mvpart_covariates <- as.data.frame(annMat)
    reg_full_tree <- mvpart::mvpart(mvpart_res ~ ., 
                            data = mvpart_covariates, 
                            control = mvpart::rpart.control(cp = 0),
                            plot.add = FALSE,
                            text.add = FALSE)
    
    if (is.null(cpTry) == TRUE){ # using cross validation if cpTry == NULL 
      
      reg_tree <- mvpart::mvpart(mvpart_res ~ ., 
                                 data = mvpart_covariates, 
                                 xv = 'min', # Selection of tree by cross-validation: "1se" - gives best tree within one SE of the overall best, "min" - the best tree, "pick" - pick the tree size interactively, "none" - no cross-validation.
                                 xval = 10, # Number of cross-validations or vector defining cross-validation groups.
                                # control = rpart.control(cp = 0),
                                 plot.add = FALSE,
                                 text.add = FALSE)
      
    } else { # prune using cpTry given
      
      reg_tree <- mvpart::prune.rpart(reg_full_tree,
                                      cp = cpTry,
                                      plot.add = FALSE,
                                      text.add = FALSE)
      
    }
    

    # rpart.plot::rpart.plot(reg_tree)
    # plotcp(reg_tree)
    piMat <- predict(reg_tree, mvpart_covariates, type="matrix")
    piMat[ piMat < 10e-10 ] <- 10e-10
    piMat[ piMat >= 1-10e-10 ] <- 1-10e-10
    # piMat[ piMat < 0.01 ] <- 0.01
    # piMat[ piMat > 0.99 ] <- 0.99
    
    
    # Incomplete log-Likelihood
    
    lic <- sum(log(rowSums(g1Mat)))
    
    # Complete data log-likelihood
    
    lc <- sum(rowSums(ziMat * log(g1Mat)))
    
    # if ( verbose == TRUE ) {
    #   print(paste('iter_cur: ', iter_cur))
    #   print(paste('cpTry: ', cpTry))
    #   print(paste('lic: ', lic))
    #   print(paste('lic in previous iteration: ', licList[iter_cur-1]))
    #   print(paste('lc: ', lc))
    #   print(paste('lc in previous iteration: ', lcList[iter_cur-1]))
    # }
    
    
    # if necessary, repeat M Step to make sure lic[cur.iter] > lic[cur.iter-1]
    
    if (lic < licList[iter_cur-1]) {
      
      # print('inside the loop')
      lc_list_from_try_cp <- rep(NA, length(cp_log10_scale_list))
      
      for (t in 1:length(cp_log10_scale_list)) {
        
        # try the different cp in the list of cp to try
        reg_tree_try <- mvpart::prune.rpart(reg_full_tree, 
                                    cp = cp_log10_scale_list[t], 
                                    plot.add = FALSE, 
                                    text.add = FALSE)
        
        # update piMat
        piMat_try <- predict(reg_tree_try, mvpart_covariates, type="matrix")
        
        # Complete data log-likelihood
        g1Mat_try <- piMat_try * gMat
        
        ziMat_try <- g1Mat_try/rowSums(g1Mat_try)
        colnames(ziMat_try) <- combVec
        lc_list_from_try_cp[t] <- sum(rowSums(ziMat_try * log(g1Mat_try)))
        
        # print(paste('t = ', t))
        # print(paste('prune_cp', cp_log10_scale_list[t]))
        # print(paste('lc = ', lc_list_from_try_cp[t]))
      }
      
      
      # find cp that maximizes the complete likelihood from the list of cps tried
      cp_prune_try <- cp_log10_scale_list[which(lc_list_from_try_cp == max(lc_list_from_try_cp))[1]]
      
      # using the cp that maximizes the lc, get new estimates for the parameters
      reg_tree_try <- mvpart::prune.rpart(reg_full_tree, cp = cp_prune_try, plot.add = FALSE, text.add = FALSE)
      
      # update piMat
      
      piMat_try <- predict(reg_tree_try, mvpart_covariates, type="matrix")
      
      # Incomplete log-Likelihood
      g1Mat_try <- piMat_try * gMat
      
      ziMat_try <- g1Mat_try/rowSums(g1Mat_try)
      colnames(ziMat_try) <- combVec
      
      lic_try <- sum(log(rowSums(g1Mat_try)))
      
      # Complete data log-likelihood
      
      lc_try <- sum(rowSums(ziMat_try * log(g1Mat_try)))
      
      
      # if ( verbose == TRUE ) {
      #   print(paste('iter_cur: ', iter_cur))
      #   print(paste('prune.cp: ', cp_prune_try))
      #   print(paste('lc: ', lc_try))
      #   print(paste('lc_in_previous_stored_iter: ', lcList[iter_cur-1]))
      #   print(paste('lic: ', lic_try))
      #   print(paste('lic_in_previous_stored_iter: ', licList[iter_cur-1]))
      # }
      
      if (lic_try < licList[iter_cur-1]){
        iter_cur <- iter_cur-1
        # print('EM converged (cannot find cp that increases incomplete data log-likelihood)')
        break 
        } else {
          cpTry <- cp_prune_try
          reg_tree <- reg_tree_try
          piMat <- piMat_try
          # piMat[ piMat < 0.01 ] <- 0.01
          # piMat[ piMat > 0.99 ] <- 0.99
          
          g1Mat <- g1Mat_try
          ziMat <- ziMat_try
          # ziMat[ ziMat < 0.01 ] <- 0.01
          # ziMat[ ziMat > 0.99 ] <- 0.99
          
          lic <- lic_try
          lc <- lc_try
        }
    }
    
    # store values
    licList[iter_cur] <- lic
    lcList[iter_cur] <- lc
  
    # stopping rules
    
    if ( abs( licList[iter_cur] - licList[iter_cur-1] ) <= 1e-4) {
      # print('EM converged based on incomplete log-likelihood')
      # print('EM converged based on incomplete log-likelihood and alpha')
      break
    }
  
    iter_cur <- iter_cur + 1
    
  }
  
  
  # one extra step to make sure final zi and pi do not have restrictions
  iter_cur <- iter_cur + 1
  
  piMat <- predict(reg_tree, mvpart_covariates, type="matrix")
  g1Mat <- piMat * gMat
  ziMat <- g1Mat/rowSums(g1Mat)
  colnames(ziMat) <- combVec
  lic <- sum(log(rowSums(g1Mat)))
  lc <- sum(rowSums(ziMat * log(g1Mat)))
  licList[iter_cur] <- lic
  lcList[iter_cur] <- lc

  if ( nrow(reg_tree$frame) > 1){
    vars_in_tree <- as.character(unique(rpart.utils::rpart.subrules.table(reg_tree)$Variable))
    vars_in_tree <- paste(vars_in_tree[order(vars_in_tree)], sep = '', collapse = ", ")
  } else{
    vars_in_tree <- 'none'
  }
  
  
  
  # marginal Z
  
  Zmarg <- matrix(NA, nrow = M, ncol = nGWAS)
  for (i in 1:nGWAS) {
    which_row <- which(binaryMat[, i] == 1)
    Zmarg[, i] = rowSums(ziMat[, which_row])
  }
  
  colnames(piMat) <- colnames(ziMat)
  colnames(Zmarg) <- paste('Zmarg_', colnames(gwasPval), sep = '')
  
  if (is.null(rownames(gwasPval))) {
    
    rownames(Zmarg) <- paste('SNP_', 1:M, sep = '')
    
  } else {
    
    rownames(Zmarg) <- rownames(gwasPval)
    
  }
  # list of variables to output
  message('Info: Incomplete log-likelihood converged in Stage 2 of multiGPATree.')
  message( "Info: Stage 2 of multiGPATree completed.")
  
  EM_out <- list()
  EM_out$numIterConvergence <- iter_cur
  EM_out$fit <- reg_tree
  EM_out$fitSelectVar <- vars_in_tree
  EM_out$Z <- ziMat
  EM_out$Zmarg <- Zmarg
  EM_out$pi <- piMat
  EM_out$licVec <- licList[1:iter_cur]
  EM_out$lcVec <- lcList[1:iter_cur]
  EM_out$annMat <- annMat
  EM_out$gwasPval <- gwasPval
                 
  return(EM_out)
  
}


# # Example test
# 
# alpha <- c(0.2, 0.2)
# test <- GPATree::aim2_sim_data (M = 500,
#                                nGWAS = 2,
#                                nAnn = 10,
#                                percent_ones_ann = 0.20,
#                                percent_overlap_ann = 0.5,
#                                percent_overlap_gwas = 0.75,
#                                trueAlphaVec = alpha)
# # checkheat <- as.matrix(cbind(test$zMat, test$annMat))
# # heatmap(checkheat, Colv = NA, Rowv = NA, scale = 'column')
# colSums(test$zMat)
# par(mfrow=c(1,2))
# for (i in 1:length(alpha)) {
#   hist(test$gwasPval[, i], breaks = 100, xlab = paste('P', i, sep = ''),
#        main = paste('Beta(', alpha[i], ', 1)', sep = ''))
# }
# 
# par(mfrow=c(1,1))
# testres <- aim2_EM(gwasPval = test$gwasPval,
#                    annMat = test$annMat,
#                    initAlpha = 0.1,
#                    cpTry = 0.001)
# testres$fit
# rpart.plot::rpart.plot(testres$fit, extra = 103)
# colSums(test$Zmat)
