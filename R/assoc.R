#' Association mapping
#'
#' This function will implement association mapping for the multiGPATree model.
#'
#' @author  Aastha Khatiwada
#'
#' @param object An object of class list
#' @param FDR FDR level. Value has to be between 0 and 1.
#' @param fdrControl Method to control FDR. Possible values are "global" (global FDR control) and "local" (local FDR control).
#'
#' @return Returns a MX4 matrix where the row represents SNPs, the fist column indicates the association between each SNP and phenotype 1, and second column indicates the leaf in which the SNP falls
#' @examples
#' \dontrun{
#' library(multiGPATree)
#'
#' # load multiGPATree example data
#' data(simdata)
#'
#' #fitting the GPATree model
#' fit <- multiGPATree(simdata$gwasPval, simdata$annMat)
#'
#' # association mapping for the multiGPATree model fit
#' assoc.fit <- assoc(fit, FDR = 0.05, fdrControl = "global")
#' }
#' @export
#' @name assoc
#' @aliases assoc,multiGPATree-method

setMethod(
  f="assoc",
  signature="multiGPATree",
  definition=function( object, FDR=0.05, fdrControl="global") {
    
    # object <- res2@fit$P2_P3
    # class(object)
    
    # start here
  
    # check arguments
    
    if ( fdrControl != "global" & fdrControl != "local" ) {
      stop( "Invalid value for 'fdrControl' argument! It should be either 'global' or 'local'." )
    }
    
    if ( FDR < 0 | FDR > 1 ) {
      stop( "Invalid value for 'FDR' argument! FDR should be between zero and one." )
    }
    
    # association mapping
    
    
    if (ncol(object@fit$Zmarg) == 1) {
      fdrmat <- 1 - object@fit$Zmarg # for marginal P1 add zmat columns 00 and 01; for P2 add zmat columns 00 and 10
      amat <- matrix( 0, nrow(fdrmat), ncol(fdrmat))
      colnames(amat) <- paste('P', 1:ncol(fdrmat), sep = '')
    }
    
    if (ncol(object@fit$Zmarg) > 1) {
      fdrmarg <- 1 - object@fit$Zmarg # for marginal P1 add zmat columns 00 and 01; for P2 add zmat columns 00 and 10
      fdrjoint <- 1 - object@fit$Z[, 4] # for joint P1 and P2, add Zmat columns 00, 10, 01
      fdrmat <- cbind(fdrmarg, fdrjoint)
      colnames(fdrmat)[3] <- paste(colnames(fdrmarg), collapse = '_')
      amat <- matrix( 0, nrow(fdrmat), ncol(fdrmat))
      colnames(amat) <- colnames(fdrmat)
    }
    
   rownames(amat) <- rownames(object@fit$Zmarg)
    
    if ( fdrControl == "local" ) {
      message( "Info: Association mapping based on the local FDR control at level ", FDR, "." )
      amat[ fdrmat <= FDR ] <- 1
      } else if ( fdrControl == "global" ) {
        message( "Info: Association mapping based on the global FDR control at level ", FDR, "." )
        # direct approach for FDR control
        for ( j in 1:ncol(amat) ) {
          pp <- fdrmat[,j]
          pp.ordered <- sort(pp)
          pp.cum <- cumsum( pp.ordered ) / c(1:length(pp))
          cutoff <- max( pp.ordered[ pp.cum <= FDR ], 0 )
          amat[ pp <= cutoff, j ] <- 1
        }
      }
   # leaf mapping ####
    
    if (ncol(amat) >= 1){
      amat <- as.data.frame(amat)
      leaf_frame <- object@fit$fit$frame[object@fit$fit$frame$var == '<leaf>', c('var', 'n')]
      leaf_frame$var <- paste('LEAF', 1:nrow(leaf_frame))
      tree_rules <- decTree(object = object)
      leaf_frame$rule <- rep(NA, nrow(leaf_frame))
      
      for (i in 1:length(tree_rules$CART_PIs)) {  
        leaf_frame$rule[i] <- tree_rules$CART_PIs[[i]]
      }
      whichnodes <- rownames(object@fit$fit$frame)[object@fit$fit$where]
      amat$leaf <- leaf_frame$var[match(whichnodes, rownames(leaf_frame))]
    }
   return(amat)
  }
  )






