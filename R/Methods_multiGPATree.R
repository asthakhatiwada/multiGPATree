#' @import graphics
#' @import methods
#' @importFrom rpart.plot rpart.plot
# generic methods for "multiGPATree" class
setMethod(
  f="show",
  signature="multiGPATree",
  definition=function( object ) {

    options(digits = 4)
    # object <- fit
    
    # constants

    M <- nrow(object@gwasPval)
    nAnn <- ncol(object@annMat)
    nGWAS<- ncol(object@gwasPval)
    l <- leaf(object)
    # estimates

    est_alpha = object@fit$alpha
    cat( "Summary: multiGPATree model results (class: multiGPATree)\n" )
    cat( "--------------------------------------------------\n" )
    cat( "Data summary:\n" )
    cat( "\tNumber of GWAS data: ", nGWAS , "\n", sep="" )
    cat( "\tNumber of Annotations: ", nAnn , "\n", sep="" )
    cat( "\tNumber of SNPs: ", M , "\n", sep="" )
    if(nAnn!=0) {
    cat( "\tAlpha estimates: ", paste(as.character(round(est_alpha, 4)), collapse=", ") , "\n", sep="" )
    cat( "Functional annotation tree description:\n" )
    if (nrow(l) <= 15){
      print( l )
    } else {
      cat("\tNumber of leaves (terminal nodes): ", nrow(l), "\n", sep="")
    }
    } else {}
    cat( "--------------------------------------------------\n" )
  }
)

