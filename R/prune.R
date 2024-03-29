#' Prune multiGPATree model fit
#'
#' This function will prune the multiGPATree model fit using the given cp value.
#'
#' @author  Aastha Khatiwada
#'
#' @param object An object of class multiGPATree.
#' @param cp The cp parameter to be used for pruning. cp must be between 0 and 1.
#'
#' @return multiGPATree model output.
#' @examples
#' \dontrun{
#' library(multiGPATree)
#'
#' # load multiGPATree example data
#' data(simdata)
#'
#' #fitting the multiGPATree model
#' fit <- multiGPATree(simdata$gwasPval, simdata$annMat)
#'
#' # pruning the multiGPATree model fit
#' pruned.fit <- prune(fit, cp = 0.005)
#' }
#' @importFrom mvpart prune prune.rpart
#' @export
#' @name prune
#' @aliases prune,multiGPATree-method
setMethod(
  f="prune",
  signature="multiGPATree",
  definition=function(object, cp=0.001) {
    
    # cp = 0.245
    # object <- res@fit$P1_P2
    
  # start here
    
  
  if ( cp < 0 | cp >1 ) {
    stop( "Inappropriate value for 'cp' argument. It should be between zero and one." )
  }
  
  
  pruned_tree <- object@fit
  # pruned_tree$fit
  
  if (ncol(object@fit$Zmarg) == 1) {
    pruned_tree$fit <- mvpart::prune(object@fit$fit, cp = cp)
  }
  
  
  if (ncol(object@fit$Zmarg) > 1) {
    pruned_tree$fit <- mvpart::prune.rpart(object@fit$fit, 
                                           cp = cp, 
                                           plot.add = FALSE, 
                                           text.add = FALSE)
  }
  
  pruned_tree$fitSelectVar <- pruned_tree$fit$frame$var
  pruned_tree$fitSelectVar <- pruned_tree$fitSelectVar[pruned_tree$fitSelectVar != '<leaf>']
  pruned_tree$fitSelectVar <- ifelse(length(pruned_tree$fitSelectVar) > 1,
                                         paste(as.character(sort(unique(pruned_tree$fitSelectVar))), 
                                               collapse = ', '),
                                         ifelse(length(pruned_tree$fitSelectVar) == 1,
                                                as.character(unique(pruned_tree$fitSelectVar)),
                                                'none')
                                     )
  
  # pruned_tree$fitSelectVar
  
  message( "Info: multiGPATree model pruned based on the complexity parameter (cp) value of ", round(cp, 6), "." )
  return_list <- pruned_tree
  
  new( Class = "multiGPATree",
       fit = return_list)
  
  }

  
)

 
  
  
