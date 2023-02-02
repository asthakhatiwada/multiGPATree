#' multiGPATree selected decision tree
#'
#' This function will extract the combinations of functional annotations selected by the multivariate decision tree in the stage 2 of multiGPATree method.
#'
#' @author  Aastha Khatiwada
#' @param object An object of class multiGPATree.
#' 
#' @return A list containing variables in combinations selected by the multivariate decision tree (CART_PIs), the combination (CART_PIs_comb) and the predicted proportions for the selected PIs(assoc_pred).
#' @import mvpart
#' @importFrom dplyr mutate
#' @importFrom dplyr %>%
#' @name decTree
#' @aliases decTree,multiGPATree-method
setMethod(
  f="decTree",
  signature="multiGPATree",
  definition=function( object ) {
  
  
  # start here
  
  CARTmod <- object@fit$fit
  minPredictedProb <- 0
  
  # if (is.null(minPredictedProb)){
  #   minPredictedProb <- 0
  # }
  
  if ( nrow(CARTmod$frame) == 1){
    
    if (ncol(CARTmod$frame) >= 8 ){
      yvals <- CARTmod$frame$yval
      ans <- list(CART_PIs = 'none',
                  CART_PIs_comb = 'none',
                  assoc_pred = CARTmod$frame$yval)
    }
    
    
    
  }
  
  if ( nrow(CARTmod$frame) > 1){
    
    term_nodes <- as.numeric( rownames(CARTmod$frame)[c( which(CARTmod$frame[, 1] == "<leaf>") )  ] )
    pth <- mvpart::path.rpart(tree = CARTmod, nodes = term_nodes, pretty = 0, print.it = F) #list of the PIs for CART
    
    if (ncol(CARTmod$frame) >= 8 ){
      yvals <- (CARTmod$frame$yval)[c(which(CARTmod$frame[, 1] == "<leaf>"))]
      y_ids <- 1:length(yvals)
      # y_ids <- order(yvals, decreasing = T) # [1:2] # selecting the cases with top two predicted y
      yvals <- yvals[c(y_ids)]
      num_CARTpis <- length(yvals) #number of PIs identified by CART
      CART_PIs <- vector("list", num_CARTpis)
      
      for (i in 1:num_CARTpis)  {
        curr_PI <- setdiff(pth[[y_ids[i]]], "root") #PI in CART format
        # curr_PI <- curr_PI[order(curr_PI)]
        elems <- length(curr_PI) #number elements in current PI
        CPM_nm <- c()
        
        for (j in 1:elems) {#NOTE 1st element always =n root
          parts <- strsplit(curr_PI[j], "=")
          part_which <- parts[[1]][2]
          
          if (part_which == '0') {#Note these variables are compliments
            el <- strsplit(curr_PI[j], "=")[[1]][1]
            elem <- paste("!", el, sep = "")
            CPM_nm <- append(CPM_nm, elem)
          }
          
          if (part_which == '1') {
            elem <- strsplit(curr_PI[j], "=")[[1]][1]
            CPM_nm <- append(CPM_nm, elem)
          }
        }
        
        CART_PIs[i] <- paste('(',  paste(CPM_nm, collapse = ' & ', sep = ''), ')', sep = '')
        
      }
      
      paste(unlist(CART_PIs), collapse = ' OR ')
      
      
      # combining only selected combinations
      yvals_select <- which(yvals > minPredictedProb)
      
      if (length(yvals_select) > 0){
        comb_select <- paste(unlist(CART_PIs[yvals_select]), collapse = ' OR ')
      } else {
        comb_select <-  'none'
      }
      
      ans <- list(CART_PIs = CART_PIs[yvals > minPredictedProb],
                  CART_PIs_comb = comb_select,
                  assoc_pred = yvals[yvals > minPredictedProb])
    }
  }
  return(ans)
  }
)
  


