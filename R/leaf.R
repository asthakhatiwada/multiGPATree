#' Functional annotation tree returned by the multiGPATree model.
#'
#' This function will provide the annotation combinations relevant to risk-associated SNPs.
#'
#' @author  Aastha Khatiwada
#'
#' @param object An object of class multiGPATree.
#'
#' @return Returns a matrix where each row corresponds to a leaf from the multiGPATree model fit and contains information regarding the local FDR for SNPs in the leaf, and also information regarding annotations that are enriched (1) or not enriched (0) for the leaf.
#' @examples
#' \dontrun{
#' library(multiGPATree)
#'
#' # load multiGPATree example data
#' data(simdata)
#'
#' #fitting the multiGPATree model
#' fit <- multiGPATree(simdata$gwasPval, simdata$annMat)
#' leaf(fit)
#' }
#' @export
#' @name leaf
#' @aliases leaf,multiGPATree-method
setMethod(
  f="leaf",
  signature="multiGPATree",
  definition=function( object ) {
  
  # start here
    options(digits = 4)
  CARTmod <- object@fit$fit
  
  if ( nrow(CARTmod$frame) == 1){
    if (ncol(CARTmod$frame) >= 8 ){
      mean_localFDR <- colMeans(1-object@fit$Zmarg, na.rm = TRUE)
      select_annMat <- as.data.frame(matrix(c(mean_localFDR, "No annotations selected"), 
                                            nrow = 1, 
                                            ncol = length(mean_localFDR) + 1)
                                     )
      rownames(select_annMat) <- 'LEAF 1'
      colnames(select_annMat) <- c(paste('local FDR', names(mean_localFDR)), "Note")
      # select_annMat$`local FDR` <- as.numeric(select_annMat$`local FDR`)
      for (i in 1:(ncol(select_annMat)-1)) {
        select_annMat[, i] <- as.numeric(select_annMat[, i])
      }
    }
  }
    
  if ( nrow(CARTmod$frame) > 1){
    
    term_nodes <- as.numeric( rownames(CARTmod$frame)[c( which(CARTmod$frame[, 1] == "<leaf>") )  ] )
    pth <- mvpart::path.rpart(tree = CARTmod, nodes = term_nodes, pretty = 0, print.it = F) #list of the PIs for CART
    anns <- setdiff(unlist(pth), 'root')
    parts <- strsplit(anns, "=")
    select_ann <- rep(NA, length(parts))
    
    for (i in 1:length(parts)) { select_ann[i] <- parts[[i]][1]    }
    
    select_ann <- unique(select_ann)
    select_annMat <- matrix(NA, nrow = length(pth), ncol = length(select_ann)+ncol(object@fit$Zmarg))
    
    if (ncol(object@fit$Zmarg) == 1){
      colnames(select_annMat) <- c('local FDR', select_ann)
    }
    
    if (ncol(object@fit$Zmarg) > 1){
      colnames(select_annMat) <- c(paste('local FDR P', 
                                         1:ncol(object@fit$Zmarg), sep = ''), 
                                   select_ann)
    }
    rownames(select_annMat) <- paste('LEAF ', 1:length(pth), sep = '')
    cartmod_frame_ind <- sort(unique(CARTmod$where))
    
    if (ncol(CARTmod$frame) >= 8 ){
      
      for (i in 1:length(pth))  {
        select_annMat[i, 1:ncol(object@fit$Zmarg)] <- colMeans(1-as.matrix(object@fit$Zmarg[which(CARTmod$where == cartmod_frame_ind[i]), ]), na.rm = TRUE)
        
        curr_PI <- setdiff(pth[[i]], "root") #PI in CART format
        
        for (j in 1:length(curr_PI)) {#NOTE 1st element always =n root
          parts <- strsplit(curr_PI[j], "=")
          select_annMat[i, which(colnames(select_annMat) == parts[[1]][1])] <- as.numeric(parts[[1]][2])
        }
      }
    }
  }
  
  select_annMat <- as.data.frame(select_annMat)
  select_annMat[is.na(select_annMat)] <- '-'
  return(select_annMat)
  }
)

