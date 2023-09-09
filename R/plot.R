#' Plot the functional annotation tree obtained from the multiGPATree model fit.
#'
#' This function will plot the functional annotation tree for the multiGPATree model fit.
#'
#' @author  Aastha Khatiwada
#'
#' @param x An object of class multiGPATree.
#' @param y missing (not required).
#' @param ... ...
#'
#' @return Returns a plot for the functional annotation tree from the multiGPATree model fit.
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
#' # plotting the multiGPATree model fit
#' plot(fit)
#' }
#' @name plot
#' @aliases plot,multiGPATree,missing-method
#' @export
setMethod(
  f="plot",
  signature=c(x = "multiGPATree", y = "missing"),
  definition=function( x, y, ... ) {
    
    # x <- res@fit$P1_P2
    
    #start here
    out <- x@fit$fit
    leaves_index <- which(out$frame$var=='<leaf>')
    out$frame$leaf <- rep(NA, nrow(out$frame))
    out$frame$leaf[leaves_index] <- paste('LEAF ', 1:length(leaves_index), sep = '')
    localFDR <- vector("list", ncol(x@fit$Zmarg))
    
    
    for (j in 1:length(localFDR)) {
      localFDR[[j]] <- rep(NA, nrow(out$frame))
      for (i in 1:length(leaves_index)) {
        SNPindex <- which(out$where == leaves_index[i])
        localFDR[[j]][leaves_index[i]] <- round(mean(1 - x@fit$Zmarg[SNPindex, j], na.rm = TRUE), 4)
      }
      out$frame$v1 <- localFDR[[j]]
      colnames(out$frame)[ncol(out$frame)] <- paste('local FDR P', j, sep = '')
    }
   
    
    node.fun1 <- function(x, labs, digits, varlen) {
    
      if (ncol(out$frame) > 10) {
        p1 <- paste(out$frame$leaf, "\n \n N = ", out$frame$n, sep = '')
        n_GWAS  <- log(ncol(out$y))/log(2)
        for (i in 1:n_GWAS) {
          p1 <- paste(p1, 
                      "\n \n local FDR P", i, " = ", 
                      formatC(round(out$frame[, 10+i], 3), 3, format = 'f'),
                      sep = '')
        }
      }
      
    
    if (ncol(out$frame) == 10) {
      colnames(out$frame)[ncol(out$frame)] <- 'local FDR'
      p1 <- paste(out$frame$leaf, "\n \n N = ", out$frame$n, 
                  "\n \n local FDR = ", formatC(round(out$frame$`local FDR P1`, 3), 3, format = 'f'),
                  sep = '')
      }
      return(p1)
    }
    rpart.plot::rpart.plot(out, 
                           node.fun = node.fun1,
                           type = 5, # type = 5, 3 works
                           extra = 1,
                           font = 2,
                           split.font = 2,
                           col = 'black',
                           under.col=10,
                           compress = TRUE,
                           ycompress = TRUE,
                           roundint=FALSE) #roundint = FALSE silences the warning message
  }
)


