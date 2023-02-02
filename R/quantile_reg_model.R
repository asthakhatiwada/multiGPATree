#' quantile_reg_model result
#'
#' Quantile regression model to predict complexity parameter in Stage 2 of the GPA-Tree method
#'
#' @format An object of class 'rq'. Contains 14 elements:
#' \describe{
#' \itemize{
#' \item{coefficients: } {coefficient of the quantile regression model}
#' \item{x: } {provides the x side of the regression model}
#' \item{y: } {provides the y side of the regression model}
#' \item{residuals: } { the residuals from the fit.}
#' \item{dual: }{the vector dual variables from the fit}
#' \item{fitted.values: }{}
#' \item{formula: }{ formula used to fit the quantile regression model.}
#' \item{terms: }{terms of the model}
#' \item{xlevels: }{}
#' \item{call: }{ function call}
#' \item{tau: }{ percentile used in the quantile regression}
#' \item{rho: }{ The value(s) of objective function at the solution.}
#' \item{method: }{the algorithmic method used to compute the fit. There are several options: The default method is the modified version of the Barrodale and Roberts algorithm for l1-regression, used by l1fit in S, and is described in detail in Koenker and d'Orey(1987, 1994), default = "br". This is quite efficient for problems up to several thousand observations, and may be used to compute the full quantile regression process. It also implements a scheme for computing confidence intervals for the estimated parameters, based on inversion of a rank test described in Koenker(1994). For larger problems it is advantageous to use the Frisch–Newton interior point method "fn". And for very large problems one can use the Frisch–Newton approach after preprocessing "pfn". Both of the latter methods are described in detail in Portnoy and Koenker(1997), this method is primarily well-suited for large n, small p problems where the parametric dimension of the model is modest. For large problems with large parametric dimension it is usually advantageous to use method "sfn" which uses the Frisch-Newton algorithm, but exploits sparse algebra to compute iterates. This is typically helpful when the model includes factor variables that, when expanded, generate design matrices that are very sparse. A sixth option "fnc" that enables the user to specify linear inequality constraints on the fitted coefficients; in this case one needs to specify the matrix R and the vector r representing the constraints in the form Rb ≥q r. See the examples. Finally, there are two penalized methods: "lasso" and "scad" that implement the lasso penalty and Fan and Li's smoothly clipped absolute deviation penalty, respectively. These methods should probably be regarded as experimental.}
#' \item{model: }{ optionally the model frame, if model=TRUE.}
#' }
#' }
#' @source Quantile Regression Model fitted using simulated data
#'
#' @seealso\link[quantreg]{rq}
"quantile_reg_model"
