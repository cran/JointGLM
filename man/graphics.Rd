\name{plot.joint}
\alias{plot.joint}
\alias{obs.vs.model}
\alias{rstand.vs.linpred}
\alias{res.vs.explvar}
\alias{absres.vs.fitted}
\alias{adjvar.vs.linpred}
\alias{qqglm}
\title{Model Checking}
\description{Graphicals representation to help in model checking.}

\usage{
\method{plot}{joint}(x, comp = "mean", var = NULL, res = 'standard',
which = 1:5, ask = nb.fig < length(which) && dev.interactive(), smooth =
TRUE, \dots)
obs.vs.model(x.joint, plot.disp = FALSE, \dots)
rstand.vs.linpred(x, smooth = TRUE, \dots)
res.vs.explvar(x, var, res = 'standard', smooth = TRUE, \dots)
absres.vs.fitted(x, res = 'standard', smooth = TRUE, \dots)
adjvar.vs.linpred(x, smooth = TRUE, \dots)
qqglm(x, \dots)
}
\arguments{
  \item{x,x.joint}{A fitted GLM/GAM object or a fitted joint model.}
  \item{comp}{The component for the joint model. Must be one of
    \code{"mean"} or \code{"disp"}.}
  \item{var}{Character. The name of the explanatory variable selected.}
  \item{res}{Character. Should be \emph{standard}, \emph{student} of
    \emph{brut}. The residual type considered.}
  \item{which}{Numeric vector. Which plot must be produced.}
  \item{ask}{Logical. If \code{TRUE}, user is asked before a new plot is
    produced.}
  \item{smooth}{Logical. If \code{TRUE} - the default, a smoothing curve
    is plotted thanks to the \code{\link{lowess}} fucntion.}
  \item{\dots}{Optional parameters to be passed to \code{\link{plot}},
    \code{\link{abline}} and \code{\link{qqnorm}} functions.}
  \item{plot.disp}{Logical. Should the +/- standard deviation error bar
    be plottted around the point selected by the \code{\link{identify}}
 function.} 
}
\value{Return a (several) graphic window(s).}
\details{The function \code{\link{plot.joint}} is a special function
  to compute all the checking plots. Checking plots implemented are :
  \enumerate{
    \item Observations functions of Fitted values
    \item Standardized Residuals functions of Linear predictor
    \item Residuals functions of an explanatory variable
    \item Absolute Residuals functions of Fitted values
    \item Adjusted Dependent Variable functions of Linear predictor
    \item \emph{QQ-plot} of the residuals
    }
  Several conclusions can be established thanks to these
  plots. One can refer to \cite{McCullagh and Nelder, 1987} for
  interpretation.}
\examples{
X <- rnorm(500, 5)
mu <- 2*X + 1
sigma2 <- exp(-X)
Y <- rnorm(500, mu, sqrt(sigma2))
data.frm <- data.frame(Y=Y, X=X)
ajust <- fitjoint("glm", 'Y~X', 'd~X', data = data.frm)
layout(matrix(c(1,1,2,2,3,3,4,4,0,5,5,0),3,byrow=TRUE))
plot(ajust, var = 'X')
}
\references{
  McCullagh P. and Nelder J. A. (1987). \emph{Generalized Linear Models
    Second Edition.} London: Chapman and Hall. ISBN 0-412-31760-5.}
\keyword{hplot}
\author{Mathieu Ribatet and Bertrand Iooss}
