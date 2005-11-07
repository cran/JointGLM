\name{Joint Modelling}
\alias{summplot.glm}
\alias{obs.vs.model}
\alias{rstand.vs.linpred}
\alias{res.vs.explvar}
\alias{absres.vs.fitted}
\alias{adjvar.vs.linpred}
\alias{qqglm}
\title{Model Checking}
\description{Graphicals representation to help in model checking.}

\usage{summplot.glm(glm, var = NULL, res = 'standard',
which = 1:5, ask = nb.fig < length(which) &&
 dev.interactive(), smooth = TRUE, \dots)
obs.vs.model(glm.joint, plot.disp = FALSE, \dots)
rstand.vs.linpred(glm, smooth = TRUE, \dots)
res.vs.explvar(glm, var, res = 'standard', smooth = TRUE, \dots)
absres.vs.fitted(glm, res = 'standard', smooth = TRUE, \dots)
adjvar.vs.linpred(glm, smooth = TRUE, \dots)
qqglm(glm, \dots)
}
\arguments{
  \item{glm}{Object of class \emph{glm}.}
  \item{glm.joint}{A list with two component. Each one are
    \emph{glm} object corresponding to the mean and dispersion. The
    dispersion object is necessary only if plot.dispt is set to
    \code{TRUE}.}
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
}
\value{Return a (several) graphic window(s).}
\details{The function \code{\link{summplot.glm}} is a special function
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
ajust <- joint.glm('Y~X', 'd~X', data = data.frm)
layout(matrix(c(1,1,2,2,3,3,4,4,0,5,5,0),3,byrow=TRUE))
summplot.glm(ajust$glm.mean, var = 'X')
}
\references{
  McCullagh P. and Nelder J. A. (1987). \emph{Generalized Linear Models
    Second Edition.} London: Chapman and Hall. ISBN 0-412-31760-5.}
\keyword{hplot}
\author{Mathieu Ribatet and Bertrand Iooss}