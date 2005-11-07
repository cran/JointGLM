\name{Joint Modelling}
\alias{joint.glm}
\title{Joint Modelling Fitting}
\description{Fit two interlinked GLM's to model dependently
  the Mean and the Dispersion.}

\usage{joint.glm(form.mean, form.disp, data,
family.mean = gaussian, family.disp = Gamma(link='log'),
eps = 10^(-6), iter.max = 100, maxit = 50, reml = TRUE)}
\arguments{
  \item{form.mean, form.disp}{R \emph{formula} put in a character
    type way which specifies the linear model for the mean
    and the dispersion respectively.}
  \item{data}{A \emph{data.frame} for the data. The first column
    corresponds to the response observations, while others to the
    explicative variables.}
  \item{family.mean, family.disp}{The \emph{family} for the mean and
    dispersion. See \code{\link{family}}.}
  \item{eps}{Optional numeric. The precision for the convergence
    test.}
  \item{iter.max, maxit}{Numerics. The maximum number of iteration in
    the Joint Modelling algorithm and \code{\link{glm.fit}} function
    respectively.}
  \item{reml}{Logical. If \emph{TRUE} - the default, the Restricted
    Maximum Likelihood Estimation is used.}
}
\value{This function returns a list with four components. The first two
  are \emph{glm ocbjects} corresponding to the final estimate of the
  mean and dispersion GLM's component. A component \emph{iterations}
  corresponding to the number of iterations in the Joint Modelling
  fitting algorithm. And a component \emph{eqd} corresponding to the
  Extended Quasi-Deviance for this Joint Model Object.}
\details{The estimation procedure of the Joint Model is based on the
  \emph{see-saw} algorithm. This procedure is explained in
  \cite{McCullagh and Nelder, 1987}. The
  pair of models is fitting thanks to the Extended Quasi-Deviance
  criterion. We start by setting the dispersion parameter for the mean
  to 1 and fitting the model for the mean. We then form \emph{d} the
  response variable for dispersion and fit the dispersion GLM. From the
  fitted values of the dispersion GLM, we form prior weights given by
  the inverse of these fitted values and refit the model for the
  mean. The alternation of fits continues until stability is
  reached. We say that stability is reached when the EQD variation - in
  percent - between two iteration in lower than an \eqn{\epsilon}.}
\examples{
X <- rnorm(500, 5)
mu <- 2*X + 1
sigma2 <- exp(-X)
Y <- rnorm(500, mu, sqrt(sigma2))
data.frm <- data.frame(Y=Y, X=X)
ajust <- joint.glm('Y~X', 'd~X', data = data.frm)
summary(ajust$glm.mean)
summary(ajust$glm.disp)
}
\references{
  Lee, Y and Nelder, J. A. (2003). \emph{Robust Design via Generalized
    Linear Models.} Journal of Quality Technology, 35:2-12.
  McCullagh P. and Nelder J. A. (1987). \emph{Generalized Linear Models
    Second Edition.} London: Chapman and Hall. ISBN 0-412-31760-5.}
\keyword{htest}
\author{Mathieu Ribatet and Bertrand Iooss}

