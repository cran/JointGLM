\name{eql}
\alias{eql}
\title{Compute the Extended Quasi-Likelihood and the Extended
  Quasi-Deviance}
\description{Compute the Extended Quasi-Likelihood and the Extended
  Quasi-Deviance for a Joint Model.}

\usage{eql(mod.mean, mod.disp, df.adj = TRUE)}
\arguments{
  \item{mod.mean, mod.disp}{Two fitted GLMs/GAMs correponding to
    the mean and the dispersion components respectively.}
  \item{df.adj}{Logical. If \emph{TRUE} - the default, adjustment on
    degree of freedom is applied for the EQD.}
}
\value{
  This function returns the Extended Quasi-Likelihood and the Extended
  Quasi-Deviance for a joint model.}
\details{
  The Extended Quasi-Likelihood (EQL) generalizes the Quasi-Likelihood
  as it behaves as a likelihood function on both mean and dispersion
  parameters. The Extended Quasi-Deviance is defined as
  \deqn{\sum_i \frac{d_i}{\phi_i} + \sum_i \log(2\pi \phi_i
  V(y_i))}{\sum_i d_i + \sum_i phi_i log( 2\pi phi_i V(y_i) )}
  where \eqn{d_i} is the deviance contribution of the i-th
  observation, \eqn{\phi_i} is the fitted values of dispersion
  parameter for the i-th observation, \eqn{V} is the variance
  function of the mean GLM and \eqn{y_i} is the i-th observation.
}
\author{Mathieu Ribatet and Bertrand Iooss}
\references{
  McCullagh P. and Nelder J. A. (1987). \emph{Generalized Linear Models
    Second Edition.} London: Chapman and Hall. ISBN 0-412-31760-5.}
\keyword{htest}
