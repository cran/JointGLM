\name{Joint Modelling}
\alias{hat.glm}
\title{Hat values for the Joint model}
\description{Computes hat values for the joint model.}

\usage{hat.glm(glm)}
\arguments{
  \item{glm}{A object of the class \emph{glm}.}
}
\value{This function returns the \emph{hat values} i.e. the
  diagonal of the \emph{hat matrix}.
}
\details{
  The expression of the genralized form for the \emph{hat matrix} is
  given by
  \deqn{H = W^{1/2}X \left(X^TWX\right)^{-1}X^TW^{1/2}}{ W^(1/2) X (
  t(X) W X )^(-1) t(X) W^(1/2)}
  where \eqn{W} is diagonal matrix representing the weigths,
  \eqn{X} a matrix corresponding to explicative variables
  realisations.
}
\keyword{htest}
\author{Mathieu Ribatet and Bertrand Iooss}
