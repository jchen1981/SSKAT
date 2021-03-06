% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SSKAT.R
\name{mSKAT2}
\alias{mSKAT2}
\title{Small-sample adjusted SKAT for a multivariate (continuous) outcome ('m' stands for 'multivariate') without de-correlation. 
De-correlation could be problematic when the number of outcome variables is not small compared to the sample size,
due to the large estimation error for the covariance matrix. In such case, de-correlation is not recommended unless
the correlation is very high.}
\usage{
mSKAT2(
  formula.H0 = NULL,
  data = NULL,
  Y,
  K,
  acc = 1e-05,
  lim = 10000,
  tol = 1e-10
)
}
\arguments{
\item{formula.H0}{a one-sided linear formula object under the null, indicating the variables to adjust.}

\item{data}{a data frame containing the variables named in formula. It is required if the formula is not null or the intercept model (~ 1).}

\item{Y}{a matrix for the multivariate phenotypes, row - samples, col - phenotypes.}

\item{K}{the kernel matrix, which quantifies the similarities between samples.}

\item{acc}{error bound for \code{davies} method.}

\item{lim}{maximum number of integration terms for \code{davies} method.}

\item{tol}{the eigenvalue cutoff, below which is considered to be 0. This is used to reduce the computation burden.}
}
\value{
\describe{A list containing
  \item{p.value}{ association p-value}
  \item{Q.adj}{ adjusted score statistic}
}
}
\description{
Compute the adjusted score statistic and p-value for a multivariate continuous outcome without de-correlation.
}
\examples{
set.seed(123)
L <- matrix(rnorm(1000), 100, 10)  # Latent factor to induce correlation in Y
Y <- scale(L \%*\% matrix(rnorm(10 * 40), 10, 40)) + matrix(rnorm(100 * 40), 100, 40)
G <- scale(Y \%*\% matrix(rnorm(40 * 10), 40, 10)) * 0.2 + matrix(rnorm(100 * 10), 100, 10)
K <- G \%*\% t(G)
mSKAT2(Y = Y, K = K)  # No de-correlation
mSKAT(Y = Y, K = K)   # De-correlation         
}
\references{
Zhan X, et al. (2017)  A small-sample multivariate kernel machine test for microbiome association studies.
Genet Epidemiol. 41(3):210-220.
}
\author{
Jun Chen
}
\keyword{SKAT}
