\name{symmetric.KL}
\alias{symmetric.KL}

\title{
Symmetrized Kullback-Leibler divergence
}
\description{
Compute the Symmetrized Kullback-Leibler divergence between a pair of normally distributed clusters.
}
\usage{
symmetric.KL(mean1, mean2, cov1, cov2)
}

\arguments{
  \item{mean1}{mean vector of length \code{p} for cluster 1, where \code{p} is the dimension of the clusters.}
  \item{mean2}{mean vector of length \code{p} for cluster 2.}
  \item{cov1}{ \code{pxp} covariance matrix for cluster 1.}
  \item{cov2}{\code{pxp} covariance matrix for cluster 2.}
}
\details{
Consider two \code{p}-dimensional, normally distributed clusters with centers \eqn{\mu1, \mu2} and covariance matrices \eqn{\Sigma1, \Sigma2}. We compute the KL divergence \code{d12} between the clusters as follows:
\deqn{d12 = 1/4 * ( t(\mu2 - \mu1) * ( \Sigma1^(-1) + \Sigma2^(-1) ) * (\mu2 - \mu1) +
    trace(\Sigma1/\Sigma2 + \Sigma2/\Sigma1) + 2p ) }
The dimension of the clusters must be same.

Note that KL-divergence is not symmetric in its original form. We converted it symmetric by averaging both way KL divergence. The symmetrized KL-divergence is not a metric because it does not satisfy triangle inequality.
}
\value{
\code{symmetric.KL} returns a numeric value measuring the Symmetrized Kullback-Leibler divergence between a pair of normally distributed clusters. 
}
\references{
Abou--Moustafa, Karim T and De La Torre, Fernando and Ferrie, Frank P (2010) Designing a Metric for the Difference between Gaussian Densities; Brain, Body and Machine, 57--70. 
}
\author{Ariful Azad}


\seealso{
 \code{\link{mahalanobis.dist}}, \code{\link{dist.cluster}}
}
\examples{
## ------------------------------------------------
## load data and retrieve a sample
## ------------------------------------------------

library(healthyFlowData)
data(hd)
sample = exprs(hd.flowSet[[1]])


## ------------------------------------------------
## cluster sample using kmeans algorithm
## ------------------------------------------------

km = kmeans(sample, centers=4, nstart=20)
cluster.labels = km$cluster

## ------------------------------------------------
## Create ClusteredSample object  
## and compute mahalanobis distance between two clsuters
## ------------------------------------------------

clustSample = ClusteredSample(labels=cluster.labels, sample=sample)
mean1 = get.center(get.clusters(clustSample)[[1]])
mean2 = get.center(get.clusters(clustSample)[[2]])
cov1 = get.cov(get.clusters(clustSample)[[1]])
cov2 = get.cov(get.clusters(clustSample)[[2]])
n1 = get.size(get.clusters(clustSample)[[1]])
n2 = get.size(get.clusters(clustSample)[[2]])
symmetric.KL(mean1, mean2, cov1, cov2)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ multivariate }
\keyword{ distance }
