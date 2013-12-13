\name{mahalanobis.dist}
\alias{mahalanobis.dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mahanalobis Distance
}
\description{
Compute the Mahalanobis distance between a pair of normally distributed clusters.
}
\usage{
mahalanobis.dist(mean1, mean2, cov1, cov2, n1, n2)
}

\arguments{
  \item{mean1}{mean vector of length \code{p} for cluster 1, where \code{p} is the dimension of the clusters.}
  \item{mean2}{mean vector of length \code{p} for cluster 2.}
  \item{cov1}{\code{pxp} covariance matrix for cluster 1.}
  \item{cov2}{\code{pxp} covariance matrix for cluster 2.}
  \item{n1}{number of cells (points) in cluster 1.}
  \item{n2}{number of cells (points) in cluster 2.}
}
\details{
Consider two \code{p}-dimensional, normally distributed clusters with centers \eqn{\mu1, \mu2} and covariance matrices \eqn{\Sigma1, \Sigma2}. Assume the size of the clusters are \code{n1} and \code{n2} respectively.
We compute the Mahalanobis distance \code{d12} between the clusters as follows:
\deqn{\Sigma = ( (n1-1) * \Sigma1 + (n2-1) * \Sigma2) / (n1+n2-2) }
  \deqn{  d12 = sqrt( t(\mu1-\mu2) * \Sigma^(-1) * (\mu1-\mu2)) }
The dimension of the clusters must be same.
}
\value{
\code{mahalanobis.dist} returns a numeric value measuring the Mahalanobis distance between a pair of normally distributed clusters. 
}
\references{
McLachlan, GJ (1999) Mahalanobis distance; Journal of Resonance 4(6), 20--26.
}
\author{Ariful Azad}


\seealso{
 \code{\link{symmetric.KL}}, \code{\link{dist.cluster}}
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
mahalanobis.dist(mean1, mean2, cov1, cov2, n1, n2)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ multivariate }
\keyword{ distance }

