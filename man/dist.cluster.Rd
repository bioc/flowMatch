\name{dist.cluster}
\alias{dist.cluster}
\title{
Dissimilarity between a pair of clusters
}
\description{
Calculate the dissimilarity between a pair of cell populations (clusters) from the distributions of the clusters.
}
\usage{
dist.cluster(cluster1,cluster2, dist.type = 'Mahalanobis')
}

\arguments{
  \item{cluster1 }{ an object of class \code{\link{Cluster}} representing the distribution parameters of the first cluster.}
  \item{cluster2 }{ an object of class \code{\link{Cluster}} representing the distribution parameters of the second cluster.}
  \item{dist.type}{ character, indicating the method with which the dissimilarity between a pair of clusters is computed. Supported dissimilarity measures are:  'Mahalanobis', 'KL' and 'Euclidean'.}
}

\value{
\code{dist.cluster} returns a numeric value denoting the dissimilarities between a pair of cell populations (clusters).
}

\details{
  Consider two \code{p}-dimensional, normally distributed clusters with centers \eqn{\mu1, \mu2} and covariance matrices \eqn{\Sigma1, \Sigma2}. Assume the size of the clusters are \code{n1} and \code{n2} respectively. We compute the dissimilarity \code{d12} between the clusters as follows:
  \enumerate{

  \item{If dist.type='Mahalanobis':  we compute the dissimilarity \code{d12} with the Mahalanobis distance between the distributions of the clusters. 
  \deqn{\Sigma = ( (n1-1) * \Sigma1 + (n2-1) * \Sigma2) / (n1+n2-2) }
  \deqn{  d12 = sqrt( t(\mu1-\mu2) * \Sigma^(-1) * (\mu1-\mu2)) }

    }
    \item{If dist.type='KL':  we compute the dissimilarity \code{d12} with the Symmetrized Kullback-Leibler divergence between the distributions of the clusters. Note that KL-divergence is not symmetric in its original form. We converted it symmetric by averaging both way KL divergence. The symmetrized KL-divergence is not a metric because it does not satisfy triangle inequality.
  \deqn{d12 = 1/4 * ( t(\mu2 - \mu1) * ( \Sigma1^(-1) + \Sigma2^(-1) ) * (\mu2 - \mu1) +
    trace(\Sigma1/\Sigma2 + \Sigma2/\Sigma1) + 2p ) }
    }
  \item{If dist.type='Euclidean':  we compute the dissimilarity \code{d12} with the Euclidean distance between the centers of the clusters. 
  \deqn{d12 =sqrt(\sum(\mu1-\mu2)^2 ) }
    }
  }
  
The dimension of the clusters must be same. 

}
\references{
McLachlan, GJ (1999) Mahalanobis distance; Journal of Resonance 4(6), 20--26.

Abou--Moustafa, Karim T and De La Torre, Fernando and Ferrie, Frank P (2010) Designing a Metric for the Difference between Gaussian Densities; Brain, Body and Machine, 57--70.
}
\author{Ariful Azad}


\seealso{
\code{\link{mahalanobis.dist}, \link{symmetric.KL}, \link{dist.matrix}}
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
clust1 = get.clusters(clustSample)[[1]]
clust2 = get.clusters(clustSample)[[2]]
dist.cluster(clust1, clust2, dist.type='Mahalanobis')
dist.cluster(clust1, clust2, dist.type='KL')
dist.cluster(clust1, clust2, dist.type='Euclidean')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ multivariate }
\keyword{ distance }
