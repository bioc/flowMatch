\name{dist.sample}
\alias{dist.sample}
\title{
Dissimilarity between a pair of clustered FC samples
}
\description{
Compute the dissimilarity between a pair of clustered FC samples by using Mixed Edge Cover (MEC) algorithm.
}
\usage{
dist.sample(clustSample1, clustSample2, dist.type='Mahalanobis', unmatch.penalty=999999)
}

\arguments{
  \item{clustSample1 }{ an object of class \code{\link{ClusteredSample}} containing cell populations from sample 1.}
  \item{clustSample2 }{ an object of class \code{\link{ClusteredSample}} containing cell populations from sample 2.}
  \item{dist.type}{character, indicating the method with which the dissimilarity between a pair of clusters is computed. Supported dissimilarity measures are:  'Mahalanobis', 'KL' and 'Euclidean'.}
  \item{unmatch.penalty }{A numeric value denoting the penalty for leaving a cluster unmatched. This parameter should be already known or be estimated empirically estimated from data (see the reference for a discussion). Default is set to a very high value so that no cluster remains unmatched. }
}
\details{
We used a robust version of matching called Mixed Edge Cover (MEC) to match clusters across a pair of samples.
MEC allows a cluster to be matched with zero, one or more than one clusters in a paired sample. 
The cost of an MEC solution is equal to the summation of dissimilarities of the matched clusters and penalty for the unmatched clusters.
The MEC algorithm finds an optimal solution by minimizing the cost of MEC, which is then used as dissimilarity between a pair of samples. 
}
\value{
\code{dist.sample} returns a numeric value representing dissimilarity between a pair of samples. This value is equal to the summation of dissimilarities of the matched clusters and penalty for the unmatched clusters.
}
\references{
Azad, Ariful and Langguth, Johannes and Fang, Youhan and Qi, Alan and Pothen, Alex (2010), Identifying rare cell populations in comparative flow cytometry; Algorithms in Bioinformatics, Springer, 162-175.
}
\author{Ariful Azad}


\seealso{
  \code{\link{ClusteredSample}}, \code{\link{match.clusters}}
}
\examples{
## ------------------------------------------------
## load data and retrieve two samples
## ------------------------------------------------

library(healthyFlowData)
data(hd)
sample1 = exprs(hd.flowSet[[1]])
sample2 = exprs(hd.flowSet[[2]])

## ------------------------------------------------
## cluster sample using kmeans algorithm
## ------------------------------------------------

clust1 = kmeans(sample1, centers=4, nstart=20)
clust2 = kmeans(sample2, centers=4, nstart=20)
cluster.labels1 = clust1$cluster
cluster.labels2 = clust2$cluster

## ------------------------------------------------
## Create ClusteredSample object  
## and compute dissimilarity between two clustered samples
## using the mixed edge cover algorithm
## ------------------------------------------------

clustSample1 = ClusteredSample(labels=cluster.labels1, sample=sample1)
clustSample2 = ClusteredSample(labels=cluster.labels2, sample=sample2)
D = dist.sample(clustSample1, clustSample2, dist.type='Mahalanobis', unmatch.penalty=999999)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ multivariate }
\keyword{ matching }
