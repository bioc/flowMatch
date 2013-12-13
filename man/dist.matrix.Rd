\name{dist.matrix}
\alias{dist.matrix}
\title{
Dissimilarity matrix between each pair of clusters/meta-clusters across a pair of samples/templates
}
\description{
Calculate a matrix storing the dissimilarities between each pair of clusters (meta-clusters) across a pair of samples (templates) \code{S1} and \code{S2}. \code{(i,j)th} entry of the matrix stores dissimilarity between \code{i-th} cluster (meta-cluster) from \code{S1} and the \code{j-th} cluster (meta-cluster) from \code{S2}.
}
\usage{
dist.matrix(object1,object2, dist.type = 'Mahalanobis')
}

\arguments{
  \item{object1 }{ an object of class \code{\link{ClusteredSample}} or \code{\link{Template}}. }
  \item{object2 }{ an object of class \code{\link{ClusteredSample}} or \code{\link{Template}}. }
  \item{dist.type}{character, indicating the method with which the dissimilarity between a pair of clusters (meta-clusters) is computed. Supported dissimilarity measures are:  'Mahalanobis', 'KL' and 'Euclidean', with the default is set to 'Mahalanobis' distance. }
}

\details{
Consider two FC samples/templates \code{S1} and \code{S2} with \code{k1} and \code{k2} clusters/meta-clusters. The dissimilarity between each pair of cluster (meta-clusters) across \code{S1} and \code{S2} is computed and stored in a (\code{k1 x k2}) matrix. The dissimilarity between \code{i-th} cluster (meta-cluster) from \code{S1} and \code{j-th} cluster (meta-cluster) from \code{S2} is computed using function \code{\link{dist.cluster}}. 
}

\value{
\code{dist.matrix} function returns a (\code{k1 x k2}) matrix where \code{k1} and \code{k2} are the number of clusters (meta-clusters) in the first and the second samples (templates) respectively. \code{(i,j)}th entry of the matrix contains the dissimilarity between the \code{i-th} cluster (meta-cluster) from sample1 (template1) and the \code{j-th} cluster (meta-cluster) from sample2 (template2).
}

\author{Ariful Azad}

\seealso{
\code{\link{dist.cluster}}
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
## and compute the Mahalanobis distance between 
## each pair of clsuters and save it in a matrix
## ------------------------------------------------

clustSample1 = ClusteredSample(labels=cluster.labels1, sample=sample1)
clustSample2 = ClusteredSample(labels=cluster.labels2, sample=sample2)
## compute the dissimilarity matrix
DM = dist.matrix(clustSample1, clustSample2, dist.type='Mahalanobis')
print(DM)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ multivariate }
\keyword{ distance }
