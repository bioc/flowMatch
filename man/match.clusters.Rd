\name{match.clusters}
\alias{match.clusters}
\alias{match.clusters.dist}
\title{
Matching of clusters/meta-clusters across FC samples/templates
}
\description{
This function computes a matching of cluster/meta-clusters across a pair of FC samples/templates.
A cluster (meta-cluster) from a sample (template) can match to zero, one or more than one clusters (meta-clusters) in another sample (template).
}
\usage{
match.clusters(object1, object2, dist.type='Mahalanobis', unmatch.penalty=999999)

match.clusters.dist(d.matrix,unmatch.penalty=999999)
}

\arguments{
  \item{object1 }{ an object of class \code{\link{ClusteredSample}} or \code{\link{Template}}.}
  \item{object2 }{ an object of class \code{\link{ClusteredSample}} or \code{\link{Template}}.}
  \item{dist.type}{character, indicating the method with which the dissimilarity between a pair of clusters (meta-clusters) is computed. Supported dissimilarity measures are:  'Mahalanobis', 'KL' and 'Euclidean', with the default is set to 'Mahalanobis' distance.}
  \item{d.matrix }{ a matrix used only in the second defination (\code{match.clusters.dist}) of the function. \code{d.matrix} stores the dissimilarities between each pair of clusters (meta-clusters) across a pair of samples (templates) \code{S1} and \code{S2}. \code{(i,j)} entry of the matrix stores dissimilarity between \code{i-th} cluster (meta-cluster) from \code{S1} and the \code{j-th} cluster (meta-cluster) from \code{S2}. \code{d.matrix} can be computed using function \code{\link{dist.matrix}}}
 \item{unmatch.penalty }{  numeric value denoting the penalty for leaving a cluster (meta-cluster) unmatched. This parameter should be already known or be estimated empirically estimated from data (see the reference for a discussion). Default is set to a very high value so that no cluster (meta-cluster) remains unmatched. }
}
\details{
We used a robust version of matching called Mixed Edge Cover (MEC) to match clusters  (meta-clusters) across a pair of samples (templates).
MEC allows a cluster (meta-cluster) to be matched with zero, one or more than one clusters (meta-clusters) across a pair of samples (template). 
The cost of an MEC solution is equal to the summation of dissimilarities of the matched clusters (meta-clusters) and penalty for the unmatched clusters (meta-clusters).
The MEC algorithm finds an optimal solution by minimizing the cost of MEC. 
}
\value{
\code{match.clusters} returns an object of class \code{\link{ClusterMatch}} representing matching of clusters (meta-clusters) across a pair of FC samples (templates). A cluster (meta-cluster) from a sample (template) can match to zero, one or more than one cluster (meta-clusters) in another sample (template).
}
\references{
Azad, Ariful and Langguth, Johannes and Fang, Youhan and Qi, Alan and Pothen, Alex (2010), Identifying rare cell populations in comparative flow cytometry; Algorithms in Bioinformatics, Springer, 162-175.
}
\author{Ariful Azad}


\seealso{
  \code{\link{dist.matrix}}, \code{\link{ClusteredSample}}, \code{\link{Template}}
}
\examples{
## ------------------------------------------------
## load data and retrieve two samples
## ------------------------------------------------

library(healthyFlowData)
data(hd)

## **********************************************************************
## ************** first matching clusters across samples ****************
## **********************************************************************

## ------------------------------------------------
## retrieve and cluster two samples using kmeans algorithm
## ------------------------------------------------
sample1 = exprs(hd.flowSet[[1]])
sample2 = exprs(hd.flowSet[[2]])

clust1 = kmeans(sample1, centers=4, nstart=20)
clust2 = kmeans(sample2, centers=4, nstart=20)
cluster.labels1 = clust1$cluster
cluster.labels2 = clust2$cluster

## ------------------------------------------------
## Create ClusteredSample object  
## and compute mahalanobis distance between two clsuters
## ------------------------------------------------

clustSample1 = ClusteredSample(labels=cluster.labels1, sample=sample1)
clustSample2 = ClusteredSample(labels=cluster.labels2, sample=sample2)
## compute the dissimilarity matrix
DM = dist.matrix(clustSample1, clustSample2, dist.type='Mahalanobis')

## ------------------------------------------------
## Computing matching of clusteres  
## An object of class "ClusterMatch" is returned 
## ------------------------------------------------

## directly from the ClusteredSample objects: approach 1
mec = match.clusters(clustSample1, clustSample2, dist.type="Mahalanobis", unmatch.penalty=99999)
## from the dissimilarity matrix: approach 2
mec = match.clusters.dist(DM, unmatch.penalty=99999)
## show the matching 
summary(mec)


## **********************************************************************
## ************** Now matching meta-clusters across templates ***********
## **********************************************************************

## ------------------------------------------------
## Retrieve each sample, clsuter it and store the
## clustered samples in a list
## ------------------------------------------------

cat('Clustering samples: ')
clustSamples = list()
for(i in 1:10) # read 10 samples and cluster them
{
  cat(i, ' ')
  sample1 = exprs(hd.flowSet[[i]])
  clust1 = kmeans(sample1, centers=4, nstart=20)
  cluster.labels1 = clust1$cluster
  clustSample1 = ClusteredSample(labels=cluster.labels1, sample=sample1)
  clustSamples = c(clustSamples, clustSample1)
}

## ------------------------------------------------
## Create two templates each from five samples
## ------------------------------------------------

template1 = create.template(clustSamples[1:5])
template2 = create.template(clustSamples[6:10])

## ------------------------------------------------
## Match meta-clusters across templates
## ------------------------------------------------

mec = match.clusters(template1, template2, dist.type="Mahalanobis", unmatch.penalty=99999)
summary(mec)

## ------------------------------------------------
## Another example of matching meta-clusters & clusters
## across a template and a sample
## ------------------------------------------------

mec = match.clusters(template1, clustSample1, dist.type="Mahalanobis", unmatch.penalty=99999)
summary(mec)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ multivariate }
\keyword{ matching }
