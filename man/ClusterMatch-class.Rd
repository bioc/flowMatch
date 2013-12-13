\name{ClusterMatch-class}
\docType{class}
\alias{ClusterMatch-class}
\alias{ClusterMatch}
\alias{show,ClusterMatch-method}
\alias{summary,ClusterMatch-method}
\alias{get.match12}
\alias{get.match12,ClusterMatch-method}
\alias{get.match21}
\alias{get.match21,ClusterMatch-method}
\alias{get.matching.cost}
\alias{get.matching.cost,ClusterMatch-method}
\alias{get.unmatch.penalty}
\alias{get.unmatch.penalty,ClusterMatch-method}


\title{ClusterMatch: A class representing matching of cluster/meta-clusters across a pair of FC samples/templates}

\description{An object of class "\code{ClusterMatch}" represents matching of cluster/meta-clusters across a pair of FC samples/templates. A cluster (meta-cluster) from a sample (template) can match to zero, one or more than one cluster (meta-clusters) in another sample (template).}

\section{Creating Object}{
  An object of class "\code{ClusterMatch}" is usually created by calling the function \code{\link{match.clusters}}: 
  
  \code{match.clusters(object1, object2, dist.type='Mahalanobis', unmatch.penalty=999999)}.
  
  Here, \code{object1} and \code{object2} are two objects of class \code{\link{ClusteredSample}} or \code{\link{Template}} between which the clusters or meta-clusters are matched. See the example section and also the \code{\link{match.clusters}} function for more details.
  
  Unless you know exactly what you are doing, creating an object of class "\code{ClusterMtach}" using \code{new} or using the constructor is discouraged.
}

\section{Slots}{
  \describe{
      Let \code{S1} and \code{S2} be two FC samples or templates with \code{k1} and \code{k2} clusters or meta-clusters respectively. Then the matching of clusters (meta-clusters) across \code{S1} and \code{S2} is represented by an object of class "\code{ClusterMatch}" that contains the following slots:

      \item{\code{match12}:}{A list of length \code{k1} where \code{match12[[i]]} stores the indices of clusters (meta-clusters)  from \code{S2} matched to the \code{i-th} clustrer (meta-cluster) of \code{S1}. \code{match12[[i]]=NULL} if the \code{i-th} cluster (meta-cluster) of \code{S1} remains unmatched, otherwise, it stores a vector of matched clusters (meta-clusters) from \code{S2}.}

      \item{\code{match21}:}{A list of length \code{k2} where \code{match21[[i]]} stores the indices of clusters (meta-clusters) from \code{S1} matched to the \code{i-th} clustrer (meta-cluster) of \code{S2}. \code{match21[[i]]=NULL} if the \code{i-th} cluster (meta-cluster) of \code{S2} remains unmatched, otherwise, it stores a vector of matched clusters (meta-clusters) from \code{S1}.}

      \item{\code{matching.cost}:}{The cost of matching clusters (meta-clusters) across the samples. It is equal to the summation of dissimilarities of the matched clusters (meta-clusters) and penalty for the unmatched clusters (meta-clusters).}

      \item{\code{unmatch.penalty}:}{A numeric value denoting the penalty for leaving a cluster (meta-cluster) unmatched. If we set it to a a very large value then no cluster (meta-cluster) remains unmatched giving an edge cover solution.}

  }
}



\section{Accessors}{

\describe{
    All the slot accessor functions take an object of class \code{ClusterMatch}. I show usage of the first accessor function. Other functions can be called similarly.
   
      \item{\code{get.match12}:}{Returns the matching from cluster in sample 1 to clusters in sample 2. See the slot description for details.
      \emph{Usage:} \code{   get.match12(object)}
     
      here \code{object} is a \code{ClusterMatch} object.
     }

      \item{\code{get.match21}:}{Returns the matching from cluster in sample 2 to clusters in sample 1. See the slot description for details.}

      \item{\code{get.matching.cost}:}{Returns the total cost of matching clusters (meta-clusters) across the pair samples/templates.}

      \item{\code{get.unmatch.penalty}:}{Returns the penalty for leaving a cluster (meta-cluster) unmatched.}

    }
}


\section{Methods}{
  \describe{
  \item{show}{

      Display details about the \code{ClusterMatch} object.
      
    }

    \item{summary}{Return descriptive summary of the matching of clusters (meta-clusters) across a pair of samples (templates). Shows both list and matrix format.

      \emph{Usage:} \code{   summary(ClusterMatch)}
    }
  }
}
\author{ Ariful Azad }

\seealso{
  
  \code{\link{match.clusters}}, \code{\link{ClusteredSample}}, \code{\link{Template}}
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

mec = match.clusters(clustSample1, clustSample2, dist.type="Mahalanobis", unmatch.penalty=99999)
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

\keyword{ cluster }
\keyword{ multivariate }
\keyword{ matching }

