\name{Cluster-class}
\docType{class}
\alias{Cluster}
\alias{Cluster-class}
\alias{show,Cluster-method}
\alias{summary,Cluster-method}
\alias{get.size}
\alias{get.size,Cluster-method}
\alias{get.center}
\alias{get.center,Cluster-method}
\alias{get.cov}
\alias{get.cov,Cluster-method}
\alias{get.cluster.id}
\alias{get.cluster.id,Cluster-method}
\alias{get.sample.id}
\alias{get.sample.id,Cluster-method}
\alias{sample.id<-}
\alias{sample.id<-,Cluster-method}

\title{Cluster: A class representing a cell population in FC}

\description{ An object of class "\code{Cluster}" represents a cluster or a cell population. We model a cluster with a normal distribution. An object of class "\code{Cluster}" therefore represents a cluster with a mean vector, a covariance matrix and the size of the cluster.}

\section{Creating Object}{
An object of class \code{Cluster} is usually created when constructing an object of class \code{\link{ClusteredSample}}.
Unless you know exactly what you are doing, creating an object of class "\code{Cluster}" using \code{new} or using the constructor is discouraged.

An object of class "\code{Cluster}" can be created using the following constructor\cr \cr
  \code{Cluster(size, center, cov, cluster.id = NA_integer_, sample.id=NA_integer_)} \cr \cr
The arguments of the constructor bear usual meaning as described in the value section above.
}

\section{Slots}{
  \describe{

    An object of class "\code{Cluster}" contains the following slots:
    \item{\code{size}:}{An integer denoting the number of points (cells) present in the cluster.}
    \item{\code{center}:}{A numeric vector denoting the center of the cluster.}
    \item{\code{cov}:}{A matrix denoting the covariances of the underlying normal distribution of the cluster.}
    \item{\code{cluster.id}:}{The index of the cluster (relative to other clusters in same sample). Default is \code{NA_integer_} .}
    \item{\code{sample.id}:}{The index of sample in which the cluster belongs to. Default is \code{NA_integer_} .}
  }
}


\section{Accessors}{

\describe{
    All the slot accessor functions take an object of class \code{Cluster}. I show usage of the first accessor function. Other functions can be called similarly.
    \item{\code{get.size}:}{Returns the number of cells in the cluster. 
     
     \emph{Usage:} \code{   get.size(object)}
     
     here \code{object} is a \code{Cluster} object.
     }
   
    \item{\code{get.center}:}{Returns the center of the cluster.}
    \item{\code{get.cov}:}{Returns the covariances matrix of the cluster.}
    \item{\code{get.cluster.id}:}{Returns  the index of the cluster (relative to other clusters in same sample).}
    \item{\code{get.sample.id}:}{Returns the index of sample in which the cluster belongs to.}
    \item{\code{sample.id<-}:}{Set the index of sample in which the cluster belongs to.}
    }
}

\section{Methods}{
  \describe{
  \item{show}{

      Display details about the \code{Cluster} object.
      
    }

    \item{summary}{Return descriptive summary for each \code{Cluster} object.

      \emph{Usage:} \code{   summary(Cluster)}
    }
  }
}



\author{ Ariful Azad }

\seealso{
  
  \code{\link{ClusteredSample}}
  }
  

\examples{
## An object of class "Cluster"" is usually created when constructing a "ClusteredSample". 
## Unless you know exactly what you are doing, creating an object of class "Cluster" 
## using new or using the constructor is discouraged.

## ------------------------------------------------
## load data and retrieve a sample
## ------------------------------------------------

library(healthyFlowData)
data(hd)
sample = exprs(hd.flowSet[[1]])

## ------------------------------------------------
## cluster sample using kmeans algorithm
## and retrive the parameters of  the first cluster
## ------------------------------------------------

km = kmeans(sample, centers=4, nstart=20)
center1 = km$centers[1,]
# compute the covariance matrix of the first cluster
cov1 = cov(sample[km$cluster==1,])
size1 = length(which(km$cluster==1))

## ------------------------------------------------
## Create an object of class "Cluster" 
## and show summary
## ------------------------------------------------

clust = Cluster(size=size1, center=center1, cov=cov1)
summary(clust)
}

\keyword{ cluster }
\keyword{ multivariate }

