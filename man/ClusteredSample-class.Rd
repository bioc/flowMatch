\name{ClusteredSample-class}
\docType{class}
\alias{ClusteredSample-class}
\alias{ClusteredSample}
\alias{plot,flowFrame,ClusteredSample-method}
\alias{plot,ANY,ClusteredSample-method}
\alias{show,ClusteredSample-method}
\alias{summary,ClusteredSample-method}

\alias{get.size,ClusteredSample-method}
\alias{get.num.clusters}
\alias{get.num.clusters,ClusteredSample-method}
\alias{get.clusters}
\alias{get.clusters,ClusteredSample-method}
\alias{get.dimension}
\alias{get.dimension,ClusteredSample-method}
\alias{get.sample.id,ClusteredSample-method}
\alias{get.labels}
\alias{get.labels,ClusteredSample-method}



\title{ClusteredSample: A class representing a clustered FC Sample}

\description{An object of class "\code{ClusteredSample}" represents a partitioning of a sample into clusters. We model a flow cytometry sample with a mixture of cell populations where a cell population is a normally distributed cluster. An object of class "\code{ClusteredSample}" therefore stores a list of clusters and other necessary parameters.}

\section{Creating Object}{
  An object of class "\code{ClusteredSample}" can be created using the following constructor\cr \cr
  \code{ClusteredSample(labels, centers=list(), covs=list(), sample=NULL, sample.id=NA_integer_)}

\itemize{
    \item{\code{labels } {A vector of integers (from \code{1:num.clusters}) indicating the cluster to which each point is allocated. This is usually obtained from a clustering algorithm.}}
    \item{\code{centers } {A list of length \code{num.clusters} storing the centers of the clusters. The ith entry of the list \code{centers[[i]]} stores the center of the ith cluster. If not specified, the constructor estimates \code{centers} from \code{sample}.}}
    \item{\code{covs }  {A list of length \code{num.clusters} storing the covariance matrices of the clusters. The ith entry of the list \code{cov[[i]]} stores the covariance matrix of the ith cluster. If not specified, the constructor estimates \code{cov} from \code{sample}.}}
    \item \code{sample } {A matrix, data frame of observations, or object of class \code{flowFrame}.  Rows correspond to observations and columns correspond to variables. It must be passed to the constructor if either \code{centers} or \code{cov} is unspecified; then \code{centers} or  \code{cov} is estimated from \code{sample}.}    
    \item{\code{sample.id  }} {The index of the sample (relative to other samples of a cohort).}
  }
}

\section{Slots}{
  \describe{
    An object of class "\code{ClusteredSample}" contains the following slots:
    \item{\code{num.clusters }}{The number of clusters in the sample.}
    \item{\code{labels }}{A vector of integers (from range \code{1:num.clusters}) indicating the cluster to which each point is assigned to. For example, \code{labels[i]=j} means that the ith element (cell) is assigned to the jth cluster.}
    \item{\code{dimension}}{Dimensionality of the sample (number of columns in data matrix).}
    \item{\code{clusters}}{A list of length \code{num.clusters} storing the cell populations. Each cluster is stored as an object of class \code{\link{Cluster}}.} 
    \item{\code{size }}{Number of cells in the sample (summation of all cluster sizes). }
    \item{\code{sample.id}}{integer, denoting the index of the sample (relative to other samples of a cohort). Default is \code{NA_integer_}}
}
}


\section{Accessors}{

\describe{
    All the slot accessor functions take an object of class \code{ClusteredSample}. I show usage of the first accessor function. Other functions can be called similarly.
    \item{\code{get.size}:}{Returns the number of cells in the sample (summation of all cluster sizes). 
     
     \emph{Usage:} \code{   get.size(object)}
     
     here \code{object} is a \code{ClusteredSample} object.
     }
   
    \item{\code{get.num.clusters }}{Returns the number of clusters in the sample.}
    \item{\code{get.labels }}{Returns the cluster labels for each cell. For example, \code{labels[i]=j} means that the ith element (cell) is assigned to the jth cluster.}
    \item{\code{get.dimension}}{Returns the dimensionality of the sample (number of columns in data matrix).}
    \item{\code{get.clusters}}{Returns the list of clusters in this sample. Each cluster is stored as an object of class \code{\link{Cluster}}.} 
    \item{\code{get.sample.id}}{Returns the index of the sample (relative to other samples of a cohort). }

    }
}


\section{Methods}{
  \describe{
  \item{show}{

      Display details about the \code{ClusteredSample} object.
      
    }

    \item{summary}{Return descriptive summary for the \code{ClusteredSample} object.

      \emph{Usage:} \code{   summary(ClusteredSample)}
    }
    \item{plot}{We plot a sample by bivariate scatter plots where different clusters are shown in different colors. \cr 

     \emph{Usage:} 
     
     \code{   plot(sample, ClusteredSample, ...)}
      

    the arguments of the plot function are:
    \itemize{
      \item{\code{sample: } {A matrix, data.frame or an object of class  \code{flowFrame} representing an FC sample.}}
      \item{\code{ClusteredSample: } { An object of class \code{ClusteredSample} storing the clustering of the sample.}}
      \item{\code{... } {Other usual plotting related parameters.}} 
    }
  }
}
}


\author{ Ariful Azad }

\seealso{
  
  \code{\link{Cluster}}
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
## Create ClusteredSample object  (Option 1 )
## without specifying centers and covs
## we need to pass FC sample for paramter estimation
## ------------------------------------------------

clustSample = ClusteredSample(labels=cluster.labels, sample=sample)

## ------------------------------------------------
## Create ClusteredSample object  (Option 2)
## specifying centers and covs 
## no need to pass the sample
## ------------------------------------------------

centers = list()
covs = list()
num.clusters = nrow(km$centers)
for(i in 1:num.clusters)
{
  centers[[i]] = km$centers[i,]
  covs[[i]] = cov(sample[cluster.labels==i,])
}
# Now we do not need to pass sample
ClusteredSample(labels=cluster.labels, centers=centers, covs=covs)

## ------------------------------------------------
## Show summary and plot a clustered sample
## ------------------------------------------------

summary(clustSample)
plot(sample, clustSample)

}

\keyword{ cluster }
\keyword{ multivariate }


