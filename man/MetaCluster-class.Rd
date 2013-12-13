\name{MetaCluster-class}
\docType{class}
\alias{MetaCluster-class}
\alias{MetaCluster}
\alias{summary,MetaCluster-method}
\alias{show,MetaCluster-method}
\alias{plot,MetaCluster-method}
\alias{plot,MetaCluster,missing-method}

\alias{get.size,MetaCluster-method}
\alias{get.num.clusters,MetaCluster-method}
\alias{get.clusters,MetaCluster-method}
\alias{get.center,MetaCluster-method}
\alias{get.cov,MetaCluster-method}

\title{MetaCluster: An S4 class representing a meta-cluster (collection of biologically similar clusters). }

\description{An object of class "\code{MetaCluster}" represents a collection of biologically similar clusters from a set of FC samples. A meta-cluster is formed by matching clusters across samples and merging the matched clusters. An object of class "\code{ClusteredSample}" stores the estimated parameter of the whole meta-cluster as well as a list of clusters participating in the meta-cluster.}



\section{Creating Object}{
An object of class \code{MetaCluster} is usually created when constructing an object of class \code{\link{Template}}.
Unless you know exactly what you are doing, creating an object of class "\code{MetaCluster}" using \code{new} or using the constructor is discouraged.

An object of class "\code{MetaCluster}" can be created using the following constructor\cr \cr
  \code{MetaCluster(clusters)} where the argument "clusters" is  a list of object of class \code{\link{Cluster}} from which the meta-cluster is created.
}


\section{Slots}{
  \describe{
%%  ~Describe the value returned
An object of class "\code{MetaCluster}" contains the following slots:
\item{\code{num.clusters}}{The number of clusters in the meta-cluster.}
\item{\code{clusters}}{A list of length \code{num.clusters} storing the clusters (cell populations) participating in this meta-cluster. Each cluster is stored as an object of class \code{\link{Cluster}}.}
\item{\code{size }}{Number of cells in the meta-cluster (summation of all cluster sizes). }
\item{\code{center}}{A numeric vector denoting the center of the meta-cluster.}
\item{\code{cov}}{A matrix denoting the covariances of the underlying normal distribution of the meta-cluster.}
}
}


\section{Accessors}{

\describe{
    All the slot accessor functions take an object of class \code{MetaCluster}. I show usage of the first accessor function. Other functions can be called similarly.
    \item{\code{get.size}:}{The number of cells in the meta-cluster(summation of all cluster sizes). 
     
     \emph{Usage:} \code{   get.size(object)}
     
     here \code{object} is a \code{MetaCluster} object.
     }
    \item{\code{get.num.clusters}}{Returns the number of clusters in the meta-cluster.}
    \item{\code{get.clusters}}{Returns the list of clusters (cell populations) participating in this meta-cluster. Each cluster is stored as an object of class \code{\link{Cluster}}.}
    \item{\code{get.size }}{Returns the number of cells in the meta-cluster (summation of all cluster sizes). }
    \item{\code{get.center}}{Returns the center of the meta-cluster.}
    \item{\code{get.cov}}{Returns the covariances matrix of the meta-cluster.}

    }
}


\section{Methods}{
  \describe{
  \item{show}{

      Display details about the \code{Metacluster} object.
      
    }

    \item{summary}{Return descriptive summary for the \code{MetaCluster} object.

      \emph{Usage:} \code{   summary(MetaCluster)}
    }
    \item{plot}{We plot a meta-cluster as a contour plot of the distribution of the underlying clusters or the combined meta-cluster. 
We consider cells in clusters or in the meta-cluster are normally distributed and represent the distribution with ellipsoid.
The axes of an ellipsoid is estimated from the eigen values and eigen vectors of the covariance matrix ("Applied Multivariate Statistical Analysis" by R. Johnson and D. Wichern, 5th edition, Prentice hall).
We then plot the bi-variate projection of the ellipsoid as 2-D ellipses. \cr 

     \emph{Usage:} 
     
     \code{plot(mc, alpha=.05, plot.mc=FALSE, ...)}

    the arguments of the plot function are:
    \itemize{
      \item{\code{mc } {An object of class  \code{MetaCluster} for which the plot function is invoked.}}
      \item{\code{alpha } { (1-alpha)*100\% quantile of the distribution of the clusters or meta-cluster is plotted.}}
      \item{\code{plot.mc }  { TRUE/FALSE, when TRUE the functions draws contour of the combined meta-cluster and when FALSE the function draws the contours of the individual clusters.}}
      \item{\code{... } {Other usual plotting related parameters.}} 
    }
  }
}
}

\author{ Ariful Azad }

\references{
  Azad, Ariful and Pyne, Saumyadipta and Pothen, Alex (2012), Matching phosphorylation response patterns of antigen-receptor-stimulated T cells via flow cytometry; BMC Bioinformatics, 13 (Suppl 2), S10.
}

\seealso{
  
  \code{\link{Cluster}}, \code{\link{Template}}
  }
  

\examples{
## ------------------------------------------------
## load data
## ------------------------------------------------

library(healthyFlowData)
data(hd)

## ------------------------------------------------
## Retrieve each sample, cluster it and store the
## clustered samples in a list
## ------------------------------------------------

cat('Clustering samples: ')
clustSamples = list()
for(i in 1:length(hd.flowSet))
{
  cat(i, ' ')
  sample1 = exprs(hd.flowSet[[i]])
  clust1 = kmeans(sample1, centers=4, nstart=20)
  cluster.labels1 = clust1$cluster
  clustSample1 = ClusteredSample(labels=cluster.labels1, sample=sample1)
  clustSamples = c(clustSamples, clustSample1)
}

## ------------------------------------------------
## Create a template from the list of clustered samples and retrieve the meta-clusters
## ------------------------------------------------

template = create.template(clustSamples)
#retrieve meta-clusters from template
mc = get.metaClusters(template)[[1]]
summary(mc)
# plot all participating cluster in this meta-cluster
plot(mc) 
# plot the outline of the combined  meta-cluster
plot(mc, plot.mc=TRUE)
}

\keyword{ cluster }
\keyword{ metacluster }
\keyword{ multivariate }
