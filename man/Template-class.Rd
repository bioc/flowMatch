\name{Template-class}
\docType{class}
\alias{Template-class}
\alias{Template}
\alias{plot,Template-method}
\alias{plot,Template,ANY-method}
\alias{plot,Template,missing-method}
\alias{show,Template-method}
\alias{summary,Template-method}

\alias{get.size,Template-method}
\alias{get.num.metaclusters}
\alias{get.num.metaclusters,Template-method}
\alias{get.dimension,Template-method}
\alias{get.template.id}
\alias{get.template.id,Template-method}
\alias{get.metaClusters}
\alias{get.metaClusters,Template-method}
\alias{get.tree}
\alias{get.tree,Template-method}

%\alias{template.tree,Template-method}

\title{Template: An S4 class representing a template of a group of FC Samples.}

\description{An object of class "\code{Template}" summarizes a group of samples belonging to same biological-class with a class-template. A template is represented by a collection of meta-clusters (\code{\link{MetaCluster}}) created from samples of same class. An object of class "\code{Template}" therefore stores a list of \code{\link{MetaCluster}} objects and other necessary parameters.}

\section{Creating Object}{
  An object of class "\code{Template}" can be created using the function \code{\link{create.template}} : \cr \cr
  \code{create.template(clustSamples, dist.type = "Mahalanobis", unmatch.penalty=999999, template.id = NA_integer_)}.
  
  The arguments to the \code{create.template} function is described below:

\itemize{
    \item{\code{clustSamples: }} {A list of \code{\link{ClusteredSample}} objects from which the template is created. The working examples describe how this objects are created by clustering FC samples.}
    \item{\code{dist.type: }}{character, indicating the method with which the dissimilarity between a pair of clusters is computed. Supported dissimilarity measures are:  'Mahalanobis', 'KL' and 'Euclidean'. If this argument is not passed then 'Mahalanobis' distance is used by default.}
    \item{\code{unmatch.penalty: }}{A numeric value denoting the penalty for leaving a cluster unmatched. This parameter should be already known or be estimated empirically estimated from data (see the reference for a discussion). Default is set to a very high value so that no cluster remains unmatched. }
    \item{\code{template.id: }}{integer, denoting the index of the template (relative to other template). Default is \code{NA_integer_}}
  }
}

\section{Slots}{
  \describe{
    \item{\code{num.metaclusters}:}{The number of meta-clusters in the template.}
    \item{\code{metaClusters}:}{A list of length \code{num.metaclusters} storing the meta-clusters. Each meta-cluster is stored as an object of class \code{\link{MetaCluster}}.}
    \item{\code{dimension}:}{Dimensionality of the samples from which the template is created.}
    \item{\code{size}:}{Number of cells in the template (summation of all meta-cluster sizes). }
    \item{\code{tree}:}{ A list (similar to an \code{hclust} object) storing the hierarchy of the samples in a template.}
    \item{\code{template.id}:}{integer, denoting the index of the template (relative to other templates). Default is \code{NA_integer_}}
}
}


\section{Accessors}{

\describe{
    All the slot accessor functions take an object of class \code{Template}. I show usage of the first accessor function. Other functions can be called similarly.
    \item{\code{get.size}:}{Number of cells in the template (summation of all meta-cluster sizes). 
     
     \emph{Usage:} \code{   get.size(object)}
     
     here \code{object} is a \code{Template} object.
     }
    \item{\code{get.num.metaclusters}:}{Returns the number of meta-clusters in the template.}
    \item{\code{get.metaClusters}:}{Returns a list of length \code{num.metaclusters} storing the meta-clusters. Each meta-cluster is stored as an object of class \code{\link{MetaCluster}}.}
    \item{\code{get.dimension}:}{Returns the  dimensionality of the samples from which the template is created.}
    \item{\code{get.tree}:}{Returns a \code{hclust} object storing the hierarchy of the samples in a template.}
    \item{\code{get.template.id}:}{Returns the index of the template (relative to other templates). }
    }
}

\section{Methods}{
  \describe{
  \item{show}{

      Display details about the \code{Template} object.
      
    }

    \item{summary}{Return descriptive summary for each \code{MetaCluster} of a Template.

      \emph{Usage:} \code{   summary(Template)}
    }
    
    \item{plot}{We plot a template as a collection of bivariate contour plots of its meta-clusters. 
To plot each meta-cluster we consider the clusters within the meta-cluster normally distributed and represent each cluster with an ellipsoid. The axes of an ellipsoid is estimated from the eigen values and eigen vectors of the covariance matrix of a cluster ("Applied Multivariate Statistical Analysis" by R. Johnson and D. Wichern, 5th edition, Prentice hall). We then plot the bivariate projection of the ellipsoid as 2-D ellipses. \cr \cr


      \emph{Usage:} 
      
      \code{plot(template, alpha=.05, plot.mc=FALSE, color.mc=NULL, colorbysample=FALSE, ...)}
      
      the arguments of the plot function are:
      \itemize{
        \item{\code{template: }}{An object of class  \code{Template} for which the plot function is invoked.}
        \item{\code{alpha: }}{ (1-alpha)*100\% quantile of the distribution of the clusters or meta-cluster is plotted.}
        \item{\code{plot.mc: }}{ TRUE/FALSE, when TRUE the functions draws contour of the combined meta-cluster and when FALSE the function draws the contours of the individual clusters.}
        \item{\code{color.mc: }}{ A character vector of length \code{num.metaclusters} denoting the colors to be used to draw the contours. The ith color of this vector is used to draw the ellipses denoting clusters in the ith meta-cluster or the combined ith meta-cluster (depending on the argument \code{plot.mc}). By default an empty vector is passed and then an arbitrary color is used to draw each meta-cluster.}
        \item{\code{colorbysample: }}{ TRUE/FALSE, when TRUE the functions draws clusters from same samples in a single color and when FALSE the function draws meta-clusters in a single color.}
    
        \item{\code{... : }}{Other usual plotting related parameters.}
      }
    }
    
    \item{template.tree}{Plot the hierarchy of samples established while creating the template-tree. See \code{\link{template.tree}}
    }
  }
}

\author{ Ariful Azad }
\references{
  Azad, Ariful and Pyne, Saumyadipta and Pothen, Alex (2012), Matching phosphorylation response patterns of antigen-receptor-stimulated T cells via flow cytometry; BMC Bioinformatics, 13 (Suppl 2), S10.
}

\seealso{
  
  \code{\link{MetaCluster}}, \code{\link{ClusteredSample}}, \code{\link{create.template}}, \code{\link{template.tree}}
  }
  

\examples{
## ------------------------------------------------
## load data
## ------------------------------------------------

library(healthyFlowData)
data(hd)

## ------------------------------------------------
## Retrieve each sample, clsuter it and store the
## clustered samples in a list
## ------------------------------------------------
set.seed(1234) # for reproducable clustering 
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
## Create a template from the list of clustered samples and plot functions
## ------------------------------------------------

template = create.template(clustSamples)
summary(template)

## plot the tree denoting the hierarchy of the samples in a template
tree = template.tree(template)

## plot the template in terms of the meta-clusters
## option-1 (default): plot contours of each cluster of the meta-clusters
plot(template)

## option-2: plot contours of each cluster of the meta-clusters with defined color
plot(template, color.mc=c('blue','black','green3','red'))

## option-3: plot contours of the meta-clusters with defined color
plot(template, plot.mc=TRUE, color.mc=c('blue','black','green3','red'))


## option-4: plot contours of each cluster of the meta-clusters with different colors for different samples
plot(template, colorbysample=TRUE)
}

\keyword{ cluster }
\keyword{ metacluster }
\keyword{ template }
\keyword{ multivariate }

