\name{create.template}
\alias{create.template}
\title{
  Creating a template of a collection of FC samples
}
\description{
  Create an object of class \code{\link{Template}} summarizes a group of samples belonging to same biological-class with a class-template. A template is represented by a collection of meta-clusters (\code{\link{MetaCluster}}) created from samples of same class. An object of class \code{\link{Template}} therefore stores a list of \code{\link{MetaCluster}} objects and other necessary parameters.
}
\usage{
  create.template(clustSamples, dist.type = "Mahalanobis", unmatch.penalty=999999, template.id = NA_integer_)
}

\arguments{
    \item{clustSamples }{A list of \code{\link{ClusteredSample}} objects from which the template is created. The working examples describe how this objects are created by clustering FC samples.}
    \item{dist.type }{character, indicating the method with which the dissimilarity between a pair of clusters is computed. Supported dissimilarity measures are:  'Mahalanobis', 'KL' and 'Euclidean'. If this argument is not passed then 'Mahalanobis' distance is used by default.}
    \item{unmatch.penalty }{A numeric value denoting the penalty for leaving a cluster unmatched. This parameter should be already known or be estimated empirically estimated from data (see the reference for a discussion). Default is set to a very high value so that no cluster remains unmatched. }
    \item{template.id }{integer, denoting the index of the template (relative to other template). Default is NA_integer_}
}
\details{
  An object of class \code{\link{Template}} summarizes a group of samples belonging to same biological-class with a class-specific template. A template is represented by a collection of meta-clusters (\code{\link{MetaCluster}}) created from samples of same class. An object of class \code{\link{Template}} therefore stores a list of \code{\link{MetaCluster}} objects and other necessary parameters.
}
\value{
  \code{dist.sample} returns a numeric value representing dissimilarity between a pair of samples. This value is equal to the summation of dissimilarities of the matched clusters and penalty for the unmatched clusters.
}
\references{
  Azad, Ariful and Pyne, Saumyadipta and Pothen, Alex (2012), Matching phosphorylation response patterns of antigen-receptor-stimulated T cells via flow cytometry; BMC Bioinformatics, 13 (Suppl 2), S10.
}
\author{Ariful Azad}


\seealso{
  \code{\link{Template}}, \code{\link{MetaCluster}}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ multivariate }
\keyword{ matching }
