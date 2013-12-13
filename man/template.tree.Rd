\name{template.tree}
\alias{template.tree}
\alias{template.tree,Template-method}
\title{
  Plot the hierarchy of samples established while creating a template-tree
}
\description{
  All samples within a template are organized as binary tree.  This function plots the hierarchy of samples established while creating a template-tree.
}
\usage{
  template.tree(object, ...)
}

\arguments{
    \item{object }{An object of class \code{\link{Template}}. The working examples describe how a template is created from a collection of FC samples.}
    \item{... }{Other usual plotting related parameters.} 
   }

\value{
  Returns a tree object of class \code{hclust} storing the hierarchy of the samples in the template.
}
\references{
  Azad, Ariful and Pyne, Saumyadipta and Pothen, Alex (2012), Matching phosphorylation response patterns of antigen-receptor-stimulated T cells via flow cytometry; BMC Bioinformatics, 13 (Suppl 2), S10.
}
\author{Ariful Azad}


\seealso{
  \code{\link{Template}}, \code{\link{create.template}}
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

  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ multivariate }
\keyword{ matching }
