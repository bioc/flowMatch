\name{flowMatch-package}
\alias{flowMatch-package}
\alias{flowMatch}
\docType{package}
\title{
Matching cell populations and building meta-clusters and templates from a collection of FC samples.
}
\description{
Matching cell populations and building meta-clusters and templates from a collection of FC samples.
}
\details{
\tabular{ll}{
Package: \tab flowMatch\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-08-01\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

}
\author{
Ariful Azad <aazad@purdue.edu>
}


\references{
  Azad, Ariful and Pyne, Saumyadipta and Pothen, Alex (2012), Matching phosphorylation response patterns of antigen-receptor-stimulated T cells via flow cytometry; BMC Bioinformatics, 13 (Suppl 2), S10.
  
  Azad, Ariful and Langguth, Johannes and Fang, Youhan and Qi, Alan and Pothen, Alex (2010), Identifying rare cell populations in comparative flow cytometry; Algorithms in Bioinformatics, Springer, 162-175.
}

\keyword{ nonparametric }
\keyword{ multivariate }
\keyword{ cluster }
\keyword{ matching }
\keyword{ meta-cluster }
\keyword{ template }


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
}
