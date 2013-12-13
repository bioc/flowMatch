\name{dist.template}
\alias{dist.template}
\title{
Dissimilarity between a pair of FC templates
}
\description{
Compute the dissimilarity between a pair of FC templates by using Mixed Edge Cover (MEC) algorithm.
}
\usage{
dist.template(template1, template2, dist.type='Mahalanobis', unmatch.penalty=999999)
}

\arguments{
  \item{template1 }{ an object of class \code{\link{Template}} containing cell populations from template 1.}
  \item{template2 }{ an object of class \code{\link{Template}} containing cell populations from template 2.}
  \item{dist.type}{character, indicating the method with which the dissimilarity between a pair of meta-clusters is computed. Supported dissimilarity measures are:  'Mahalanobis', 'KL' and 'Euclidean'.}
  \item{unmatch.penalty }{A numeric value denoting the penalty for leaving a meta-cluster unmatched. This parameter should be already known or be estimated empirically estimated from data (see the reference for a discussion). Default is set to a very high value so that no meta-cluster remains unmatched. }
}
\details{
We used a robust version of matching called Mixed Edge Cover (MEC) to match meta-clusters across a pair of templates.
MEC allows a meta-cluster to be matched with zero, one or more than one meta-clusters in a paired template. 
The cost of an MEC solution is equal to the summation of dissimilarities of the matched meta-clusters and penalty for the unmatched meta-clusters.
The MEC algorithm finds an optimal solution by minimizing the cost of MEC, which is then used as dissimilarity between a pair of templates. 
}
\value{
\code{dist.template} returns a numeric value representing dissimilarity between a pair of templates. This value is equal to the summation of dissimilarities of the matched meta-clusters and penalty for the unmatched meta-clusters.
}
\references{
Azad, Ariful and Langguth, Johannes and Fang, Youhan and Qi, Alan and Pothen, Alex (2010), Identifying rare cell populations in comparative flow cytometry; Algorithms in Bioinformatics, Springer, 162-175.
}
\author{Ariful Azad}


\seealso{
  \code{\link{Template}}, \code{\link{match.clusters}}
}
\examples{
## ------------------------------------------------
## load data and retrieve two templates
## ------------------------------------------------

library(healthyFlowData)
data(hd)

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

D = dist.template(template1, template2, dist.type='Mahalanobis', unmatch.penalty=999999)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ multivariate }
\keyword{ matching }
