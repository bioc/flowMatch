## ========================================================================================
## calculate the Mahalanobis distance between two Gaussian distributions (clusters)
## mean1 and cov1 are the parameters for cluster1 and
## mean2 and cov2 are for cluster2
## n1 and n2 are the number of cells (points) cluster1 and cluster2 respectively.
## ========================================================================================

mahalanobis.dist = function(mean1, mean2, cov1, cov2, n1, n2)
{
  ## ref (mclachlan1999mahalanobis, Springer 1999)
  if(length(mean1) != length(mean2) )
    stop("funtion mahalanobis.dist(): dimensions of the means do not match\n")
  if( nrow(cov1) != nrow(cov2) || ncol(cov1) != ncol(cov2) )
    stop("funtion mahalanobis.dist(): dimensions of the covariances do not match\n")
  sigmaCombined = ( cov1*(n1-1) + cov2*(n2-1) ) / (n1+n2-2);
  mdist = sqrt(t(mean1-mean2) %*%  solve(sigmaCombined) %*% (mean1-mean2));
  mdist = as.numeric (mdist)
  if(mdist<0)
    mdist = 0
  return (mdist);
}

## ========================================================================================
# calculate the symmetrized Kullback-Leibler divergence between two Gaussian distributions 
# mean1 and cov1 are the parameters for cluster1 and
# mean2 and cov2 are for cluster2
# note that the dimension of the two Gaussian must be same.
## ========================================================================================

symmetric.KL = function(mean1, mean2, cov1, cov2)
{
  
  if(length(mean1) != length(mean2) )
    stop("funtion KL.symmetric(): dimensions of the means do not match\n")
  if( nrow(cov1) != nrow(cov2) || ncol(cov1) != ncol(cov2) )
    stop("funtion KL.symmetric(): dimensions of the covariances do not match\n")
  # definition 1
  # the log(det) part cancell out after summation 
  #KL1 = .5* ( log(det(sigma2) / det (sigma1)) + sum(diag( solve(sigma2) %*% sigma1)) + t(mu2 - mu1) %*% solve(sigma2) %*% (mu2 - mu1) - length(mu1) );
  #KL2 = .5* ( log(det(sigma1) / det (sigma2)) + sum(diag( solve(sigma1) %*% sigma2)) + t(mu1 - mu2) %*% solve(sigma1) %*% (mu1 - mu2) - length(mu2) );
  # KL = (KL1+KL2) /2  
  #definition 2, direct calculation
  KL = .25 * ( t(mean2 - mean1) %*% (solve(cov1) + solve(cov2)) %*% (mean2 - mean1) +
                 sum(diag( ( solve(cov2) %*% cov1) + (solve(cov1) %*% cov2) )) - 2*length(mean1) )
  KL = as.numeric(KL)   
  if(KL<0)
    KL=0
  return (KL)
}

## ========================================================================================
## Compute dissimilarity between a pair of clusters or metaclusters
## Arguments:
##  cluster1, cluster2: objects of class "Cluster" or "MetaClusters"
##    dist.type: The method with which the dissimilarity is computed between clusters
##                Mahalanobis, KL and Euclidean are supported
## Return value: a dissimilarity measure between cluster1 and cluster2
## ========================================================================================
dist.cluster = function(cluster1, cluster2, dist.type = 'Mahalanobis')
{
  if(!is(cluster1,'Cluster') && !is(cluster1,'MetaCluster'))
  {
    stop('cluster1 must be objects of class Cluster or MetaCluster')
  }
  if(!is(cluster2,'Cluster') && !is(cluster2,'MetaCluster'))
  {
    stop('cluster2 must be objects of class Cluster or MetaCluster')
  }
  mu1 = get.center(cluster1)
  sigma1 = get.cov(cluster1)
  n1 = get.size(cluster1)
  
  mu2 = get.center(cluster2)
  sigma2 = get.cov(cluster2)
  n2 = get.size(cluster2)

  if(dist.type=='Mahalanobis')
  {
    dist = mahalanobis.dist(mu1,mu2,sigma1,sigma2,n1,n2);
  } 
  else if(dist.type=='KL')
  {
    dist = symmetric.KL(mu1,mu2,sigma1,sigma2);
  }
  else if(dist.type=='Euclidean') 
  {
    dist = sqrt(sum((mu1-mu2)^2 ))
  }
  else
  {
    stop('Allowable distance types are: Euclidean, Mahalanobis and KL ')
  }
  
  return (dist)
}

###########################################################################################
###########################################################################################

## calculate dissimilarity matrix between two clustered samples or two templates
## or a clustered sample and a template
## (i,j) entry of the matrix is the dissimilarity
## between the ith cluster from sample 1 and the jth cluster from sample 2
## input:
##    object1: An object of class "ClusteredSample" or "Template" 
##    object2: An object of class "ClusteredSample" or "Template" 
##    dist.type: The method with which the dissimilarity is computed between clusters
##                Mahalanobis, KL and Euclidean are supported

###########################################################################################
###########################################################################################

dist.matrix = function(object1, object2, dist.type = 'Mahalanobis')
{
  n1=NULL
  n2=NULL
  if(is(object1, 'ClusteredSample'))
    n1 = get.num.clusters(object1)
  else if(is(object1, 'Template'))
    n1 = get.num.metaclusters(object1)
  else
    stop('object1 must be an object of class ClusteredSample or Template')
  if(is(object2, 'ClusteredSample'))
    n2 = get.num.clusters(object2)
  else if(is(object2, 'Template'))
    n2 = get.num.metaclusters(object2)
  else
    stop('object2 must be an object of class ClusteredSample or Template')

  DM = matrix(0, nrow=n1, ncol=n2)
  for(c1 in 1:n1)
  {
    if(is(object1, 'ClusteredSample'))
      cls1 = get.clusters(object1)[[c1]]
    else if(is(object1, 'Template'))
      cls1 = get.metaClusters(object1)[[c1]]
    for(c2 in 1:n2)
    {
      if(is(object2, 'ClusteredSample'))
        cls2 = get.clusters(object2)[[c2]]
      else if(is(object2, 'Template'))
        cls2 = get.metaClusters(object2)[[c2]]
      DM[c1,c2] = dist.cluster(cls1, cls2, dist.type)
    }
  }
  return (DM);
}



###########################################################################################
###########################################################################################

## calculate dissimilarity between two clustered samples by Mixed edge cover
## input:
##    clustSample1: An object of class "ClusteredSample" 
##    clustSample2: An object of class "ClusteredSample" 
##    dist.type: The method with which the dissimilarity is computed between clusters
##                Mahalanobis, KL and Euclidean are supported

###########################################################################################
###########################################################################################

dist.sample = function(clustSample1, clustSample2, dist.type = 'Mahalanobis', unmatch.penalty=999999)
{
  if(!is(clustSample1,'ClusteredSample'))
  {
    stop('clustSample1 must be an object of class ClusteredSample')
  }
  if(!is(clustSample2,'ClusteredSample') )
  {
    stop('clustSample2 must be an object of class ClusteredSample')
  }
  
  D = match.clusters(clustSample1, clustSample2, dist.type = dist.type, unmatch.penalty=unmatch.penalty)
  
  return (get.matching.cost(D));
}



###########################################################################################
###########################################################################################

## calculate dissimilarity between two clustered samples by Mixed edge cover
## input:
##    template1: An object of class "Template" 
##    template2: An object of class "Template" 
##    dist.type: The method with which the dissimilarity is computed between meta-clusters
##                Mahalanobis, KL and Euclidean are supported

###########################################################################################
###########################################################################################

dist.template = function(template1, template2, dist.type = 'Mahalanobis', unmatch.penalty=999999)
{
  if(!is(template1,'Template'))
  {
    stop('template1 must be an object of class Template')
  }
  if(!is(template2,'Template') )
  {
    stop('template2 must be an object of class Template')
  }
  
  D = match.clusters(template1, template2, dist.type = dist.type, unmatch.penalty=unmatch.penalty)
  
  return (get.matching.cost(D));
}
