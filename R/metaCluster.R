create.template <- function(clustSamples, dist.type = "Mahalanobis", unmatch.penalty=999999, template.id = NA_integer_)
{
  DMs = list();
  
  num.samples = length(clustSamples);
  k=1
  for(i in 1:num.samples)
  {
    clustSample1 = clustSamples[[i]]
    if(!is(clustSample1,'ClusteredSample'))
    {
      stop('Each element of list clustSamples must be an object of class ClusteredSample')
    }
    for(j in i:num.samples)
    {
      clustSample2 = clustSamples[[j]]
      if(!is(clustSample2,'ClusteredSample'))
      {
        stop('Each element of list clustSamples must be an object of class ClusteredSample')
      }
      #DMs[,,j,i]=  dist.matrix(clusteredSamples[[i]], clusteredSamples[[j]], dist.type) ;
      DM = dist.matrix(clustSample1, clustSample2, dist.type) 
      DMs[[k]] = DM
      k = k+1
    }
  }
  
  consistency = 0
  beta=-1
  lambda = unmatch.penalty
  template1 = .Call( "createTemplate", DMs, consistency, lambda, beta, PACKAGE = "flowMatch" )
  #class(template1$tree) = 'hclust';
  #tree.dendrogram <- as.dendrogram(template$tree)
  #template$tree =  tree.dendrogram;
  
  ## now construct an object of class Template
  meta.clusters = list()
  for(i in 1:length(template1$template$metaClusters))
  {
    clusters = list() # list of clusters in a meta-cluster
    mc = template1$template$metaClusters[[i]]
    for(j in 1:length(mc$clusters))
    {
      sample.id = mc$samples[j] + 1 ## indexing issue in c++ and R
      cluster.id = mc$clusters[j] + 1
      clustSample = clustSamples[[sample.id]]
      cluster = get.clusters(clustSample)[[cluster.id]]
      sample.id(cluster) <- as.integer(sample.id) ### set sample id here ...........
      clusters = c(clusters, cluster)
    }
    meta.cluster = MetaCluster(clusters)
    meta.clusters = c(meta.clusters, meta.cluster)
  }
  template = Template(meta.clusters, tree = template1$tree, template.id=template.id)
  return (template);
}


