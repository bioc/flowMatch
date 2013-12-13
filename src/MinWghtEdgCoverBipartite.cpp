#include <iostream>
#include <cstdlib> // for atoi, atof.
#include <cstring> // for strcmp, strcpy.
#include <map>
#include <string>
#include <Rcpp.h>
#include "BipartiteGraph.h"


//#define STAT
#define SUCCESS 1
#define FAILURE 0

using namespace std;
using namespace Rcpp ;

void mergeTemplate(struct templatePair tPair,int t1, int t2,struct classTemplate& mergedTemplate, int treeIndex );
void sampleFiles(std::vector<string>& fileNamesVec);
void sampleFiles1(std::vector<string>& fileNamesVec);

//classTemplate buildTemplate(BipartiteGraph::Weight lambda, int stimulation,string base[] ,string time) ;
double findMinCostPair(vector<vector<struct templatePair> > templatePairs,int* template1, int* template2 );
void rearrange(vector<vector<struct templatePair> >& templatePairs,classTemplate mergedTemplate, int t1, int t2,BipartiteGraph::Weight lambda);
void rearrange_reEstimate(vector<vector<struct templatePair> >& templatePairs,classTemplate mergedTemplate,int t1, int t2,BipartiteGraph::Weight lambda);

void printTemplate(classTemplate ct);
void computeEdgCover(templatePair& newPair,BipartiteGraph::Weight lambda);
BipartiteGraph::Weight metaMetaDist(int mcSource, struct metaCluster mcGroup, int s, int t1, int t2, 
		vector<vector<BipartiteGraph::Weight> > weightMatrixVecVec1,
		vector<vector<BipartiteGraph::Weight> > weightMatrixVecVec2);

BipartiteGraph::Weight metaMetaDist_basic(metaCluster mc1, metaCluster mc2, vector<vector<BipartiteGraph> > bGraphs);


void initAlign_reEstimate(int numSample, vector<vector<BipartiteGraph> >& bGraphs, vector<vector<struct templatePair> >& templatePairs, double lambda, SEXP tmixes);
void initAlign(int numSample, vector<vector<BipartiteGraph> >& bGraphs, vector<vector<struct templatePair> >& templatePairs, double lambda);

classTemplate buildTemplate(vector<vector<struct templatePair> >& templatePairs, double lambda, struct mclust& tree);
classTemplate buildTemplate_reEstimate(vector<vector<struct templatePair> >& templatePairs, double lambda, struct mclust& tree, SEXP samples, SEXP tmixes, double fraction);


double computeDegConsistency(vector<vector<struct templatePair> >& templatePairs, vector<vector<BipartiteGraph> >& bGraphs, double beta);
void computeDegConsistencyWeighted(vector<vector<struct templatePair> >& templatePairs, vector<vector<BipartiteGraph> >& bGraphs, double beta);
void computeDegConsistencyNew(vector<vector<struct templatePair> >& templatePairs, vector<vector<BipartiteGraph> >& bGraphs, double beta, double lambda);
double samples2TemplateCost(vector<vector<struct templatePair> >& templatePairs, vector<vector<BipartiteGraph> >& bGraphs, classTemplate rootTemplate, BipartiteGraph::Weight lambda);
void mySort(vector<int>& samples, vector<int>& clusters);
void update_height(mclust& tree, vector<vector<struct templatePair> > tPairs, vector<vector<BipartiteGraph> > bGraphs, double lambda);


struct metaCluster
{
	std::vector<int> cluster; // history
	std::vector<int> originSample; //history
	std::vector<int> metaCluster; // immediate  
	std::vector<int> originTemplate; // immediate
	std::vector<int> metaClusterSize; //immediate

	//tmix parameters
	std::vector<double> mean;
	std::vector<std::vector<double> > cov;
	double nu;

};

struct classTemplate
{
	//int sampleCount;
	std::vector<int> samplesVec; // history ... list the samples in the template
	int myTreeIndex; // singleton template always get -negative values
	// for non-sigleton template these values are useful
	int leftChildTreeIndex;
	int rightChildTreeIndex;
	//int metaClusterCount;
	std::vector<struct metaCluster> metaClustersVec;
	Rcpp::List tmix; // saves the tmix associated with this template (only used with the reastimated version)
};

struct templatePair
{
	struct classTemplate template1;
	struct classTemplate template2;
	std::vector<std::vector<BipartiteGraph::Vertex> > sCoverVecVec;
	std::vector<std::vector<BipartiteGraph::Vertex> > tCoverVecVec;
	BipartiteGraph::Weight coverWeight;
	vector<vector<BipartiteGraph::Weight> > weightMatrixVecVec; // save the weight matrix for this pair
	// only used when consistency option is used
	vector<vector<BipartiteGraph::Weight> > weightMatrixVecVec_original; // also keep the weight matrix in the original scale if consistency is used
	// only used when used consistency+reestimate together
	vector<vector<double> > alpha; // save alpha for this pair, needs to update for reestimated version

};

struct treeNode
{
	struct classTemplate cTemplate;
	struct treeNode* left;
	struct treeNode* right;
	double leftDist;
	double rightDist;
};

struct mclust
{
	//an n-1 by 2 matrix. Row i of merge describes the merging of clusters at step i of the clustering.
	//If an element j in the row is negative, then observation -j was merged at this stage.
	//If j is positive then the merge was with the cluster formed at the (earlier) stage j of the algorithm.
	//Thus negative entries in merge indicate agglomerations of singletons, and positive entries indicate agglomerations of non-singletons.
	std::vector<std::vector<int> > merge;
	std::vector<double> height;
	std::vector<struct classTemplate> templates;
};



void cover(vector<vector<BipartiteGraph> > & bGraphs);
classTemplate getTemplate_reEstimate(vector<vector<BipartiteGraph> > & bGraphs, mclust& tree, int consistency, SEXP samples, SEXP tmixes, double fraction, double lambda);
classTemplate getTemplate(vector<vector<BipartiteGraph> > & bGraphs, mclust& tree, int consistency, double lambda, double beta);
Rcpp::List createResult(classTemplate ct, mclust tree);
Rcpp::List listify_template(classTemplate ct);




//##### for testing
vector<vector<BipartiteGraph> > bGraphs_global;


/******************************************************************************
 ** This function creates a template from samples already clustered
 ** ( no estimation of parameters in the internal nodes)
 ** Parameters:
 ** 	dist: a list of distances for all pair of samples
 ** 	con : Will consistency be used. 0: do not use , 1: unweighted, 2: weighted
 ******************************************************************************/
RcppExport SEXP createTemplate(SEXP dist, SEXP con, SEXP lda, SEXP b)
{

	int consistency = as<int>(con);
	double lambda = as<double>(lda);
	double beta = as<double>(b);
	Rcpp::List lDist (dist);
	int numSamples = (-1 + sqrt(1+8 * lDist.length()))/2;
	// Rcpp::Rcout<< "Number of samples: "<< numSamples << endl;
	//Rcpp::Rcout<< "number of elements = " << lDist.length() << endl;
	int curMatrix = 0;

	vector<vector<BipartiteGraph> > bGraphs;
	bGraphs.resize(numSamples);
	for(int i=0; i<numSamples; i++)
	{
		bGraphs[i].resize(numSamples);
		for(int j=i; j<numSamples; j++)
		{
			SEXP tmp = lDist[curMatrix++];
			Rcpp::NumericMatrix nmDist(tmp);
			int nrows = nmDist.nrow();
			int ncolumns = nmDist.ncol();
			// create a 2D vector to store a distance matrix
			vector<vector<double> >  distVecVec;
			distVecVec.resize(nrows);
			for(int m=0; m<nrows; m++)
			{
				distVecVec[m].resize(ncolumns);
				for(int n=0; n<ncolumns; n++)
				{
					distVecVec[m][n] = nmDist(m,n);
				}
			}

			// now create a bipartite graph from the matrix
			BipartiteGraph graph(distVecVec);

			bGraphs[i][j] = graph;
		}
	}

	bGraphs_global = bGraphs;   /// ##############

	mclust tree;
	classTemplate ct = getTemplate(bGraphs,tree, consistency, lambda, beta);
	return(Rcpp::wrap(createResult(ct,tree)));
	// prepare a wrapper for the return template
	 //Rcpp::List a(1);
	 //return(a);
}



Rcpp::List listify_hclust(mclust tree)
{
	Rcpp::NumericVector height( tree.height.begin(), tree.height.end());
	int numSamples = tree.height.size() + 1;
	Rcpp::NumericMatrix merge(numSamples-1, 2);
	Rcpp::NumericVector labels(numSamples);
	Rcpp::NumericVector order(numSamples);

	for(int i=0; i<(numSamples-1); i++)
	{
		merge(i,0) = tree.merge[i][0];
		merge(i,1) = tree.merge[i][1];
		labels(i) = i+1;
		order(i) = i+1;
	}
	order(numSamples-1) = numSamples;
	labels(numSamples-1) = numSamples;

	return Rcpp::List::create( Rcpp::Named("height") = height ,
			Rcpp::Named("merge") = merge,
			Rcpp::Named("labels") = labels,
			Rcpp::Named("order") = order);

}

Rcpp::List listify_template(classTemplate ct)
{

	int numMetaCluster = ct.metaClustersVec.size();
	//Rcpp::Rcout<< "Number of Meta Cluster = " << numMetaCluster << endl;
	Rcpp::List allmc(numMetaCluster); // it worked
	for(int i=0; i<numMetaCluster; i++)
	{
		metaCluster mc = ct.metaClustersVec[i];
		// why not working ???
		/*
		for(i=0; i<mc.cluster.size(); i++)
		{
			mc.cluster[i] = mc.cluster[i] + 1;
			mc.originSample[i] = mc.originSample[i] + 1;
		}*/
		mySort(mc.originSample, mc.cluster);
		// make it 1-base indexing for R

		Rcpp::NumericVector nvCluster( mc.cluster.begin(), mc.cluster.end()); // history
		Rcpp::NumericVector nvOriginSample(mc.originSample.begin(),mc.originSample.end()); //history
		Rcpp::NumericVector nvMetaCluster (mc.metaCluster.begin(), mc.metaCluster.end()); // immediate
		Rcpp::NumericVector nvOriginTemplate(mc.originTemplate.begin(), mc.originTemplate.end()); // immediate
		Rcpp::NumericVector nvMetaClusterSize(mc.metaClusterSize.begin(), mc.metaClusterSize.end()); //immediate

		// create the meta-cluster as a list
		Rcpp::List mcl = Rcpp::List::create( Rcpp::Named("clusters") = nvCluster , // original clsuter id
				Rcpp::Named("samples") = nvOriginSample, // original sample
				Rcpp::Named("metaCluster") = nvMetaCluster,
				Rcpp::Named("originTemplate") = nvOriginTemplate,
				Rcpp::Named("metaClusterSize") = nvMetaClusterSize);

		allmc[i] = mcl;

	}

	//int sampleCount;
	Rcpp::NumericVector nvSamplesVec(ct.samplesVec.begin(), ct.samplesVec.end()); // history ... list the samples in the template
	//Rcpp::NumericVector nv(ct.templateVec.begin(), ct.templateVec.end());  // immediate

	return Rcpp::List::create(Rcpp::Named("samples") = nvSamplesVec,
			Rcpp::Named("metaClusters") = allmc);

}


Rcpp::List createResult(classTemplate ct, mclust tree)
{
	Rcpp::List t = listify_template(ct);
	Rcpp::List tr = listify_hclust(tree);
	return Rcpp::List::create(Rcpp::Named("template") = t,
			Rcpp::Named("tree") = tr);
}


void update_height(mclust& tree, vector<vector<struct templatePair> > tPairs, vector<vector<BipartiteGraph> > bGraphs, double lambda)
{
	for(int i=0; i<tree.height.size(); i++)
	{
		int left = tree.merge[i][0];
		int right = tree.merge[i][1];
		classTemplate t1;
		classTemplate t2;
		if(left >= 0)
			t1 = tree.templates[left-1];
		else
			t1 = tPairs[-left-1][-left-1].template1;
		if(right >= 0)
			t2 = tree.templates[right-1];
		else
			t2 = tPairs[-right-1][-right-1].template1;


		int sNumVertices = t1.metaClustersVec.size();
		int tNumVertices = t2.metaClustersVec.size();

		//Rcpp::Rcout<< sNumVertices << "  *** " << tNumVertices << endl;
		vector<vector<BipartiteGraph::Weight> > newWeightMatrixVecVec;
		newWeightMatrixVecVec.reserve(sNumVertices);
		newWeightMatrixVecVec.resize(sNumVertices);

		//for every meta-cluster of ith template
		for(int j=0; j<sNumVertices; j++)
		{
			newWeightMatrixVecVec[j].reserve(tNumVertices);
			newWeightMatrixVecVec[j].resize(tNumVertices);
			// for every meta-cluster of new template
			for(int k=0; k<tNumVertices; k++)
			{
				newWeightMatrixVecVec[j][k] = metaMetaDist_basic(t1.metaClustersVec[j], t2.metaClustersVec[k], bGraphs);
			}
		}

		BipartiteGraph bgraph(newWeightMatrixVecVec);
		std::vector<std::vector<int> > sCoverVecVec;
		std::vector<std::vector<int> > tCoverVecVec;
		double coverWeight;
		bgraph.MinWghtGenEdgCover(sCoverVecVec,tCoverVecVec,&coverWeight,lambda);
		tree.height[i] = coverWeight/2;
	}

}

/*
 *  Create a template with or without the consistency
 *  consistency = 0 : no concictency
 *  consistency = 1 : unweighted consistency
 *  consistency = 2 : weighted consistency
 */
classTemplate getTemplate(vector<vector<BipartiteGraph> > & bGraphs, mclust& tree, int consistency, double lambda, double beta)
{
	//BipartiteGraph::Weight lambda = 3.0;
	double avgSample2TemplateCost;
	int numSample = bGraphs.size();
	classTemplate rootTemplate;

	vector<vector<struct templatePair> > templatePairs;
	initAlign(numSample, bGraphs, templatePairs, lambda);

	if(consistency == 0)
	{
		//Rcpp::Rcout<< "*********** Without consistency **********\n";

		vector<vector<struct templatePair> > templatePairs1 = templatePairs; // do not change the original one
		rootTemplate = buildTemplate(templatePairs1, lambda, tree);
		update_height(tree, templatePairs, bGraphs, lambda);
		printTemplate(rootTemplate);
	}

	else if(consistency == 1)
	{
		//Rcpp::Rcout<< "*********** with consistency **********\n";
		vector<vector<BipartiteGraph> > bGraphs_updated = bGraphs;
		//double beta = 0;
		if (beta<0)
			beta = (numSample-2)/(double)numSample;
		//beta = .5;
		// Rcpp::Rcout<< "Beta = " << beta << endl;
		computeDegConsistency(templatePairs, bGraphs_updated, beta);
		//lambda = 1;//lambda * (1-alpha_med);
		//Rcpp::Rcout<< "The rescaled lambda is: " << lambda << endl;
		//computeDegConsistencyWeighted(templatePairs, bGraphs_updated, beta);  /// check for problem
		vector<vector<struct templatePair> > templatePairs_updated;
		initAlign(numSample, bGraphs_updated, templatePairs_updated, lambda);
		rootTemplate = buildTemplate(templatePairs_updated, lambda, tree);
		//update_height(tree, templatePairs, bGraphs, lambda);

		printTemplate(rootTemplate);
		
	}

	else if(consistency == 2)
	{
		//Rcpp::Rcout<< "*********** with consistency **********\n";
		vector<vector<BipartiteGraph> > bGraphs_updated = bGraphs;
		if (beta<0)
			beta = (numSample-2)/(double)numSample;

		//computeDegConsistency(templatePairs, bGraphs_updated, beta);
		computeDegConsistencyWeighted(templatePairs, bGraphs_updated, beta);  /// check for problem
		vector<vector<struct templatePair> > templatePairs_updated;
		initAlign(numSample, bGraphs_updated, templatePairs_updated, lambda);
		rootTemplate = buildTemplate(templatePairs_updated, lambda, tree);
		//update_height(tree, templatePairs, bGraphs, lambda);

		printTemplate(rootTemplate);
		
	}
	else if(consistency == 3)
		{
			//Rcpp::Rcout<< "*********** with new consistency **********\n";
			vector<vector<BipartiteGraph> > bGraphs_updated = bGraphs;
			//double beta = 0;
			if (beta<0)
				beta = (numSample-2)/(double)numSample;
			// Rcpp::Rcout<< "Beta = " << beta << endl;
			//computeDegConsistency(templatePairs, bGraphs_updated, beta);
			computeDegConsistencyNew(templatePairs, bGraphs_updated, beta, lambda);  /// check for problem
			vector<vector<struct templatePair> > templatePairs_updated;
			initAlign(numSample, bGraphs_updated, templatePairs_updated, lambda);
			rootTemplate = buildTemplate(templatePairs_updated, lambda, tree);
			update_height(tree, templatePairs, bGraphs, lambda);

			printTemplate(rootTemplate);
			
		}
	else
	{
		// Rcpp::Rcout<< "Invalid consistency option .... returning \n";
	}


	return rootTemplate;
}




void initAlign(int numSample, vector<vector<BipartiteGraph> >& bGraphs, vector<vector<struct templatePair> >& templatePairs, double lambda)
{
	templatePairs.resize(numSample);
	for(int i=0; i<numSample; i++)
	{
		templatePairs[i].resize(numSample);
		for(int j=i; j<numSample; j++)
		{
			BipartiteGraph graph = bGraphs[i][j];
			classTemplate t1;
			t1.samplesVec.push_back(i);
			t1.myTreeIndex = -1 * (i+1); // singleton template always get -1

			int metaClusterCount1 = graph.GetSNumVertices();
			for(int k=0; k< metaClusterCount1; k++)
			{
				struct metaCluster mc;
				mc.cluster.push_back(k); // initially every meta-cluster has only one cluster
				mc.metaCluster.push_back(k);
				mc.originSample.push_back(i);
				mc.originTemplate.push_back(i);
				mc.metaClusterSize.push_back(1);
				t1.metaClustersVec.push_back(mc);
			}

			classTemplate t2;
			t2.samplesVec.push_back(j);
			t2.myTreeIndex = -1 * (j+1); // singleton template always get -1

			int metaClusterCount2 = graph.GetTNumVertices();
			for(int k=0; k<metaClusterCount2; k++)
			{
				struct metaCluster mc;
				mc.cluster.push_back(k);
				mc.metaCluster.push_back(k);
				mc.originSample.push_back(j);
				mc.originTemplate.push_back(j);
				mc.metaClusterSize.push_back(1);
				t2.metaClustersVec.push_back(mc);
			}


			templatePair tPair;
			tPair.template1 = t1;
			tPair.template2 = t2;
			tPair.weightMatrixVecVec = graph.GetSWeightMatrix(); // assumes graph edge weights are sorted based on adjacency

			graph.MinWghtGenEdgCover(tPair.sCoverVecVec,tPair.tCoverVecVec,&tPair.coverWeight,lambda);

			//printf("%d %d %lf\n",i,j,tPair.coverWeight);

			templatePairs[i][j] = tPair;
			//templatePairs[j][i] = tPair; // extra storage for ease of computation of degree of consistency
		}
	}
}


/*
 * Compute the degree of consistency for each pair of clusters weighted
 */
void computeDegConsistencyNew(vector<vector<struct templatePair> >& templatePairs, vector<vector<BipartiteGraph> >& bGraphs, double beta, double lambda)
{
	int numSample = templatePairs.size();
	for(int x=0; x<numSample; x++)
	{
		for(int y=x+1; y<numSample; y++)
		{
			BipartiteGraph graph_xy = bGraphs[x][y];
			for(int cx=0; cx< graph_xy.GetSNumVertices(); cx++)
			{
				vector<int> cover_cx_y = templatePairs[x][y].sCoverVecVec[cx];
				for(int cy=0; cy< graph_xy.GetTNumVertices(); cy++)
				{
					//------------- get direct matching delta(cx,cy) ---------------
					int delta_cx_cy  = 0;
					for(unsigned int i=0; i<cover_cx_y.size(); i++)
					{
						if(cover_cx_y[i] == cy)
						{
							delta_cx_cy = 1;
							break;
						}
					}
					// ----------- direct matchign calculation finished ---------------

					//------------ compute degree of consistency d(cx,cy) -------------
					//int d_cx_cy = 0;
					double w_cx_cy = 0;
					//double t_cx_cy = 0;

					for(int z=0; z<numSample; z++)
					{
						if(z!=x && z!=y)
						{
							int dz_cx_cy = 0;
							//double wz_cx_cy = 0;
							//double tz_cx_cy = 0;
							vector<int> cover_cx_z;
							vector<int> cover_cy_z;
							vector<double> weight_cx_z;
							vector<double> weight_cy_z;

							if(x<z) cover_cx_z = templatePairs[x][z].sCoverVecVec[cx];
							else cover_cx_z = templatePairs[z][x].tCoverVecVec[cx];
							if(y<z) cover_cy_z = templatePairs[y][z].sCoverVecVec[cy];
							else cover_cy_z = templatePairs[z][y].tCoverVecVec[cy];

							if(x<z) weight_cx_z = bGraphs[x][z].mSWeightVecVec[cx];
							else weight_cx_z = bGraphs[z][x].mTWeightVecVec[cx];
							if(y<z) weight_cy_z = bGraphs[y][z].mSWeightVecVec[cy];
							else weight_cy_z = bGraphs[z][y].mTWeightVecVec[cy];

							// find common neighbors
							double max_weight = 0;
							for(unsigned int i=0; i<cover_cx_z.size(); i++)
							{
								for(unsigned int j=0; j<cover_cy_z.size(); j++)
								{
									if(cover_cx_z[i] == cover_cy_z[j])
									{
										int cz = cover_cx_z[i];
										dz_cx_cy++;
										//wz_cx_cy += weight_cx_z[cz] + weight_cy_z[cz];
										if(max_weight < weight_cx_z[cz])
											max_weight = weight_cx_z[cz];
										if(max_weight < weight_cy_z[cz])
											max_weight = weight_cy_z[cz];
									}
								}
							}
							if(dz_cx_cy == 0)
							{
								max_weight= lambda;

							}

							w_cx_cy += max_weight;
						}
					}
					//------------ degree of consistency calculation finished-------------


					double w = graph_xy.mSWeightVecVec[cx][cy]; // the indexing is possible if graph is compplete and edges are sorted
					double w_new =  beta * w_cx_cy/(numSample-2) +  (1-beta) * w   ;

          // save the updated weight
					bGraphs[x][y].mSWeightVecVec[cx][cy] = w_new;
					bGraphs[x][y].mTWeightVecVec[cy][cx] = w_new;
					//printf("[%d][%d][%d][%d] direct matchign = %d, indirect matchign = %d, alpha=%lf\n",x,y,cx,cy,delta_cx_cy,d_cx_cy,alpha_cx_cy);
				}
			}
		}
	}
}



/*
 * Compute the degree of consistency for each pair of clusters unweighted
 */
double computeDegConsistency(vector<vector<struct templatePair> >& templatePairs, vector<vector<BipartiteGraph> >& bGraphs, double beta)
{
	int numSample = templatePairs.size();
	//vector<double> alpha;
	for(int x=0; x<numSample; x++)
	{
		for(int y=x+1; y<numSample; y++)
		{
			BipartiteGraph graph_xy = bGraphs[x][y];
			for(int cx=0; cx< graph_xy.GetSNumVertices(); cx++)
			{
				vector<int> cover_cx_y = templatePairs[x][y].sCoverVecVec[cx];
				for(int cy=0; cy< graph_xy.GetTNumVertices(); cy++)
				{
					//------------- get direct matching delta(cx,cy) ---------------
					int delta_cx_cy  = 0;
					for(unsigned int i=0; i<cover_cx_y.size(); i++)
					{
						if(cover_cx_y[i] == cy)
						{
							delta_cx_cy = 1;
							break;
						}
					}
					// ----------- direct matchign calculation finished ---------------

					//------------ compute degree of consistency d(cx,cy) -------------
					int d_cx_cy = 0;
					for(int z=0; z<numSample; z++)
					{
						if(z!=x && z!=y)
						{
							int dz_cx_cy = 0;
							vector<int> cover_cx_z;
							vector<int> cover_cy_z;
							if(x<z) cover_cx_z = templatePairs[x][z].sCoverVecVec[cx];
							else cover_cx_z = templatePairs[z][x].tCoverVecVec[cx];
							if(y<z) cover_cy_z = templatePairs[y][z].sCoverVecVec[cy];
							else cover_cy_z = templatePairs[z][y].tCoverVecVec[cy];
							// find common neighbors
							for(unsigned int i=0; i<cover_cx_z.size() && dz_cx_cy==0; i++)
							{
								for(unsigned int j=0; j<cover_cy_z.size() && dz_cx_cy==0 ; j++)
								{
									if(cover_cx_z[i] == cover_cy_z[j])
									{
										d_cx_cy++;
										dz_cx_cy = 1;
									}
								}
							}
						}
					}
					//------------ degree of consistency calculation finished-------------

					// now calculate alpha and update weight
					double alpha_cx_cy = beta * (double)d_cx_cy/(numSample-2) + (1-beta) * delta_cx_cy;
					double w = graph_xy.mSWeightVecVec[cx][cy]; // the indexing is possible if graph is compplete and edges are sorted
					double w_new = w * (1-alpha_cx_cy);
					
					// save the updated weight
					bGraphs[x][y].mSWeightVecVec[cx][cy] = w_new;
					bGraphs[x][y].mTWeightVecVec[cy][cx] = w_new;
					//printf("[%d][%d][%d][%d] direct matchign = %d, indirect matchign = %d, alpha=%lf\n",x,y,cx,cy,delta_cx_cy,d_cx_cy,alpha_cx_cy);
				}
			}
		}
	}
	return 1.0; // just dummy return
}



/*
 * Compute the degree of consistency for each pair of clusters weighted
 */
void computeDegConsistencyWeighted(vector<vector<struct templatePair> >& templatePairs, vector<vector<BipartiteGraph> >& bGraphs, double beta)
{
	int numSample = templatePairs.size();
	for(int x=0; x<numSample; x++)
	{
		for(int y=x+1; y<numSample; y++)
		{
			BipartiteGraph graph_xy = bGraphs[x][y];
			for(int cx=0; cx< graph_xy.GetSNumVertices(); cx++)
			{
				vector<int> cover_cx_y = templatePairs[x][y].sCoverVecVec[cx];
				for(int cy=0; cy< graph_xy.GetTNumVertices(); cy++)
				{
					//------------- get direct matching delta(cx,cy) ---------------
					int delta_cx_cy  = 0;
					for(unsigned int i=0; i<cover_cx_y.size(); i++)
					{
						if(cover_cx_y[i] == cy)
						{
							delta_cx_cy = 1;
							break;
						}
					}
					// ----------- direct matchign calculation finished ---------------

					//------------ compute degree of consistency d(cx,cy) -------------
					//int d_cx_cy = 0;
					double w_cx_cy = 0;
					double t_cx_cy = 0;

					for(int z=0; z<numSample; z++)
					{
						if(z!=x && z!=y)
						{
							int dz_cx_cy = 0;
							double wz_cx_cy = 0;
							double tz_cx_cy = 0;
							vector<int> cover_cx_z;
							vector<int> cover_cy_z;
							vector<double> weight_cx_z;
							vector<double> weight_cy_z;

							if(x<z) cover_cx_z = templatePairs[x][z].sCoverVecVec[cx];
							else cover_cx_z = templatePairs[z][x].tCoverVecVec[cx];
							if(y<z) cover_cy_z = templatePairs[y][z].sCoverVecVec[cy];
							else cover_cy_z = templatePairs[z][y].tCoverVecVec[cy];

							if(x<z) weight_cx_z = bGraphs[x][z].mSWeightVecVec[cx];
							else weight_cx_z = bGraphs[z][x].mTWeightVecVec[cx];
							if(y<z) weight_cy_z = bGraphs[y][z].mSWeightVecVec[cy];
							else weight_cy_z = bGraphs[z][y].mTWeightVecVec[cy];

							// find common neighbors
							for(unsigned int i=0; i<cover_cx_z.size(); i++)
							{
								for(unsigned int j=0; j<cover_cy_z.size(); j++)
								{
									if(cover_cx_z[i] == cover_cy_z[j])
									{
										int cz = cover_cx_z[i];
										dz_cx_cy++;
										wz_cx_cy += weight_cx_z[cz] + weight_cy_z[cz];
									}
								}
							}
							if(dz_cx_cy > 0)
							{
								wz_cx_cy /= dz_cx_cy;
								wz_cx_cy /= 2;
							}
							if(wz_cx_cy > 0)
							{
								tz_cx_cy = wz_cx_cy;
							}
							else //compute tz
							{
								double min_tz = 9999999;
								if(weight_cx_z.size() != weight_cy_z.size())
								{
									//Rcpp::Rcout<< "The size of weight_cx_z and weight_cy_z do not match.. quiting \n";
									continue;
								}
								for(unsigned int cz=0; cz < weight_cx_z.size() ; cz++)
								{
									double tz = weight_cx_z[cz] + weight_cy_z[cz];
									if(tz < min_tz) min_tz = tz;
								}
								tz_cx_cy = min_tz/2;
							}

							t_cx_cy += tz_cx_cy;
							w_cx_cy += wz_cx_cy;
						}
					}
					//------------ degree of consistency calculation finished-------------

					// now calculate alpha and update weight
					double alpha_cx_cy;
					if(t_cx_cy < .00001) // stop unstable division
						alpha_cx_cy = (1-beta) * delta_cx_cy;
					else
						alpha_cx_cy = beta * (double)w_cx_cy/t_cx_cy + (1-beta) * delta_cx_cy;
					double w = graph_xy.mSWeightVecVec[cx][cy]; // the indexing is possible if graph is compplete and edges are sorted
					double w_new = w * (1-alpha_cx_cy);

          // save the updated weight
					bGraphs[x][y].mSWeightVecVec[cx][cy] = w_new;
					bGraphs[x][y].mTWeightVecVec[cy][cx] = w_new;
					//printf("[%d][%d][%d][%d] direct matchign = %d, indirect matchign = %d, alpha=%lf\n",x,y,cx,cy,delta_cx_cy,d_cx_cy,alpha_cx_cy);
				}
			}
		}
	}
}


/*
 * UPGMA type tree building to get templates
 */
classTemplate buildTemplate(vector<vector<struct templatePair> >& templatePairs, double lambda, mclust& tree)
{

	classTemplate rootTemplate;
	int numSamples = templatePairs.size();
	int curIndex = 0;


	tree.merge.resize(numSamples-1);
	for(int i=0; i< (numSamples-1); i++)
	{
		tree.merge[i].resize(2);
	}
	tree.height.resize(numSamples-1);
	tree.templates.resize(numSamples-1);


	while(templatePairs.size() > 1)
	{
		int template1, template2;
		double minCost = findMinCostPair(templatePairs, &template1, &template2 );


		classTemplate mergedTemplate;
		mergeTemplate(templatePairs[template1][template2], template1, template2,mergedTemplate, curIndex+1 );

		tree.merge[curIndex][0] = mergedTemplate.leftChildTreeIndex;
		tree.merge[curIndex][1] = mergedTemplate.rightChildTreeIndex;
		tree.height[curIndex] = minCost/2;
		tree.templates[curIndex++] = mergedTemplate;

		if(templatePairs.size() == 2)
		{
			return mergedTemplate;
		}

		rearrange(templatePairs, mergedTemplate, template1, template2,lambda);
	}

	return rootTemplate;
}





void rearrange(vector<vector<struct templatePair> >& templatePairs,classTemplate mergedTemplate,int t1, int t2,BipartiteGraph::Weight lambda)
{
	//assert(t1 < t2);
	// Rcpp::Rcout<< t1 << "+" << t2 << endl;
	// now compute distance from all meta-clusters of new template to others
	int templateCount = templatePairs.size();


	for(int i=0; i<templateCount; i++)
	{
		templatePair newPair;
		newPair.template2 = mergedTemplate;
		templatePair oldPair1,oldPair2;
		if(i == t1 || i == t2) continue;
		if(i < t1)
		{
			oldPair1 = templatePairs[i][t1];
			newPair.template1 = oldPair1.template1;
		}
		else
		{
			oldPair1 = templatePairs[t1][i];
			newPair.template1 = oldPair1.template2;
		}
		if(i < t2)
			oldPair2 = templatePairs[i][t2];
		else
			oldPair2 = templatePairs[t2][i];


		int sNumVertices = newPair.template1.metaClustersVec.size();
		int tNumVertices = newPair.template2.metaClustersVec.size();

		//Rcpp::Rcout<< sNumVertices << "  *** " << tNumVertices << endl;
		vector<vector<BipartiteGraph::Weight> > newWeightMatrixVecVec;
		newWeightMatrixVecVec.reserve(sNumVertices);
		newWeightMatrixVecVec.resize(sNumVertices);

		//for every meta-cluster of ith template
		for(int j=0; j<sNumVertices; j++)
		{
			newWeightMatrixVecVec[j].reserve(tNumVertices);
			newWeightMatrixVecVec[j].resize(tNumVertices);
			// notice, we do not need the content of this meta-cluster


			// for every meta-cluster of new template
			for(int k=0; k<tNumVertices; k++)
			{
				newWeightMatrixVecVec[j][k] = metaMetaDist(j, mergedTemplate.metaClustersVec[k], i, t1, t2, oldPair1.weightMatrixVecVec, oldPair2.weightMatrixVecVec);

				//double temp = metaMetaDist1(newPair.template1.metaClustersVec[j], mergedTemplate.metaClustersVec[k]);
				//if(temp != newWeightMatrixVecVec[j][k])
				//	Rcpp::Rcout<< "not matched !!" << newWeightMatrixVecVec[j][k] << " & " << temp << endl;

			}
		}

		// erase t1 and t2 columns from ith row
		templatePairs[i].erase(templatePairs[i].begin() + t1);
		templatePairs[i].erase(templatePairs[i].begin() + t2 -1); // -1 since t1 is already deleated
		newPair.weightMatrixVecVec = newWeightMatrixVecVec;
		computeEdgCover(newPair,lambda);
		templatePairs[i].push_back(newPair);
	}

	// remove the row corresponding t1 and t2
	templatePairs.erase(templatePairs.begin() + t1);
	templatePairs.erase(templatePairs.begin() + t2 -1);

	int newSize = templatePairs.size();
	templatePairs.reserve(newSize+1);
	templatePairs.resize(newSize+1);
	templatePairs[newSize].reserve(newSize+1);
	templatePairs[newSize].resize(newSize+1);
}



/* 
calculate distance from a meta-clsuter to a group of meta-clusters
The group of meta cluster comes from either template t1 or template t2

@@@@@ I have verified that this version and metaMetaDist1 calculate same distance @@@@@@

We used UPGMA formula - 7.6 from the red book
Parameters:

s - index of source template
t1, t2 - index of other two templates that are newly merged
mcSource - index of metacluster from template s
mcGroup - a metacluster from recently formed from template t1 and t2
weightMatrixVecVec1 - weight matrix of s-t1 pair
weightMatrixVecVec2 - weight matrix of s-t2 pair
 */
BipartiteGraph::Weight metaMetaDist(int mcSource, metaCluster mcGroup, int s, int t1, int t2, 
		vector<vector<BipartiteGraph::Weight> > weightMatrixVecVec1,
		vector<vector<BipartiteGraph::Weight> > weightMatrixVecVec2)
{

	int mcGroupSize = mcGroup.metaCluster.size();//immediate size;
	int clusterCount = 0;
	BipartiteGraph::Weight w = 0;
	// for each meta-cluster in the group
	for(int k=0; k<mcGroupSize; k++)
	{
		int mcIndex = mcGroup.metaCluster[k];
		int tIndex = mcGroup.originTemplate[k]; // the value is either t1 or t2
		int metaClusterSize = mcGroup.metaClusterSize[k];
		clusterCount +=metaClusterSize;
		//Rcpp::Rcout<< clusterCount << "...\n";
		if (tIndex == t1) // use weightMatrixVecVec1
		{
			//decide which indexing of weightMatrixVecVec1 will be used 
			if(s<t1)
			{
				w += metaClusterSize * weightMatrixVecVec1[mcSource][mcIndex];
			}
			else
			{
				w += metaClusterSize * weightMatrixVecVec1[mcIndex][mcSource];

			}

		}
		else if(tIndex == t2) //use weightMatrixVecVec2
		{
			if(s<t2)
			{
				w += metaClusterSize * weightMatrixVecVec2[mcSource][mcIndex];
			}
			else
			{
				w += metaClusterSize * weightMatrixVecVec2[mcIndex][mcSource];

			}

		}
		else
		{
			return(999999);
		}

	}
	if((unsigned int)clusterCount != mcGroup.cluster.size())
	{
		Rcpp::Rcout<< "error in metaMetaDist(): " << clusterCount << " " << mcGroup.cluster.size() << endl;
	}
	return (w/clusterCount) ;
}



BipartiteGraph::Weight metaMetaDist_basic(metaCluster mc1, metaCluster mc2, vector<vector<BipartiteGraph> > bGraphs)
{
	int mc1Size = mc1.originSample.size();
	int mc2Size = mc2.originSample.size();

	double dist = 0;
	for(int i=0; i<mc1Size; i++)
	{
		int s1 = mc1.originSample[i];
		int c1 = mc1.cluster[i];
		for(int j=0; j<mc2Size; j++)
		{
			int s2 = mc2.originSample[j];
			int c2 = mc2.cluster[j];
			//Rcpp::Rcout<<"(" << s1 << ","<< c1 << " & "<< s2 << ","<< c2 << ") \n";
			if(s1<s2)
			{
				dist += bGraphs[s1][s2].mSWeightVecVec[c1][c2];
				// If you want to use this version please pass bGraphs to this function
			}
			else
			{
				dist += bGraphs[s2][s1].mSWeightVecVec[c2][c1];
			}
		}
	}

	dist = dist/(mc1Size * mc2Size);

	return dist;
}


void computeEdgCover(templatePair& newPair,BipartiteGraph::Weight lambda)
{

	BipartiteGraph graph(newPair.weightMatrixVecVec);
	graph.MinWghtGenEdgCover(newPair.sCoverVecVec,newPair.tCoverVecVec,&newPair.coverWeight,lambda);

}



void mergeTemplate(templatePair tPair,int t1, int t2,classTemplate& mergedTemplate, int myTreeIndex )
{

	vector<vector<BipartiteGraph::Vertex> > sCoverVecVec = tPair.sCoverVecVec;
	vector<vector<BipartiteGraph::Vertex> > tCoverVecVec = tPair.tCoverVecVec;
	classTemplate template1 =  	tPair.template1;
	classTemplate template2 =  	tPair.template2;

	//printTemplate(template1);
	//printTemplate(template2);

	//classTemplate mergedTemplate;
	mergedTemplate.samplesVec = template1.samplesVec;
	mergedTemplate.samplesVec.insert(mergedTemplate.samplesVec.end(),template2.samplesVec.begin(), template2.samplesVec.end());

	mergedTemplate.myTreeIndex = myTreeIndex;
	mergedTemplate.leftChildTreeIndex = template1.myTreeIndex;
	mergedTemplate.rightChildTreeIndex = template2.myTreeIndex;
	//Rcpp::Rcout<< "start merging S side " << endl;
	// form group by merging the stars
	for(unsigned int i=0; i<sCoverVecVec.size(); i++)
	{

		metaCluster mc;
		int mateCount = sCoverVecVec[i].size();

		mc.metaCluster.push_back(i);
		mc.originTemplate.push_back(t1);
		mc.metaClusterSize.push_back(template1.metaClustersVec[i].cluster.size());
		mc.cluster = template1.metaClustersVec[i].cluster;
		mc.originSample = template1.metaClustersVec[i].originSample;

		for(int j=0; j<mateCount; j++)
		{
			int mate = sCoverVecVec[i][j];
			mc.metaCluster.push_back(mate);
			mc.originTemplate.push_back(t2);
			mc.metaClusterSize.push_back(template2.metaClustersVec[mate].cluster.size());

			mc.cluster.insert(mc.cluster.end(),template2.metaClustersVec[mate].cluster.begin(), template2.metaClustersVec[mate].cluster.end());
			mc.originSample.insert(mc.originSample.end(),template2.metaClustersVec[mate].originSample.begin(), template2.metaClustersVec[mate].originSample.end());

		}
		if(mateCount !=1 || tCoverVecVec[sCoverVecVec[i][0]].size() == 1) ///
		{
			mergedTemplate.metaClustersVec.push_back(mc);
		}

	}

	//Rcpp::Rcout<< "done merging S side " << endl;
	// merge star rooted in the T side
	for(unsigned int i=0; i<tCoverVecVec.size(); i++)
	{
		metaCluster mc;
		int mateCount = tCoverVecVec[i].size();
		if(mateCount != 1)
		{

			mc.metaCluster.push_back(i);
			mc.originTemplate.push_back(t2); 	
			mc.metaClusterSize.push_back(template2.metaClustersVec[i].cluster.size());

			mc.cluster = template2.metaClustersVec[i].cluster;
			mc.originSample = template2.metaClustersVec[i].originSample;

			for(int j=0; j<mateCount; j++)
			{
				int mate = tCoverVecVec[i][j];
				mc.metaCluster.push_back(mate);
				mc.originTemplate.push_back(t1);
				mc.metaClusterSize.push_back(template1.metaClustersVec[mate].cluster.size());

				mc.cluster.insert(mc.cluster.end(),template1.metaClustersVec[mate].cluster.begin(), template1.metaClustersVec[mate].cluster.end());
				mc.originSample.insert(mc.originSample.end(),template1.metaClustersVec[mate].originSample.begin(), template1.metaClustersVec[mate].originSample.end());

			}

			mergedTemplate.metaClustersVec.push_back(mc);
		}

	}
}

double findMinCostPair(vector<vector<struct templatePair> > templatePairs,int* template1, int* template2 )
{
	int templateCount = templatePairs.size();
	BipartiteGraph::Weight minWeight = 999999.0;
	for(int i=0; i<templateCount; i++)
	{
		for(int j=i+1; j<templateCount; j++)
		{
			BipartiteGraph::Weight weight = templatePairs[i][j].coverWeight;
			if(minWeight > weight)
			{
				minWeight = weight;
				*template1 = i;
				*template2 = j;
			}
		}
	}
	return minWeight;
}



void mySort(vector<int>& samples, vector<int>& clusters)
{
	int len = samples.size();
	int j;
	for(int i=1; i<len; i++)
	{
		int key = samples[i];
		int key1 = clusters[i];
		for(j=i-1; j>=0 && samples[j] > key ;  j--)
		{
			samples[j+1] = samples[j];
			clusters[j+1] = clusters[j];
		}
		samples[j+1] = key;
		clusters[j+1] = key1;

	}
}

void printTemplate(classTemplate ct)
{

	int numMetaCluster = ct.metaClustersVec.size();
	Rcpp::Rcout<< "Number of Meta Cluster = " << numMetaCluster << endl;

	for(int i=0; i<numMetaCluster; i++)
	{
		metaCluster mc = ct.metaClustersVec[i];
		mySort(mc.originSample, mc.cluster);
		int numCluster = mc.cluster.size();
		Rcpp::Rcout<< "mc" << i+1 << " = [";
		for(int j=0; j<numCluster; j++)
		{
			Rcpp::Rcout<<  mc.originSample[j]+1 << " " << mc.cluster[j]+1 ;
			if(j == numCluster-1)
				Rcpp::Rcout<< "];";
			else
				Rcpp::Rcout<<"; ";

		}
		Rcpp::Rcout<< endl;
	}
}

