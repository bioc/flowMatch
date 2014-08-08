#include <Rcpp.h>
#include <string>
#include "BipartiteGraph.h"


// store the results from a mixed edge cover
struct mixedEdgeCover
{
  std::vector<std::vector<BipartiteGraph::Vertex> > sCoverVecVec;
	std::vector<std::vector<BipartiteGraph::Vertex> > tCoverVecVec;
	double coverWeight;
};



Rcpp::List computeMEC(Rcpp::NumericMatrix dist, double lambda)
{

// create necessary vectors for the MEC function
  int nrow = dist.nrow();
	int ncol = dist.ncol();
	// create a 2D vector to store a distance matrix
    std::vector<std::vector<double> >  distVecVec;
	distVecVec.resize(nrow);
	for(int m=0; m<nrow; m++)
	{
		distVecVec[m].resize(ncol);
		for(int n=0; n<ncol; n++)
		{
			distVecVec[m][n] = dist(m,n);
		}
	}


	// create a bipartite graph from the matrix and compute mixed edge cover 
	BipartiteGraph bgraph(distVecVec);
	mixedEdgeCover mec;
	bgraph.MinWghtGenEdgCover(mec.sCoverVecVec,mec.tCoverVecVec,&mec.coverWeight,lambda);
  
  // create results to be returned to R
  Rcpp::List coverLeft(mec.sCoverVecVec.size());
	Rcpp::List coverRight(mec.tCoverVecVec.size());
	for(int i=0; i<mec.sCoverVecVec.size(); i++)
	{
		for(int j=0; j<mec.sCoverVecVec[i].size(); j++)
		{
			mec.sCoverVecVec[i][j] ++; // for 1 based indexing in R
		}
		Rcpp::NumericVector nvCover( mec.sCoverVecVec[i].begin(), mec.sCoverVecVec[i].end());
		coverLeft[i] =  nvCover;
	}

	for(int i=0; i<mec.tCoverVecVec.size(); i++)
	{
		for(int j=0; j<mec.tCoverVecVec[i].size(); j++)
		{
			mec.tCoverVecVec[i][j] ++; // for 1 based indexing in R
		}
		Rcpp::NumericVector nvCover( mec.tCoverVecVec[i].begin(), mec.tCoverVecVec[i].end());
		coverRight[i] =  nvCover;
	}

	return (Rcpp::List::create(Rcpp::Named("match12") = coverLeft,
			Rcpp::Named("match21") = coverRight,
			Rcpp::Named("matching.cost") = mec.coverWeight));
}



RCPP_MODULE(flowMatch_module){
	using namespace Rcpp ;

	function( "computeMEC" , &computeMEC  , "Match clusters across a pair of sample by mixed edge cover algorithm." ) ;
}





