#include "BipartiteGraph.h"
#include <Rcpp.h>
using namespace std;


const BipartiteGraph::Size BipartiteGraph::cMaxSize =
		static_cast<BipartiteGraph::Size>(0x7fffffff);
const BipartiteGraph::Size BipartiteGraph::cInfSize =
		static_cast<BipartiteGraph::Size>(0xffffffff);
const BipartiteGraph::Vertex BipartiteGraph::cNullVertex =
		static_cast<BipartiteGraph::Vertex>(0xffffffff);

const BipartiteGraph::Weight BipartiteGraph::cZero =
		static_cast<BipartiteGraph::Weight>(0.0);
const BipartiteGraph::Weight BipartiteGraph::cPositiveOne =
		static_cast<BipartiteGraph::Weight>(1.0);
const BipartiteGraph::Weight BipartiteGraph::cPositiveInfinity =
		BipartiteGraph::cPositiveOne / BipartiteGraph::cZero;
const BipartiteGraph::Weight BipartiteGraph::cNegativeInfinity =
		-BipartiteGraph::cPositiveInfinity;


BipartiteGraph::BipartiteGraph(Size sNumVertices, Size tNumVertices):
    						mSNumVertices(sNumVertices), mTNumVertices(tNumVertices) {
	vector<vector<Vertex> > tmpSVertexVecVec;
	tmpSVertexVecVec.reserve(sNumVertices);
	tmpSVertexVecVec.resize(sNumVertices);
	mSVertexVecVec.swap(tmpSVertexVecVec);
	vector<vector<Vertex> > tmpTVertexVecVec;
	tmpTVertexVecVec.reserve(tNumVertices);
	tmpTVertexVecVec.resize(tNumVertices);
	mTVertexVecVec.swap(tmpTVertexVecVec);
	vector<vector<Weight> > tmpSWeightVecVec;
	tmpSWeightVecVec.reserve(sNumVertices);
	tmpSWeightVecVec.resize(sNumVertices);
	mSWeightVecVec.swap(tmpSWeightVecVec);
	vector<vector<Weight> > tmpTWeightVecVec;
	tmpTWeightVecVec.reserve(tNumVertices);
	tmpTWeightVecVec.resize(tNumVertices);
	mTWeightVecVec.swap(tmpTWeightVecVec);
	vector<Weight> tmpSWeightVec;
	tmpSWeightVec.reserve(sNumVertices);
	tmpSWeightVec.resize(sNumVertices);
	mSWeightVec.swap(tmpSWeightVec);
	vector<Weight> tmpTWeightVec;
	tmpTWeightVec.reserve(tNumVertices);
	tmpTWeightVec.resize(tNumVertices);
	mTWeightVec.swap(tmpTWeightVec);
}

// a constructor that creates only a complete bipartite graph
BipartiteGraph::BipartiteGraph(vector<vector<double> >&  distVecVec)
{
	mSNumVertices = distVecVec.size();
	mTNumVertices = distVecVec[0].size();

	mSVertexVecVec.resize(mSNumVertices);
	mTVertexVecVec.resize(mTNumVertices);
	mSWeightVecVec.resize(mSNumVertices);
	mTWeightVecVec.resize(mTNumVertices);
	for(int i=0; i<mSNumVertices; i++)
	{
		mSVertexVecVec[i].resize(mTNumVertices);
		mSWeightVecVec[i].resize(mTNumVertices);
	}
	for(int j=0; j<mTNumVertices; j++)
	{
		mTVertexVecVec[j].resize(mSNumVertices);
		mTWeightVecVec[j].resize(mSNumVertices);
	}

	for(int i=0; i<mSNumVertices; i++)
	{
		for(int j=0; j<mTNumVertices; j++)
		{
			mSVertexVecVec[i][j] = j;
			mTVertexVecVec[j][i] = i;
			mSWeightVecVec[i][j] = distVecVec[i][j];
			mTWeightVecVec[j][i] = distVecVec[i][j];
		}
	}
}



BipartiteGraph::Error BipartiteGraph::Read(const char* fileName) {
	// TODO: check max signed and +/-1 shift.
	ifstream fin(fileName);
	if (!fin) {
		return eErrFileOpen;
	}
	// skip comment lines.
	while (fin.peek() == '%') {
		const size_t numChars = 1024;
		char line[numChars];
		fin.getline(line, numChars);
		if (fin.eof()) {
			return eErrInvalidFile;
		}
	}
	Size sNumVertices;
	fin >> sNumVertices;
	if (fin.eof()) {
		return eErrInvalidFile;
	}
	if (sNumVertices < 0) {
		return eErrInvalidSize;
	}
	Size tNumVertices;
	fin >> tNumVertices;
	if (fin.eof()) {
		return eErrInvalidFile;
	}
	if (tNumVertices < 0) {
		return eErrInvalidSize;
	}


	LongSize numEdges;
	fin >> numEdges;
	if (fin.eof()) {
		return eErrInvalidFile;
	}
	if (numEdges < 0)
	{
		return eErrInvalidSize;
	}
	vector<Vertex> sVertexVec;
	sVertexVec.reserve(numEdges);
	sVertexVec.resize(numEdges);
	vector<Vertex> tVertexVec;
	tVertexVec.reserve(numEdges);
	tVertexVec.resize(numEdges);
	vector<Weight> weightVec;
	weightVec.reserve(numEdges);
	weightVec.resize(numEdges);
	vector<Vertex>::size_type numTriples = numEdges;
	Vertex* sVertexArr = (numTriples == 0) ? NULL : &sVertexVec[0];
	Vertex* tVertexArr = (numTriples == 0) ? NULL : &tVertexVec[0];
	Weight* weightArr = (numTriples == 0) ? NULL : &weightVec[0];
	for (vector<Vertex>::size_type i = 0; i < numTriples; ++i) {
		fin >> sVertexArr[i];
		if (fin.eof()) {
			return eErrInvalidFile;
		}
		--(sVertexArr[i]);
		if ((sVertexArr[i] < 0) || (sVertexArr[i] > sNumVertices - 1)) {
			return eErrInvalidVtx;
		}
		fin >> tVertexArr[i];
		if (fin.eof()) {
			return eErrInvalidFile;
		}
		--(tVertexArr[i]);
		if ((tVertexArr[i] < 0) || (tVertexArr[i] > tNumVertices - 1)) {
			return eErrInvalidVtx;
		}
		fin >> weightArr[i];
		if (fin.eof()) {
			return eErrInvalidFile;
		}
		if (weightArr[i] < 0.0) {
			weightArr[i] = -weightArr[i];
		}
	}
	vector<vector<Vertex> > tmpSVertexVecVec;
	vector<vector<Vertex> > tmpTVertexVecVec;
	vector<vector<Weight> > tmpSWeightVecVec;
	vector<vector<Weight> > tmpTWeightVecVec;
	bool success = FormAdjacenciesFromTriples
			(sNumVertices, tNumVertices, sVertexVec, tVertexVec, weightVec,
					&tmpSVertexVecVec, &tmpTVertexVecVec,
					&tmpSWeightVecVec, &tmpTWeightVecVec);
	// adjacencies are storted there , you must keep that since I used it in the indexing in complete bipartite graphs.
	// However in non-complete graphs you can not use such indexing
	if (success == false) {
		return eErrInvalidAdj;
	}
	vector<Weight> tmpSWeightVec;
	tmpSWeightVec.reserve(sNumVertices);
	tmpSWeightVec.resize(sNumVertices);
	vector<Weight> tmpTWeightVec;
	tmpTWeightVec.reserve(tNumVertices);
	tmpTWeightVec.resize(tNumVertices);
	mSNumVertices = sNumVertices;
	mTNumVertices = tNumVertices;
	mSVertexVecVec.swap(tmpSVertexVecVec);
	mTVertexVecVec.swap(tmpTVertexVecVec);
	mSWeightVecVec.swap(tmpSWeightVecVec);
	mTWeightVecVec.swap(tmpTWeightVecVec);
	mSWeightVec.swap(tmpSWeightVec);
	mTWeightVec.swap(tmpTWeightVec);
	return eErrNone;
}



vector<vector<BipartiteGraph::Weight> > BipartiteGraph::GetSWeightMatrix() 
{
	return mSWeightVecVec; 
}

// add two distinguish vertices to the graph
bool BipartiteGraph::AddDummyVtx(Weight lambda)  
{
	Size sNumVertices = GetSNumVertices();
	Size tNumVertices = GetTNumVertices();

	mSVertexVecVec.reserve(sNumVertices+1);	
	mSVertexVecVec.resize(sNumVertices+1);
	mTVertexVecVec.reserve(tNumVertices+1);	
	mTVertexVecVec.resize(tNumVertices+1);
	mSWeightVecVec.reserve(sNumVertices+1);	
	mSWeightVecVec.resize(sNumVertices+1);
	mTWeightVecVec.reserve(tNumVertices+1);	
	mTWeightVecVec.resize(tNumVertices+1);

	mSVertexVecVec[sNumVertices].reserve(tNumVertices+1);	
	mSVertexVecVec[sNumVertices].resize(tNumVertices+1);
	mTVertexVecVec[tNumVertices].reserve(sNumVertices+1);	
	mTVertexVecVec[tNumVertices].resize(sNumVertices+1);
	mSWeightVecVec[sNumVertices].reserve(tNumVertices+1);	
	mSWeightVecVec[sNumVertices].resize(tNumVertices+1);
	mTWeightVecVec[tNumVertices].reserve(sNumVertices+1);	
	mTWeightVecVec[tNumVertices].resize(sNumVertices+1);

	// add edges to S side dummy vertex
	for(int i=0; i< tNumVertices; i++)
	{
		mSVertexVecVec[sNumVertices][i] = i;
		mSWeightVecVec[sNumVertices][i] = lambda;
		mTVertexVecVec[i].push_back(sNumVertices);
		mTWeightVecVec[i].push_back(lambda);

	}
	mSVertexVecVec[sNumVertices][tNumVertices] = tNumVertices;
	mSWeightVecVec[sNumVertices][tNumVertices] = 0;

	// add edges to T side dummy vertex
	for(int i=0; i< sNumVertices; i++)
	{
		mTVertexVecVec[tNumVertices][i] = i;
		mTWeightVecVec[tNumVertices][i] = lambda;
		mSVertexVecVec[i].push_back(tNumVertices);
		mSWeightVecVec[i].push_back(lambda);
	}
	mTVertexVecVec[tNumVertices][sNumVertices] = sNumVertices;
	mTWeightVecVec[tNumVertices][sNumVertices] = 0;

	mTNumVertices ++;
	mSNumVertices ++;

	return true;
}



bool BipartiteGraph::RemoveDummyVtx(
		vector<vector<Vertex> >& sCoverVecVec,
		vector<vector<Vertex> >& tCoverVecVec)
{
	//printEdgCover(sCoverVecVec,tCoverVecVec);

	Size sNumVertices = GetSNumVertices();
	Size tNumVertices = GetTNumVertices();

	assert((unsigned int)sNumVertices == sCoverVecVec.size());
	assert((unsigned int)tNumVertices == tCoverVecVec.size());

	for(unsigned int i=0; i<sCoverVecVec[sNumVertices-1].size();i++)
	{
		Vertex t = sCoverVecVec[sNumVertices-1][i];
		if(t < tNumVertices-1)
		{
			tCoverVecVec[t].erase(tCoverVecVec[t].begin(),tCoverVecVec[t].end()); // be careful .. check
		}
	}

	for(unsigned int i=0; i<tCoverVecVec[tNumVertices-1].size();i++)  
	{
		Vertex s = tCoverVecVec[tNumVertices-1][i];
		if(s < sNumVertices-1)
		{
			sCoverVecVec[s].erase(sCoverVecVec[s].begin(),sCoverVecVec[s].end()); // be careful .. check
		}
	}

	// remove the dummy vtx from cover vector
	sCoverVecVec.erase(sCoverVecVec.end());
	tCoverVecVec.erase(tCoverVecVec.end());
	//printEdgCover(sCoverVecVec,tCoverVecVec);
	return true;
}



// create a graph suitable for edge-cover
bool BipartiteGraph::CopyGraphForEdgeCover()  
{	
	// convert it to minimum weight problem by multiplying weight of every edge by -1
	for(int i=0; i<mSNumVertices; i++)
	{
		Size curDegree = mSWeightVecVec[i].size();
		for(int j=0; j<curDegree ; j++)
		{
			mSWeightVecVec[i][j] *= -1;
		}
	}

	for(int i=0; i<mTNumVertices; i++)
	{
		Size curDegree = mTWeightVecVec[i].size();
		for(int j=0; j<curDegree ; j++)
		{
			mTWeightVecVec[i][j] *= -1;
		}
	}

	// make replicated graph .. S vertice wll go to T side and vice versa
	//Rcpp::Rcout << "making replicated graph" << endl;

	Size mSNumVerticesCpy = mSNumVertices + mTNumVertices;
	Size mTNumVerticesCpy = mSNumVerticesCpy;

	mSVertexVecVec.reserve(mSNumVerticesCpy);
	mSVertexVecVec.resize(mSNumVerticesCpy);
	mTVertexVecVec.reserve(mSNumVerticesCpy);
	mTVertexVecVec.resize(mTNumVerticesCpy);
	mSWeightVecVec.reserve(mSNumVerticesCpy);
	mSWeightVecVec.resize(mSNumVerticesCpy);
	mTWeightVecVec.reserve(mSNumVerticesCpy);
	mTWeightVecVec.resize(mTNumVerticesCpy);

	for(int i=mSNumVertices; i<mSNumVerticesCpy; i++)
	{
		unsigned int tIndex = i - mSNumVertices;
		unsigned int tSize = mTVertexVecVec[tIndex].size();
		mSVertexVecVec[i].reserve(tSize);
		mSVertexVecVec[i].resize(tSize);
		for(unsigned int j=0; j<tSize; j++)
		{
			mSVertexVecVec[i][j] = mTVertexVecVec[tIndex][j] + mTNumVertices;
		}


		mSWeightVecVec[i].reserve(tSize);
		//mSWeightVecVec[i].resize(tSize);
		mSWeightVecVec[i].insert(mSWeightVecVec[i].begin(),mTWeightVecVec[tIndex].begin(),mTWeightVecVec[tIndex].end());

	}

	for(int i= mTNumVertices; i<mTNumVerticesCpy; i++)
	{
		unsigned int sIndex = i - mTNumVertices;
		unsigned int sSize = mSVertexVecVec[sIndex].size();
		mTVertexVecVec[i].reserve(sSize);
		mTVertexVecVec[i].resize(sSize);
		for(unsigned int j=0; j<sSize; j++)
		{
			mTVertexVecVec[i][j] = mSVertexVecVec[sIndex][j] + mSNumVertices;
		}


		mTWeightVecVec[i].reserve(sSize);
		//mTWeightVecVec[i].resize(sSize);
		mTWeightVecVec[i].insert(mTWeightVecVec[i].begin(),mSWeightVecVec[sIndex].begin(),mSWeightVecVec[sIndex].end());
	}

	// add cross-copy extra edges with 2*minWeight // since already negated before ..so get the max insdtead
	for(int i=0; i<mSNumVertices; i++)
	{
		Weight maxWeight = *(std::max_element(mSWeightVecVec[i].begin(), mSWeightVecVec[i].end()));
		mSVertexVecVec[i].push_back(mTNumVertices + i);
		mSWeightVecVec[i].push_back(2*maxWeight);
		mTVertexVecVec[mTNumVertices + i].push_back(i);
		mTWeightVecVec[mTNumVertices + i].push_back(2*maxWeight);
	} 

	for(int i=0; i<mTNumVertices; i++)
	{
		Weight maxWeight = *(std::max_element(mTWeightVecVec[i].begin(), mTWeightVecVec[i].end()));
		mTVertexVecVec[i].push_back(mSNumVertices + i);
		mTWeightVecVec[i].push_back(2*maxWeight);
		mSVertexVecVec[mSNumVertices + i].push_back(i);
		mSWeightVecVec[mSNumVertices + i].push_back(2*maxWeight);
	} 
	mSNumVertices += mTNumVertices; 
	mTNumVertices = mSNumVertices;
	return true; // do some error handling 
}

/*
 * Compute minimum weight edge cover from perfect matching in a copy graph.
 * The copy graph is already generated as described in Schinivar (Combinatorial Optimization)
 */

BipartiteGraph::Error BipartiteGraph::GetEdgeCoverFromMatching(
		Weight matchingWeight,
		vector<Vertex>& sMateVec,
		vector<Vertex>& tMateVec,
		vector<vector<Vertex> >& sCoverVecVec,
		vector<vector<Vertex> >& tCoverVecVec,
		Weight* coverWeight)
{
	Size sNumVertices = GetSNumVertices();
	Size tNumVertices = GetTNumVertices();
	//Rcpp::Rcout << 	 sNumVertices << " **" << tNumVertices << endl;
	if (static_cast<Size>(sMateVec.size()) != sNumVertices + tNumVertices) 
	{
		return eErrInvalidSize;
	}
	if (static_cast<Size>(tMateVec.size()) != sNumVertices + tNumVertices)
	{
		return eErrInvalidSize;
	}

	//MyPrintMatching(sMateVec);

	vector<vector<Vertex> > sTmpCoverVecVec;
	sTmpCoverVecVec.reserve(sNumVertices);
	sTmpCoverVecVec.resize(sNumVertices);
	vector<vector<Vertex> > tTmpCoverVecVec;
	tTmpCoverVecVec.reserve(tNumVertices);
	tTmpCoverVecVec.resize(tNumVertices);

	for (Vertex s = 0; s < sNumVertices; ++s) 
	{
		Vertex t = sMateVec[s];

		if (t == cNullVertex) // should not happen since it is a perfect matching
		{
			continue;
		}
		if(t < tNumVertices) // matched with a normal vertex
		{
			sTmpCoverVecVec[s].push_back(t);
			tTmpCoverVecVec[t].push_back(s);
		}
		else // matched with a shadow vertex
		{
			vector<Weight>::iterator it = (std::min_element(mSWeightVecVec[s].begin(), mSWeightVecVec[s].end()));
			int newMateIndex = distance(mSWeightVecVec[s].begin(),it);
			Vertex newMate = mSVertexVecVec[s][newMateIndex];
			sTmpCoverVecVec[s].push_back(newMate);
			tTmpCoverVecVec[newMate].push_back(s);
		}
	}

	for (Vertex t = 0; t < tNumVertices; ++t)
	{
		Vertex s = tMateVec[t];
		if (s == cNullVertex) // should not happen since it is a perfect matching
		{
			continue;
		}
		if(s >= sNumVertices) // a shadow vertx matched with a normal t vertex
		{

			vector<Weight>::iterator it= (std::min_element (mTWeightVecVec[t].begin(), mTWeightVecVec[t].end()));
			int newMateIndex = distance(mTWeightVecVec[t].begin(),it);
			Vertex newMate = mTVertexVecVec[t][newMateIndex];
			tTmpCoverVecVec[t].push_back(newMate);
			sTmpCoverVecVec[newMate].push_back(t);
		}
	}

	//Since there is no path of length >2 , no need to check other direction
	*coverWeight = (matchingWeight * -1) / 2;
	sCoverVecVec.swap (sTmpCoverVecVec);
	tCoverVecVec.swap (tTmpCoverVecVec);

	//printEdgeCover(sCoverVecVec,tCoverVecVec);
	return eErrNone;

}


BipartiteGraph::Error BipartiteGraph::MinWghtEdgCover(
		vector<vector<Vertex> >& sCoverVecVec,
		vector<vector<Vertex> >& tCoverVecVec,
		Weight* coverWeight)
{
	// create a temporary copy graph for computing edge cover
	BipartiteGraph copyGraph(*this); 
	copyGraph.CopyGraphForEdgeCover();
	//copyGraph.Write("test.out");

	// -----compute a maximum edge-weight perfect matching from the copy graph------

	vector<BipartiteGraph::Vertex> sMateVec;
	vector<BipartiteGraph::Vertex> tMateVec;
	BipartiteGraph::Size cardinality;
	BipartiteGraph::Weight matchingWeight;
	BipartiteGraph::Error err;

	map<string, BipartiteGraph::IndexedPriorityQueueType> priorityQueueMap;
	priorityQueueMap["list"] = BipartiteGraph::eLstIndexedPriorityQueue;
	priorityQueueMap["binheap"] = BipartiteGraph::eBinHeapIndexedPriorityQueue;
	priorityQueueMap["fibheap"] = BipartiteGraph::eFibHeapIndexedPriorityQueue;

	map<string, BipartiteGraph::IndexedPriorityQueueType>::const_iterator it = priorityQueueMap.find("list"); // alwasy use list ....
	if (it == priorityQueueMap.end())
	{
		//Rcpp::Rcout << "invalid priority queue data structure" << endl;
		return eErrNone;
	}
	BipartiteGraph::IndexedPriorityQueueType indexedPriorityQueueType = it->second;

	err = copyGraph.ComputeMaxEdgWghtPerfMatching
			(sMateVec, tMateVec, &cardinality, &matchingWeight, true,false,
					BipartiteGraph::eVecQueueStack, BipartiteGraph::eLstIndexedQueue,
					indexedPriorityQueueType);
	if (err == BipartiteGraph::eErrNoPerfMatching)
	{
		//Rcpp::Rcout << "NO PERFECT MATCHING!" << endl;
		return err;
	}

	// -----compute a maximum edge-weight perfect matching from the copy graph ends ------


	//copyGraph.MyPrintMatching(sMateVec); // for debug only

	// -----compute a min-weight edge cover from perfect matching ------
	err = GetEdgeCoverFromMatching (matchingWeight, sMateVec, tMateVec, sCoverVecVec, tCoverVecVec, coverWeight); 

	return eErrNone;
}


BipartiteGraph::Error BipartiteGraph::MinWghtGenEdgCover(
		vector<vector<Vertex> >& sCoverVecVec,
		vector<vector<Vertex> >& tCoverVecVec,
		Weight* coverWeight, Weight lambda)
{
	// create a temporary copy graph for computing edge cover
	BipartiteGraph auxGraph(*this); 
	auxGraph.AddDummyVtx(lambda);
	//auxGraph.Write("test.out");
	auxGraph.MinWghtEdgCover(sCoverVecVec,tCoverVecVec,coverWeight);
	//Rcpp::Rcout << "coverWeight= " << *coverWeight << endl;
	auxGraph.RemoveDummyVtx(sCoverVecVec,tCoverVecVec);
	//Rcpp::Rcout << "coverWeight= " << *coverWeight << endl;
	return eErrNone;
}


BipartiteGraph::Error BipartiteGraph::ComputeMaxEdgWghtPerfMatching(
		vector<Vertex>& sMateVec, vector<Vertex>& tMateVec, Size* cardinality,
		Weight* weight, bool initialize, bool reverse,
		QueueStackType queueStackType, IndexedQueueType indexedQueueType,
		IndexedPriorityQueueType indexedPriorityQueueType) const {
	if (mSNumVertices != mTNumVertices) {
		return eErrNoPerfMatching;
	}
	vector<Vertex> tmpSMateVec;
	tmpSMateVec.reserve(mSNumVertices);
	tmpSMateVec.resize(mSNumVertices);
	sMateVec.swap(tmpSMateVec);
	vector<Vertex> tmpTMateVec;
	tmpTMateVec.reserve(mTNumVertices);
	tmpTMateVec.resize(mTNumVertices);
	tMateVec.swap(tmpTMateVec);
	Vertex* sMateArr = (reverse == false) ?
			((mSNumVertices == 0) ? NULL : &sMateVec[0]) :
			((mTNumVertices == 0) ? NULL : &tMateVec[0]);
	Vertex* tMateArr = (reverse == false) ?
			((mTNumVertices == 0) ? NULL : &tMateVec[0]) :
			((mSNumVertices == 0) ? NULL : &sMateVec[0]);
	InitDualsForPerfEdgWght initDualsForPerfEdgWght(*this);

	return ComputeMaxEdgWghtPerfMatching1(sMateArr, tMateArr, cardinality, weight, initialize, reverse);

}

BipartiteGraph::Error BipartiteGraph::ComputeMaxEdgWghtPerfMatching1(
    Vertex* sMateArr, Vertex* tMateArr, Size* cardinality, Weight* weight,
    bool initialize, bool reverse) const {
  std::vector<Size> cardinalityVec;
  std::vector<Weight> weightVec;
  std::vector<Count> numGraphSearchesVec;
  std::vector<LargeCount> numVisitedVerticesVec;
  std::vector<Count> numAugOperationsVec;
  std::vector<Count> numRevOperationsVec;
  std::vector<LargeSize> aggregatedAugPathLengthVec;
  std::vector<LargeSize> aggregatedRevPathLengthVec;
  std::vector<Size> minAugPathLengthVec;
  std::vector<Size> minRevPathLengthVec;
  std::vector<Size> maxAugPathLengthVec;
  std::vector<Size> maxRevPathLengthVec;
  std::vector<Count> numAugPathLengthsVec;
  std::vector<Count> numRevPathLengthsVec;
  std::vector<Size>::size_type maxNumPhases = 3;
  cardinalityVec.reserve(maxNumPhases);
  weightVec.reserve(maxNumPhases);
  numGraphSearchesVec.reserve(maxNumPhases);
  numVisitedVerticesVec.reserve(maxNumPhases);
  numAugOperationsVec.reserve(maxNumPhases);
  numRevOperationsVec.reserve(maxNumPhases);
  aggregatedAugPathLengthVec.reserve(maxNumPhases);
  aggregatedRevPathLengthVec.reserve(maxNumPhases);
  minAugPathLengthVec.reserve(maxNumPhases);
  minRevPathLengthVec.reserve(maxNumPhases);
  maxAugPathLengthVec.reserve(maxNumPhases);
  maxRevPathLengthVec.reserve(maxNumPhases);
  numAugPathLengthsVec.reserve(maxNumPhases);
  numRevPathLengthsVec.reserve(maxNumPhases);
  Size sNumVertices = (reverse == false) ? mSNumVertices : mTNumVertices;
  Size tNumVertices = (reverse == false) ? mTNumVertices : mSNumVertices;
  std::vector<Weight> sDualVec;
  sDualVec.reserve(sNumVertices);
  sDualVec.resize(sNumVertices);
  Weight* sDualArr = (sNumVertices == 0) ? NULL : &sDualVec[0];
  std::vector<Weight> tDualVec;
  tDualVec.reserve(tNumVertices);
  tDualVec.resize(tNumVertices);
  Weight* tDualArr = (tNumVertices == 0) ? NULL : &tDualVec[0];

  // initialize duals
  InitDualsForPerfEdgWght initDualsFunction1(*this);
  initDualsFunction1(sDualArr, tDualArr, reverse);
  if (initialize == false) {
    if (sMateArr != NULL) {
      std::fill(&sMateArr[0], &sMateArr[sNumVertices], cNullVertex);
    }
    if (tMateArr != NULL) {
      std::fill(&tMateArr[0], &tMateArr[tNumVertices], cNullVertex);
    }
    *cardinality = 0;
  } else {
    Size iniCardinality;
    Weight iniWeight;
    Count iniNumVisitedVertices;
    // greedy initialization.
    InitGreedyForEdgWght(sMateArr, tMateArr, sDualArr, tDualArr,
                         &iniCardinality, &iniWeight, &iniNumVisitedVertices,
                         reverse);
    cardinalityVec.push_back(iniCardinality);
    weightVec.push_back(iniWeight);
    numGraphSearchesVec.push_back(1);
    numVisitedVerticesVec.push_back
      (static_cast<LargeCount>(iniNumVisitedVertices));
    numAugOperationsVec.push_back(static_cast<Count>(iniCardinality));
    numRevOperationsVec.push_back(0);
    aggregatedAugPathLengthVec.push_back
      (static_cast<LargeSize>(iniCardinality));
    aggregatedRevPathLengthVec.push_back(static_cast<LargeSize>(0));
    minAugPathLengthVec.push_back(1);
    minRevPathLengthVec.push_back(0);
    maxAugPathLengthVec.push_back(1);
    maxRevPathLengthVec.push_back(0);
    numAugPathLengthsVec.push_back(1);
    numRevPathLengthsVec.push_back(0);
    *cardinality = iniCardinality;
  }
  std::vector<Vertex> sPtrVec;
  sPtrVec.reserve(sNumVertices);
  sPtrVec.resize(sNumVertices);
  Vertex* sPtrArr = (sNumVertices == 0) ? NULL : &sPtrVec[0];
  if (sPtrArr != NULL) {
    std::fill(&sPtrArr[0], &sPtrArr[sNumVertices], cNullVertex);
  }
  std::vector<Vertex> tPtrVec;
  tPtrVec.reserve(tNumVertices);
  tPtrVec.resize(tNumVertices);
  Vertex* tPtrArr = (tNumVertices == 0) ? NULL : &tPtrVec[0];
  if (tPtrArr != NULL) {
    std::fill(&tPtrArr[0], &tPtrArr[tNumVertices], cNullVertex);
  }
  std::vector<Size> sLengthVec;
  sLengthVec.reserve(sNumVertices);
  sLengthVec.resize(sNumVertices);
  Size* sLengthArr = (sNumVertices == 0) ? NULL : &sLengthVec[0];
  if (sLengthArr != NULL) {
    std::fill(&sLengthArr[0], &sLengthArr[sNumVertices], cInfSize);
  }
  std::vector<Size> tLengthVec;
  tLengthVec.reserve(tNumVertices);
  tLengthVec.resize(tNumVertices);
  Size* tLengthArr = (tNumVertices == 0) ? NULL : &tLengthVec[0];
  if (tLengthArr != NULL) {
    std::fill(&tLengthArr[0], &tLengthArr[tNumVertices], cInfSize);
  }
  std::vector<Weight> sDistanceVec;
  sDistanceVec.reserve(sNumVertices);
  sDistanceVec.resize(sNumVertices);
  Weight* sDistanceArr = (sNumVertices == 0) ? NULL : &sDistanceVec[0];
  if (sDistanceArr != NULL) {
    std::fill(&sDistanceArr[0], &sDistanceArr[sNumVertices],
              cPositiveInfinity);
  }
  std::vector<Weight> tDistanceVec;
  tDistanceVec.reserve(tNumVertices);
  tDistanceVec.resize(tNumVertices);
  Weight* tDistanceArr = (tNumVertices == 0) ? NULL : &tDistanceVec[0];
  if (tDistanceArr != NULL) {
    std::fill(&tDistanceArr[0], &tDistanceArr[tNumVertices],
              cPositiveInfinity);
  }
  std::vector<Status> sStatusVec;
  sStatusVec.reserve(sNumVertices);
  sStatusVec.resize(sNumVertices);
  Status* sStatusArr = (sNumVertices == 0) ? NULL : &sStatusVec[0];
  if (sStatusArr != NULL) {
    std::fill(&sStatusArr[0], &sStatusArr[sNumVertices], eStatusIdle);
  }
  std::vector<Status> tStatusVec;
  tStatusVec.reserve(tNumVertices);
  tStatusVec.resize(tNumVertices);
  Status* tStatusArr = (tNumVertices == 0) ? NULL : &tStatusVec[0];
  if (tStatusArr != NULL) {
    std::fill(&tStatusArr[0], &tStatusArr[tNumVertices], eStatusIdle);
  }
  VecQueue<Vertex> sProcessableQue(sNumVertices);
  VecQueue<Vertex> sProcessedQue(sNumVertices);
  LstIndexedMinPriorityQueue<Vertex, Weight> tProcessableQue(tNumVertices);
  VecQueue<Vertex> tProcessedQue(tNumVertices);
  VecStack<Vertex> sProcessableStk(sNumVertices);
  LstIndexedQueue<Vertex> sExposedQue(sNumVertices);
  for (Vertex s = 0; s < sNumVertices; ++s) {
    if (sMateArr[s] == cNullVertex) {
        sExposedQue.Push(s);
    }
  }
  VecQueue<Vertex> sLastQue(sNumVertices);
  VecQueue<Vertex> tLastQue(tNumVertices);
  Count numGraphSearches = 0;
  LargeCount numVisitedVertices = static_cast<LargeCount>(0);
  Count numAugOperations = 0;
  LargeSize aggregatedAugPathLength = static_cast<LargeSize>(0);
  Size minAugPathLength = 0;
  Size maxAugPathLength = 0;
  Count numAugPathLengths = 0;
  std::vector<Count> numAugsPerAugPathLengthVec;
  numAugsPerAugPathLengthVec.reserve((sNumVertices + tNumVertices) / 2);
  numAugsPerAugPathLengthVec.resize((sNumVertices + tNumVertices) / 2);
  Count* numAugsPerAugPathLengthArr =
    (numAugsPerAugPathLengthVec.empty() == true) ? NULL :
    &numAugsPerAugPathLengthVec[0];
  if (numAugsPerAugPathLengthArr != NULL) {
    std::fill(&numAugsPerAugPathLengthArr[0],
              &numAugsPerAugPathLengthArr[numAugsPerAugPathLengthVec.size()],
              0);
  }
  do {
    Size shortestAugPathLength;
    Count currentNumVisitedVertices;
    ComputeShortestAugPathLengthEdgWght<VecQueue<Vertex>, LstIndexedQueue<Vertex> >
      (sMateArr, tMateArr, sDualArr, tDualArr, sLengthArr, tLengthArr,
       sStatusArr, tStatusArr, sProcessableQue, sProcessedQue, tProcessedQue,
       sExposedQue, &shortestAugPathLength, &currentNumVisitedVertices,
       reverse);
    ++numGraphSearches;
    numVisitedVertices += static_cast<LargeCount>(currentNumVisitedVertices);
    if (shortestAugPathLength != cInfSize) {
      FindMaximalSetShortestAugPathsEdgWght<VecQueue<Vertex>, VecStack<Vertex>, LstIndexedQueue<Vertex> >
        (sMateArr, tMateArr, sDualArr, tDualArr, sPtrArr, sLengthArr,
         tLengthArr, sStatusArr, tStatusArr, sProcessableStk, sExposedQue,
         shortestAugPathLength, sLastQue, tLastQue, &currentNumVisitedVertices,
         reverse);
      ++numGraphSearches;
      numVisitedVertices += static_cast<LargeCount>(currentNumVisitedVertices);
      while (tLastQue.Empty() == false) {
        Vertex sLast = sLastQue.Front();
        sLastQue.Pop();
        Vertex tLast = tLastQue.Front();
        tLastQue.Pop();
        aggregatedAugPathLength += Augment<LstIndexedQueue<Vertex> >
          (sMateArr, tMateArr, sPtrArr, sExposedQue, sLast, tLast);
        ++numAugOperations;
        ++(numAugsPerAugPathLengthArr[shortestAugPathLength / 2]);
        ++(*cardinality);
      }
      maxAugPathLength = shortestAugPathLength;
      ++numAugPathLengths;
    }
    while (sProcessableQue.Empty() == false) {
      Vertex s = sProcessableQue.Front();
      sProcessableQue.Pop();
      sPtrArr[s] = cNullVertex;
      sLengthArr[s] = cInfSize;
      sStatusArr[s] = eStatusIdle;
    }
    while (sProcessedQue.Empty() == false) {
      Vertex s = sProcessedQue.Front();
      sProcessedQue.Pop();
      sPtrArr[s] = cNullVertex;
      sLengthArr[s] = cInfSize;
      sStatusArr[s] = eStatusIdle;
    }
    while (tProcessedQue.Empty() == false) {
      Vertex t = tProcessedQue.Front();
      tProcessedQue.Pop();
      tLengthArr[t] = cInfSize;
      tStatusArr[t] = eStatusIdle;
    }
    if (shortestAugPathLength == cInfSize) {
      break;
    }
  } while (true);
  sExposedQue.Clear();
  for (std::vector<Count>::size_type l = 0;
       l < numAugsPerAugPathLengthVec.size(); ++l) {
    if ((minAugPathLength == 0) && (numAugsPerAugPathLengthArr[l] > 0)) {
      minAugPathLength = 2 * l + 1;
      break;
    }
  }
  cardinalityVec.push_back(*cardinality);
  weightVec.push_back(0.0);
  numGraphSearchesVec.push_back(numGraphSearches);
  numVisitedVerticesVec.push_back(numVisitedVertices);
  numAugOperationsVec.push_back(numAugOperations);
  numRevOperationsVec.push_back(0);
  aggregatedAugPathLengthVec.push_back(aggregatedAugPathLength);
  aggregatedRevPathLengthVec.push_back(static_cast<LargeSize>(0));
  minAugPathLengthVec.push_back(minAugPathLength);
  minRevPathLengthVec.push_back(0);
  maxAugPathLengthVec.push_back(maxAugPathLength);
  maxRevPathLengthVec.push_back(0);
  numAugPathLengthsVec.push_back(numAugPathLengths);
  numRevPathLengthsVec.push_back(0);
  numGraphSearches = 0;
  numVisitedVertices = static_cast<LargeCount>(0);
  numAugOperations = 0;
  aggregatedAugPathLength = static_cast<LargeSize>(0);
  minAugPathLength = 0;
  maxAugPathLength = 0;
  numAugPathLengths = 0;
  if (numAugsPerAugPathLengthArr != NULL) {
    std::fill(&numAugsPerAugPathLengthArr[0],
              &numAugsPerAugPathLengthArr[numAugsPerAugPathLengthVec.size()],
              0);
  }
  for (Vertex sFirst = 0; sFirst < sNumVertices; ++sFirst) {
    if (sMateArr[sFirst] != cNullVertex) {
      continue;
    }
    Vertex sLast;
    Vertex tLast;
    Weight dualShift;
    Count currentNumVisitedVertices;
    FindAugPathEdgWghtPerfSingleSource<VecQueue<Vertex>, LstIndexedMinPriorityQueue<Vertex, Weight> >
      (sMateArr, tMateArr, sDualArr, tDualArr, sPtrArr, tPtrArr, sDistanceArr,
       tDistanceArr, sStatusArr, tStatusArr, sProcessableQue, sProcessedQue,
       tProcessableQue, tProcessedQue, sFirst, &sLast, &tLast, &dualShift,
       &currentNumVisitedVertices, reverse);
    ++numGraphSearches;
    numVisitedVertices += static_cast<LargeCount>(currentNumVisitedVertices);
    if (tLast != cNullVertex) {
      Size augPathLength = Augment(sMateArr, tMateArr, sPtrArr, sLast, tLast);
      aggregatedAugPathLength += augPathLength;
      ++numAugOperations;
      if (maxAugPathLength < augPathLength) {
        maxAugPathLength = augPathLength;
      }
      ++(numAugsPerAugPathLengthArr[augPathLength / 2]);
      ++(*cardinality);
    } else {
      return eErrNoPerfMatching;
    }
    while (sProcessableQue.Empty() == false) {
      Vertex s = sProcessableQue.Front();
      sProcessableQue.Pop();
      sPtrArr[s] = cNullVertex;
      sDistanceArr[s] = cPositiveInfinity;
      sStatusArr[s] = eStatusIdle;
    }
    while (sProcessedQue.Empty() == false) {
      Vertex s = sProcessedQue.Front();
      sProcessedQue.Pop();
      sDualArr[s] -= dualShift - sDistanceArr[s];
      sPtrArr[s] = cNullVertex;
      sDistanceArr[s] = cPositiveInfinity;
      sStatusArr[s] = eStatusIdle;
    }
    while (tProcessableQue.Empty() == false) {
      Vertex t = tProcessableQue.Top();
      tProcessableQue.Erase(t);
      tPtrArr[t] = cNullVertex;
      tDistanceArr[t] = cPositiveInfinity;
      tStatusArr[t] = eStatusIdle;
    }
    while (tProcessedQue.Empty() == false) {
      Vertex t = tProcessedQue.Front();
      tProcessedQue.Pop();
      tDualArr[t] += dualShift - tDistanceArr[t];
      tPtrArr[t] = cNullVertex;
      tDistanceArr[t] = cPositiveInfinity;
      tStatusArr[t] = eStatusIdle;
    }
  }
  for (std::vector<Count>::size_type l = 0;
       l < numAugsPerAugPathLengthVec.size(); ++l) {
    if ((minAugPathLength == 0) && (numAugsPerAugPathLengthArr[l] > 0)) {
      minAugPathLength = 2 * l + 1;
      break;
    }
  }
  for (std::vector<Count>::size_type l = 0;
       l < numAugsPerAugPathLengthVec.size(); ++l) {
    if (numAugsPerAugPathLengthArr[l] > 0) {
      ++numAugPathLengths;
    }
  }
  *weight = 0.0;
  {
    for (Vertex s = 0; s < sNumVertices; ++s) {
      *weight += sDualArr[s];
    }
    for (Vertex t = 0; t < tNumVertices; ++t) {
      *weight += tDualArr[t];
    }
  }
  cardinalityVec.push_back(*cardinality);
  weightVec.push_back(*weight);
  numGraphSearchesVec.push_back(numGraphSearches);
  numVisitedVerticesVec.push_back(numVisitedVertices);
  numAugOperationsVec.push_back(numAugOperations);
  numRevOperationsVec.push_back(0);
  aggregatedAugPathLengthVec.push_back(aggregatedAugPathLength);
  aggregatedRevPathLengthVec.push_back(static_cast<LargeSize>(0));
  minAugPathLengthVec.push_back(minAugPathLength);
  minRevPathLengthVec.push_back(0);
  maxAugPathLengthVec.push_back(maxAugPathLength);
  maxRevPathLengthVec.push_back(0);
  numAugPathLengthsVec.push_back(numAugPathLengths);
  numRevPathLengthsVec.push_back(0);
  Statistics statistics
    (reverse, eEdgeWeight, cardinalityVec, weightVec, numGraphSearchesVec,
     numVisitedVerticesVec, numAugOperationsVec, numRevOperationsVec,
     aggregatedAugPathLengthVec, aggregatedRevPathLengthVec,
     minAugPathLengthVec, minRevPathLengthVec, maxAugPathLengthVec,
     maxRevPathLengthVec, numAugPathLengthsVec, numRevPathLengthsVec);
  //statistics.Print(std::string(__PRETTY_FUNCTION__));
  return eErrNone;
}





void BipartiteGraph::SortAdjacencies(Size sNumVertices, Size tNumVertices,
		vector<vector<Vertex> >* sVertexVecVec,
		vector<vector<Vertex> >* tVertexVecVec,
		vector<vector<Weight> >* sWeightVecVec,
		vector<vector<Weight> >* tWeightVecVec) const {
	// preconditions:
	// - consistent adjacencies (not necessarily sorted and they may contain
	// duplicate edges).
	// postconditions:
	// - consistent + sorted adjacencies (they may contain duplicate edges).
	// TODO: assertions.
	vector<Vertex>* sVertexVecArr =
			(sNumVertices == 0) ? NULL : &(*sVertexVecVec)[0];
	vector<Vertex>* tVertexVecArr =
			(tNumVertices == 0) ? NULL : &(*tVertexVecVec)[0];
	vector<Weight>* sWeightVecArr =
			(sNumVertices == 0) ? NULL : &(*sWeightVecVec)[0];
	vector<Weight>* tWeightVecArr =
			(tNumVertices == 0) ? NULL : &(*tWeightVecVec)[0];
	vector<vector<Vertex>::size_type> sIdxVec;
	sIdxVec.reserve(sNumVertices);
	sIdxVec.resize(sNumVertices);
	vector<Vertex>::size_type* sIdxArr =
			(sNumVertices == 0) ? NULL : &sIdxVec[0];
	if (sIdxArr != NULL) {
		fill(&sIdxArr[0], &sIdxArr[sNumVertices], 0);
	}
	vector<vector<Vertex>::size_type> tIdxVec;
	tIdxVec.reserve(tNumVertices);
	tIdxVec.resize(tNumVertices);
	vector<Vertex>::size_type* tIdxArr =
			(tNumVertices == 0) ? NULL : &tIdxVec[0];
	if (tIdxArr != NULL) {
		fill(&tIdxArr[0], &tIdxArr[tNumVertices], 0);
	}
	// first bucket sort step: s -> t. this sorts the t-adjacencies.
	for (Vertex s = 0; s < sNumVertices; ++s) {
		vector<Vertex>::size_type sSize = sVertexVecArr[s].size();
		const Vertex* sVertexArr = (sSize == 0) ? NULL : &sVertexVecArr[s][0];
		const Weight* sWeightArr = (sSize == 0) ? NULL : &sWeightVecArr[s][0];
		for (vector<Vertex>::size_type i = 0; i < sSize; ++i) {
			Vertex t = sVertexArr[i];
			Weight stWeight = sWeightArr[i];
			tVertexVecArr[t][tIdxArr[t]] = s;
			tWeightVecArr[t][tIdxArr[t]] = stWeight;
			++(tIdxArr[t]);
		}
	}
	// second bucket sort step: t -> s. this sorts the s-adjacencies.
	for (Vertex t = 0; t < tNumVertices; ++t) {
		vector<Vertex>::size_type tSize = tVertexVecArr[t].size();
		const Vertex* tVertexArr = (tSize == 0) ? NULL : &tVertexVecArr[t][0];
		const Weight* tWeightArr = (tSize == 0) ? NULL : &tWeightVecArr[t][0];
		for (vector<Vertex>::size_type j = 0; j < tSize; ++j) {
			Vertex s = tVertexArr[j];
			Weight tsWeight = tWeightArr[j];
			sVertexVecArr[s][sIdxArr[s]] = t;
			sWeightVecArr[s][sIdxArr[s]] = tsWeight;
			++(sIdxArr[s]);
		}
	}
}

bool BipartiteGraph::CheckAdjacenciesForDuplicates(Size sNumVertices,
		Size tNumVertices, const vector<vector<Vertex> >& sVertexVecVec,
		const vector<vector<Vertex> >& tVertexVecVec,
		const vector<vector<Weight> >& sWeightVecVec,
		const vector<vector<Weight> >& tWeightVecVec) const {
	// preconditions:
	// - consistent + sorted adjacencies (they may contain duplicate edges).
	// returns success value:
	// - <true> if the adjacencies do not contain duplicate edges.
	// - <false> otherwise.
	// TODO: assertions.
	const vector<Vertex>* sVertexVecArr =
			(sNumVertices == 0) ? NULL : &sVertexVecVec[0];
	const vector<Vertex>* tVertexVecArr =
			(tNumVertices == 0) ? NULL : &tVertexVecVec[0];
	// check the s-adjacencies again, this time for duplicate t-vertices.
	for (Vertex s = 0; s < sNumVertices; ++s) {
		if (adjacent_find(sVertexVecArr[s].begin(), sVertexVecArr[s].end())
				!= sVertexVecArr[s].end()) {
			return false;
		}
	}
	// check the t-adjacencies again, this time for duplicate s-vertices.
	for (Vertex t = 0; t < tNumVertices; ++t) {
		if (adjacent_find(tVertexVecArr[t].begin(), tVertexVecArr[t].end())
				!= tVertexVecArr[t].end()) {
			return false;
		}
	}
	return true;
}





bool BipartiteGraph::FormAdjacenciesFromTriples(Size sNumVertices,
		Size tNumVertices, const vector<Vertex>& sVertexVec,
		const vector<Vertex>& tVertexVec, const vector<Weight>& weightVec,
		vector<vector<Vertex> >* sVertexVecVec,
		vector<vector<Vertex> >* tVertexVecVec,
		vector<vector<Weight> >* sWeightVecVec,
		vector<vector<Weight> >* tWeightVecVec) const {
	// preconditions:
	// - <sNumVertices> and <tNumVertices> are nonnegative.
	// - the triple vectors have the same size.
	// - every s-vertex is in the [0, <sNumVertices> - 1] range.
	// - every t-vertex is in the [0, <tNumVertices> - 1] range.
	// - every weight is nonnegative.
	// postconditions:
	// - valid adjacencies (consistent + sorted + no duplicate edges), if the
	// triples do no contain duplicate edges.
	// TODO: check max signed.
	// TODO: assertions.
	assert(sVertexVec.size() == tVertexVec.size());
	assert(tVertexVec.size() == weightVec.size());
	// compute the adjacency sizes.
	vector<vector<Vertex>::size_type> sSizeVec;
	sSizeVec.reserve(sNumVertices);
	sSizeVec.resize(sNumVertices);
	vector<Vertex>::size_type* sSizeArr =
			(sNumVertices == 0) ? NULL : &sSizeVec[0];
	if (sSizeArr != NULL) {
		fill(&sSizeArr[0], &sSizeArr[sNumVertices], 0);
	}
	vector<vector<Vertex>::size_type> tSizeVec;
	tSizeVec.reserve(tNumVertices);
	tSizeVec.resize(tNumVertices);
	vector<Vertex>::size_type* tSizeArr =
			(tNumVertices == 0) ? NULL : &tSizeVec[0];
	if (tSizeArr != NULL) {
		fill(&tSizeArr[0], &tSizeArr[tNumVertices], 0);
	}
	vector<Vertex>::size_type numTriples = sVertexVec.size();
	const Vertex* sVertexArr = (numTriples == 0) ? NULL : &sVertexVec[0];
	const Vertex* tVertexArr = (numTriples == 0) ? NULL : &tVertexVec[0];
	const Weight* weightArr = (numTriples == 0) ? NULL : &weightVec[0];
	for (vector<Vertex>::size_type i = 0; i < numTriples; ++i) {
		Vertex s = sVertexArr[i];
		Vertex t = tVertexArr[i];
		++(sSizeArr[s]);
		++(tSizeArr[t]);
	}
	// setup the adjacencies.
	vector<vector<Vertex> > tmpSVertexVecVec;
	tmpSVertexVecVec.reserve(sNumVertices);
	tmpSVertexVecVec.resize(sNumVertices);
	vector<Vertex>* tmpSVertexVecArr =
			(sNumVertices == 0) ? NULL : &tmpSVertexVecVec[0];
	vector<vector<Vertex> > tmpTVertexVecVec;
	tmpTVertexVecVec.reserve(tNumVertices);
	tmpTVertexVecVec.resize(tNumVertices);
	vector<Vertex>* tmpTVertexVecArr =
			(tNumVertices == 0) ? NULL : &tmpTVertexVecVec[0];
	vector<vector<Weight> > tmpSWeightVecVec;
	tmpSWeightVecVec.reserve(sNumVertices);
	tmpSWeightVecVec.resize(sNumVertices);
	vector<Weight>* tmpSWeightVecArr =
			(sNumVertices == 0) ? NULL : &tmpSWeightVecVec[0];
	vector<vector<Weight> > tmpTWeightVecVec;
	tmpTWeightVecVec.reserve(tNumVertices);
	tmpTWeightVecVec.resize(tNumVertices);
	vector<Weight>* tmpTWeightVecArr =
			(tNumVertices == 0) ? NULL : &tmpTWeightVecVec[0];
	for (Vertex s = 0; s < sNumVertices; ++s) {
		tmpSVertexVecArr[s].reserve(sSizeArr[s]);
		tmpSVertexVecArr[s].resize(sSizeArr[s]);
		tmpSWeightVecArr[s].reserve(sSizeArr[s]);
		tmpSWeightVecArr[s].resize(sSizeArr[s]);
	}
	for (Vertex t = 0; t < tNumVertices; ++t) {
		tmpTVertexVecArr[t].reserve(tSizeArr[t]);
		tmpTVertexVecArr[t].resize(tSizeArr[t]);
		tmpTWeightVecArr[t].reserve(tSizeArr[t]);
		tmpTWeightVecArr[t].resize(tSizeArr[t]);
	}
	// transfer data from triples to adjacencies.
	vector<vector<Vertex>::size_type> sIdxVec;
	sIdxVec.reserve(sNumVertices);
	sIdxVec.resize(sNumVertices);
	vector<Vertex>::size_type* sIdxArr =
			(sNumVertices == 0) ? NULL : &sIdxVec[0];
	if (sIdxArr != NULL) {
		fill(&sIdxArr[0], &sIdxArr[sNumVertices], 0);
	}
	vector<vector<Vertex>::size_type> tIdxVec;
	tIdxVec.reserve(tNumVertices);
	tIdxVec.resize(tNumVertices);
	vector<Vertex>::size_type* tIdxArr =
			(tNumVertices == 0) ? NULL : &tIdxVec[0];
	if (tIdxArr != NULL) {
		fill(&tIdxArr[0], &tIdxArr[tNumVertices], 0);
	}
	for (vector<Vertex>::size_type i = 0; i < numTriples; ++i) {
		Vertex s = sVertexArr[i];
		Vertex t = tVertexArr[i];
		Weight stWeight = weightArr[i];
		tmpSVertexVecArr[s][sIdxArr[s]] = t;
		tmpSWeightVecArr[s][sIdxArr[s]] = stWeight;
		++(sIdxArr[s]);
		tmpTVertexVecArr[t][tIdxArr[t]] = s;
		tmpTWeightVecArr[t][tIdxArr[t]] = stWeight;
		++(tIdxArr[t]);
	}
	// sort the adjacencies and check them for duplicate edges.
	SortAdjacencies
	(sNumVertices, tNumVertices, &tmpSVertexVecVec, &tmpTVertexVecVec,
			&tmpSWeightVecVec, &tmpTWeightVecVec);
	bool success =
			CheckAdjacenciesForDuplicates(sNumVertices, tNumVertices,
					tmpSVertexVecVec, tmpTVertexVecVec,
					tmpSWeightVecVec, tmpTWeightVecVec);
	if (success == true) {
		sVertexVecVec->swap(tmpSVertexVecVec);
		tVertexVecVec->swap(tmpTVertexVecVec);
		sWeightVecVec->swap(tmpSWeightVecVec);
		tWeightVecVec->swap(tmpTWeightVecVec);
	}
	return success;
}


void BipartiteGraph::InitDualsForPerfEdgWght::operator()(Weight* sDualArr,
		Weight* tDualArr, bool reverse) const {
	Size sNumVertices = (reverse == false) ?
			mGraph.mSNumVertices : mGraph.mTNumVertices;
	Size tNumVertices = (reverse == false) ?
			mGraph.mTNumVertices : mGraph.mSNumVertices;
	const vector<Vertex>* tVertexVecArr = (reverse == false) ?
			((mGraph.mTNumVertices == 0) ? NULL : &mGraph.mTVertexVecVec[0]) :
			((mGraph.mSNumVertices == 0) ? NULL : &mGraph.mSVertexVecVec[0]);
	const vector<Weight>* sWeightVecArr = (reverse == false) ?
			((mGraph.mSNumVertices == 0) ? NULL : &mGraph.mSWeightVecVec[0]) :
			((mGraph.mTNumVertices == 0) ? NULL : &mGraph.mTWeightVecVec[0]);
	const vector<Weight>* tWeightVecArr = (reverse == false) ?
			((mGraph.mTNumVertices == 0) ? NULL : &mGraph.mTWeightVecVec[0]) :
			((mGraph.mSNumVertices == 0) ? NULL : &mGraph.mSWeightVecVec[0]);
	for (Vertex s = 0; s < sNumVertices; ++s) {
		sDualArr[s] = 0.0;
		vector<Weight>::size_type sSize = sWeightVecArr[s].size();
		const Weight* sWeightArr = (sSize == 0) ? NULL : &sWeightVecArr[s][0];
		for (vector<Vertex>::size_type i = 0; i < sSize; ++i) {
			Weight stWeight = sWeightArr[i];
			if (sDualArr[s] < stWeight) {
				sDualArr[s] = stWeight;
			}
		}
	}
	for (Vertex t = 0; t < tNumVertices; ++t) {
		tDualArr[t] = cNegativeInfinity;
		vector<Vertex>::size_type tSize = tVertexVecArr[t].size();
		const Vertex* tVertexArr = (tSize == 0) ? NULL : &tVertexVecArr[t][0];
		const Weight* tWeightArr = (tSize == 0) ? NULL : &tWeightVecArr[t][0];
		for (vector<Vertex>::size_type j = 0; j < tSize; ++j) {
			Vertex s = tVertexArr[j];
			Weight tsWeight = tWeightArr[j];
			if (tDualArr[t] < tsWeight - sDualArr[s]) {
				tDualArr[t] = tsWeight - sDualArr[s];
			}
		}
	}
}

void BipartiteGraph::InitGreedyForEdgWght(Vertex* sMateArr, Vertex* tMateArr,
		Weight* sDualArr, Weight* tDualArr, Size* cardinality, Weight* weight,
		Count* numVisitedVertices, bool reverse) const {
	Size sNumVertices = (reverse == false) ? mSNumVertices : mTNumVertices;
	Size tNumVertices = (reverse == false) ? mTNumVertices : mSNumVertices;
	if (sMateArr != NULL) {
		fill(&sMateArr[0], &sMateArr[sNumVertices], cNullVertex);
	}
	if (tMateArr != NULL) {
		fill(&tMateArr[0], &tMateArr[tNumVertices], cNullVertex);
	}
	*cardinality = 0;
	*weight = 0.0;
	*numVisitedVertices = 0;
	const vector<Vertex>* sVertexVecArr = (reverse == false) ?
			((mSNumVertices == 0) ? NULL : &mSVertexVecVec[0]) :
			((mTNumVertices == 0) ? NULL : &mTVertexVecVec[0]);
	const vector<Weight>* sWeightVecArr = (reverse == false) ?
			((mSNumVertices == 0) ? NULL : &mSWeightVecVec[0]) :
			((mTNumVertices == 0) ? NULL : &mTWeightVecVec[0]);
	for (Vertex s = 0; s < sNumVertices; ++s) {
		if (sMateArr[s] == cNullVertex) {
			++(*numVisitedVertices);
			vector<Vertex>::size_type sSize = sVertexVecArr[s].size();
			const Vertex* sVertexArr = (sSize == 0) ? NULL : &sVertexVecArr[s][0];
			const Weight* sWeightArr = (sSize == 0) ? NULL : &sWeightVecArr[s][0];
			for (vector<Vertex>::size_type i = 0; i < sSize; ++i) {
				Vertex t = sVertexArr[i];
				Weight stWeight = sWeightArr[i];
				if ((tMateArr[t] == cNullVertex) &&
						(tDualArr[t] == stWeight - sDualArr[s])) {
					sMateArr[s] = t;
					tMateArr[t] = s;
					++(*cardinality);
					*weight += stWeight;
					++(*numVisitedVertices);
					break;
				}
			}
		}
	}

}

inline BipartiteGraph::Size BipartiteGraph::Augment(Vertex* sMateArr,
		Vertex* tMateArr, Vertex* sPtrArr, Vertex sLast, Vertex tLast) const {
	Vertex s = sLast;
	Vertex t = tLast;
	Size length = -1;
	while (s != cNullVertex) {
		Vertex tt = sMateArr[s];
		sMateArr[s] = t;
		tMateArr[t] = s;
		s = sPtrArr[s];
		t = tt;
		length += 2;
	}
	return length;
}






