#ifndef BIPARTITE_GRAPH_H
#define BIPARTITE_GRAPH_H

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <stack>
#include <queue>
#include <functional>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <map>
#include "StackQueue.h"


class BipartiteGraph {
public:

  typedef int Size;
  typedef int LongSize;
  typedef int Vertex;
  typedef double Weight;
  typedef unsigned int Mark;
  typedef double LargeSize;
  typedef unsigned int Count;
  typedef double LargeCount;

  enum Error {
    eErrNone = 0,
    eErrFileOpen,
    eErrFileRead,
    eErrInvalidFile,
    eErrInvalidSize,
    eErrInvalidVtx,
    eErrInvalidAdj,
    eErrNoPerfMatching
  };

  enum ComputationType {
    eCardinality = 0,
    eEdgeWeight,
    eVertexWeight
  };


  class Statistics {
  public:
    Statistics(bool reverse,
               ComputationType computationType,
               const std::vector<Size>& cardinalityVec,
               const std::vector<Weight>& weightVec,
               const std::vector<Count>& numGraphSearchesVec,
               const std::vector<LargeCount>& numVisitedVerticesVec,
               const std::vector<Count>& numAugOperationsVec,
               const std::vector<Count>& numRevOperationsVec,
               const std::vector<LargeSize>& aggregatedAugPathLengthVec,
               const std::vector<LargeSize>& aggregatedRevPathLengthVec,
               const std::vector<Size>& minAugPathLengthVec,
               const std::vector<Size>& minRevPathLengthVec,
               const std::vector<Size>& maxAugPathLengthVec,
               const std::vector<Size>& maxRevPathLengthVec,
               const std::vector<Count>& numAugPathLengthsVec,
               const std::vector<Count>& numRevPathLengthsVec):
        mReverse(reverse),
        mComputationType(computationType),
        mCardinalityVec(cardinalityVec),
        mWeightVec(weightVec),
        mNumGraphSearchesVec(numGraphSearchesVec),
        mNumVisitedVerticesVec(numVisitedVerticesVec),
        mNumAugOperationsVec(numAugOperationsVec),
        mNumRevOperationsVec(numRevOperationsVec),
        mAggregatedAugPathLengthVec(aggregatedAugPathLengthVec),
        mAggregatedRevPathLengthVec(aggregatedRevPathLengthVec),
        mMinAugPathLengthVec(minAugPathLengthVec),
        mMinRevPathLengthVec(minRevPathLengthVec),
        mMaxAugPathLengthVec(maxAugPathLengthVec),
        mMaxRevPathLengthVec(maxRevPathLengthVec),
        mNumAugPathLengthsVec(numAugPathLengthsVec),
        mNumRevPathLengthsVec(numRevPathLengthsVec) {
      assert(mCardinalityVec.size() == mWeightVec.size());
      assert(mCardinalityVec.size() == mNumGraphSearchesVec.size());
      assert(mCardinalityVec.size() == mNumVisitedVerticesVec.size());
      assert(mCardinalityVec.size() == mNumAugOperationsVec.size());
      assert(mCardinalityVec.size() == mNumRevOperationsVec.size());
      assert(mCardinalityVec.size() == mAggregatedAugPathLengthVec.size());
      assert(mCardinalityVec.size() == mAggregatedRevPathLengthVec.size());
      assert(mCardinalityVec.size() == mMinAugPathLengthVec.size());
      assert(mCardinalityVec.size() == mMinRevPathLengthVec.size());
      assert(mCardinalityVec.size() == mMaxAugPathLengthVec.size());
      assert(mCardinalityVec.size() == mMaxRevPathLengthVec.size());
      assert(mCardinalityVec.size() == mNumAugPathLengthsVec.size());
      assert(mCardinalityVec.size() == mNumRevPathLengthsVec.size());
    }

    void Print(const std::string& caller) const;

  private:
    bool mReverse;
    ComputationType mComputationType;
    const std::vector<Size>& mCardinalityVec;
    const std::vector<Weight>& mWeightVec;
    const std::vector<Count>& mNumGraphSearchesVec;
    const std::vector<LargeCount>& mNumVisitedVerticesVec;
    const std::vector<Count>& mNumAugOperationsVec;
    const std::vector<Count>& mNumRevOperationsVec;
    const std::vector<LargeSize>& mAggregatedAugPathLengthVec;
    const std::vector<LargeSize>& mAggregatedRevPathLengthVec;
    const std::vector<Size>& mMinAugPathLengthVec;
    const std::vector<Size>& mMinRevPathLengthVec;
    const std::vector<Size>& mMaxAugPathLengthVec;
    const std::vector<Size>& mMaxRevPathLengthVec;
    const std::vector<Count>& mNumAugPathLengthsVec;
    const std::vector<Count>& mNumRevPathLengthsVec;
  };

  enum QueueStackType {
    eVecQueueStack = 0,
    eLstQueueStack
  };

  enum IndexedQueueType {
    eLstIndexedQueue = 0
  };

  enum IndexedPriorityQueueType {
    eLstIndexedPriorityQueue = 0,
    eBinHeapIndexedPriorityQueue,
    eFibHeapIndexedPriorityQueue
  };

  explicit BipartiteGraph(Size sNumVertices = 0,Size tNumVertices = 0);

  BipartiteGraph(std::vector<std::vector<double> >&  distVecVec);

  ~BipartiteGraph() {}

  Size GetSNumVertices(void) const { return mSNumVertices; };

  Size GetTNumVertices(void) const { return mTNumVertices; };

  Error Read(const char* fileName);
 

  Error SetAdjacency(Vertex vertex,
                     const std::vector<Vertex>& vertexVec,
                     const std::vector<Weight>& weightVec,
                     bool reverse);

  //Error SetAdjacency(std::vector<std::vector<Vertex> > sVertexVecVec,
  //                   std::vector<std::vector<Vertex> > tVertexVecVec,
  //                   std::vector<std::vector<Weight> > sWeightVecVec,
  //                   std::vector<std::vector<Weight> > tWeightVecVec);

  Error SetVtxWeights(const std::vector<Weight>& sWeightVec,
                      const std::vector<Weight>& tWeightVec);

  void GetVtxWeights(std::vector<Weight>* sWeightVec,
                     std::vector<Weight>* tWeightVec) const;

  void SetEdgWeightsFromVtxWeights(void);

  void GetVtxDegrees(std::vector<Size>* sDegreeVec,
                     std::vector<Size>* tDegreeVec) const;

  Error GetMatchedVtxWeight(const std::vector<Vertex>& sMateVec,
                            const std::vector<Vertex>& tMateVec,
                            Size cardinality,
                            Weight *sWeight,
                            Weight *tWeight,
                            Weight *weight) const;

  Error GetMatchedEdgWeight(const std::vector<Vertex>& sMateVec,
                            const std::vector<Vertex>& tMateVec,
                            Size cardinality,
                            Weight *weight) const;



//---------- header for min weight edge cover -------------
std::vector<std::vector<Weight> > GetSWeightMatrix() ;
bool AddDummyVtx(Weight lambda);
bool RemoveDummyVtx(
	std::vector<std::vector<Vertex> >& sCoverVecVec,
	std::vector<std::vector<Vertex> >& tCoverVecVec);
bool CopyGraphForEdgeCover() ;
Error MinWghtEdgCover(
	std::vector<std::vector<Vertex> >& sCoverVecVec,
	std::vector<std::vector<Vertex> >& tCoverVecVec,
	Weight* coverWeight) ;

Error MinWghtGenEdgCover(
	std::vector<std::vector<Vertex> >& sCoverVecVec,
	std::vector<std::vector<Vertex> >& tCoverVecVec,
	Weight* coverWeight, Weight lambda);
	
Error GetEdgeCoverFromMatching(
	Weight matchingWeight,
	std::vector<Vertex>& sMateVec, 
	std::vector<Vertex>& tMateVec,
	std::vector<std::vector<Vertex> >& sCoverVecVec,
	std::vector<std::vector<Vertex> >& tCoverVecVec,
	Weight* coverWeight) ;
//---------- header for min weight edge cover ends --------


  // implicitly single source.
  //
  // forced direction of graph searches. if reverse is false then the direction
  // of graph searches is from S to T, otherwise it is from T to S.
  Error ComputeMaxEdgWghtPerfMatching(std::vector<Vertex>& sMateVec,
                                      std::vector<Vertex>& tMateVec,
                                      Size* cardinality,
                                      Weight* weight,
                                      bool initialize,
                                      bool reverse,
                                      QueueStackType queueStackType =
                                        eVecQueueStack,
                                      IndexedQueueType indexedQueueType =
                                        eLstIndexedQueue,
                                      IndexedPriorityQueueType
                                        indexedPriorityQueueType =
                                        eBinHeapIndexedPriorityQueue) const;

  Error ComputeMaxEdgWghtPerfMatching1(Vertex* sMateArr,
                                       Vertex* tMateArr,
                                       Size* cardinality,
                                       Weight* weight,
                                       bool initialize,
                                       bool reverse) const;
  void Print(bool compact = true) const;


	// I made them public for global access
	 mutable std::vector<std::vector<Weight> > mSWeightVecVec;
	 mutable std::vector<std::vector<Weight> > mTWeightVecVec;
private:
  enum Status {
    eStatusIdle = 0,
    eStatusProcessable,
    eStatusProcessed,
    eStatusProcessed2,
    eStatusLast
  };

  mutable Size mSNumVertices;
  mutable Size mTNumVertices;
  mutable std::vector<std::vector<Vertex> > mSVertexVecVec;
  mutable std::vector<std::vector<Vertex> > mTVertexVecVec;
  //mutable std::vector<std::vector<Weight> > mSWeightVecVec;
  //mutable std::vector<std::vector<Weight> > mTWeightVecVec;
  std::vector<Weight> mSWeightVec;
  std::vector<Weight> mTWeightVec;

  static const Size cMaxSize;
  static const Size cInfSize;
  static const Vertex cNullVertex;
  static const Weight cZero;
  static const Weight cPositiveOne;
  static const Weight cPositiveInfinity;
  static const Weight cNegativeInfinity;

  void SortAdjacencies(Size sNumVertices,
                       Size tNumVertices,
                       std::vector<std::vector<Vertex> >* sVertexVecVec,
                       std::vector<std::vector<Vertex> >* tVertexVecVec,
                       std::vector<std::vector<Weight> >* sWeightVecVec,
                       std::vector<std::vector<Weight> >* tWeightVecVec) const;

  bool CheckAdjacenciesForDuplicates(Size sNumVertices,
                                     Size tNumVertices,
                                     const std::vector<std::vector<Vertex> >&
                                       sVertexVecVec,
                                     const std::vector<std::vector<Vertex> >&
                                       tVertexVecVec,
                                     const std::vector<std::vector<Weight> >&
                                       sWeightVecVec,
                                     const std::vector<std::vector<Weight> >&
                                       tWeightVecVec) const;

  bool FormAdjacenciesFromTriples(Size sNumVertices,
                                  Size tNumVertices,
                                  const std::vector<Vertex>& sVertexVec,
                                  const std::vector<Vertex>& tVertexVec,
                                  const std::vector<Weight>& weightVec,
                                  std::vector<std::vector<Vertex> >*
                                    sVertexVecVec,
                                  std::vector<std::vector<Vertex> >*
                                    tVertexVecVec,
                                  std::vector<std::vector<Weight> >*
                                    sWeightVecVec,
                                  std::vector<std::vector<Weight> >*
                                    tWeightVecVec) const;


  class InitDualsForPerfEdgWght {
  public:
    InitDualsForPerfEdgWght(const BipartiteGraph& graph): mGraph(graph) {}
    void operator()(Weight* sDualArr,
                    Weight* tDualArr,
                    bool reverse) const;
  private:
    const BipartiteGraph& mGraph;
  };


  void InitGreedyForEdgWght(Vertex* sMateArr,
                            Vertex* tMateArr,
                            Weight* sDualArr,
                            Weight* tDualArr,
                            Size* cardinality,
                            Weight* weight,
                            Count* numVisitedVertices,
                            bool reverse) const;

  Size Augment(Vertex* sMateArr,
               Vertex* tMateArr,
               Vertex* sPtrArr,
               Vertex sLast,
               Vertex tLast) const;

  template<class IndexedQueue>
  Size Augment(Vertex* sMateArr,
               Vertex* tMateArr,
               Vertex* sPtrArr,
               IndexedQueue& sExposedQue,
               Vertex sLast,
               Vertex tLast) const;



  // implicitly breadth-first.
  template<class Queue,
           class IndexedQueue>
  void ComputeShortestAugPathLengthEdgWght(Vertex* sMateArr,
                                           Vertex* tMateArr,
                                           Weight* sDualArr,
                                           Weight* tDualArr,
                                           Size* sLengthArr,
                                           Size* tLengthArr,
                                           Status* sStatusArr,
                                           Status* tStatusArr,
                                           Queue& sProcessableQue,
                                           Queue& sProcessedQue,
                                           Queue& tProcessedQue,
                                           IndexedQueue& sExposedQue,
                                           Size* shortestAugPathLength,
                                           Count* numVisitedVertices,
                                           bool reverse) const;


  // implicitly breadth-first driven depth-first.
  template<class Queue,
           class Stack,
           class IndexedQueue>
  void FindMaximalSetShortestAugPathsEdgWght(Vertex* sMateArr,
                                             Vertex* tMateArr,
                                             Weight* sDualArr,
                                             Weight* tDualArr,
                                             Vertex* sPtrArr,
                                             Size* sLengthArr,
                                             Size* tLengthArr,
                                             Status* sStatusArr,
                                             Status* tStatusArr,
                                             Stack& sProcessableStk,
                                             IndexedQueue& sExposedQue,
                                             Size shortestAugPathLength,
                                             Queue& sLastQue,
                                             Queue& tLastQue,
                                             Count* numVisitedVertices,
                                             bool reverse) const;

  template<class Queue,
           class IndexedPriorityQueue>
  void ProcessEdgWghtPerfSinglePath(Vertex* sMateArr,
                                    Vertex* tMateArr,
                                    Weight* sDualArr,
                                    Weight* tDualArr,
                                    Vertex* sPtrArr,
                                    Vertex* tPtrArr,
                                    Weight* sDistanceArr,
                                    Weight* tDistanceArr,
                                    Status* sStatusArr,
                                    Status* tStatusArr,
                                    Queue& sProcessableQue,
                                    Queue& sProcessedQue,
                                    IndexedPriorityQueue& tProcessableQue,
                                    Queue& tProcessedQue,
                                    Vertex* sLast,
                                    Vertex* tLast,
                                    Weight* dualShift,
                                    Count* numVisitedVertices,
                                    bool reverse) const;



  // implicitly breadth-first.
  template<class Queue,
           class IndexedPriorityQueue>
  void FindAugPathEdgWghtPerfSingleSource(Vertex* sMateArr,
                                          Vertex* tMateArr,
                                          Weight* sDualArr,
                                          Weight* tDualArr,
                                          Vertex* sPtrArr,
                                          Vertex* tPtrArr,
                                          Weight* sDistanceArr,
                                          Weight* tDistanceArr,
                                          Status* sStatusArr,
                                          Status* tStatusarr,
                                          Queue& sProcessableQue,
                                          Queue& sProcessedQue,
                                          IndexedPriorityQueue&
                                            tProcessableQue,
                                          Queue& tProcessedQue,
                                          Vertex sFirst,
                                          Vertex* sLast,
                                          Vertex* tLast,
                                          Weight* dualShift,
                                          Count* numVisitedVertices,
                                          bool reverse) const;

};

template<class IndexedQueue>
inline BipartiteGraph::Size BipartiteGraph::Augment(Vertex* sMateArr,
    Vertex* tMateArr, Vertex* sPtrArr, IndexedQueue& sExposedQue, Vertex sLast,
    Vertex tLast) const {
  Vertex s = sLast;
  Vertex t = tLast;
  Size length = -1;
  while (s != cNullVertex) {
    Vertex tt = sMateArr[s];
    sMateArr[s] = t;
    tMateArr[t] = s;
    if (sPtrArr[s] == cNullVertex) {
      sExposedQue.Erase(s);
    }
    s = sPtrArr[s];
    t = tt;
    length += 2;
  }
  return length;
}



template<class Queue, class IndexedQueue>
inline void BipartiteGraph::ComputeShortestAugPathLengthEdgWght(
    Vertex* sMateArr, Vertex* tMateArr, Weight* sDualArr, Weight* tDualArr,
    Size* sLengthArr, Size* tLengthArr, Status* sStatusArr, Status* tStatusArr,
    Queue& sProcessableQue, Queue& sProcessedQue, Queue& tProcessedQue,
    IndexedQueue& sExposedQue, Size* shortestAugPathLength,
    Count* numVisitedVertices, bool reverse) const {
  *shortestAugPathLength = cInfSize;
  *numVisitedVertices = 0;
  const std::vector<Vertex>* sVertexVecArr = (reverse == false) ?
    ((mSNumVertices == 0) ? NULL : &mSVertexVecVec[0]) :
    ((mTNumVertices == 0) ? NULL : &mTVertexVecVec[0]);
  const std::vector<Weight>* sWeightVecArr = (reverse == false) ?
    ((mSNumVertices == 0) ? NULL : &mSWeightVecVec[0]) :
    ((mTNumVertices == 0) ? NULL : &mTWeightVecVec[0]);
  Vertex sFirst = sExposedQue.First();
  while (sFirst != cNullVertex) {
    sLengthArr[sFirst] = 0;
    sProcessableQue.Push(sFirst);
    sStatusArr[sFirst] = eStatusProcessable;
    ++(*numVisitedVertices);
    sFirst = sExposedQue.Next(sFirst);
  }
  Size sCurrentLevel = -1;
  while (sProcessableQue.Empty() == false) {
    Vertex s = sProcessableQue.Front();
    sProcessableQue.Pop();
    sProcessedQue.Push(s);
    sStatusArr[s] = eStatusProcessed;
    if (sCurrentLevel < sLengthArr[s] / 2) {
      if (*shortestAugPathLength != cInfSize) {
        break;
      }
      ++sCurrentLevel;
    }
    std::vector<Vertex>::size_type sSize = sVertexVecArr[s].size();
    const Vertex* sVertexArr = (sSize == 0) ? NULL : &sVertexVecArr[s][0];
    const Weight* sWeightArr = (sSize == 0) ? NULL : &sWeightVecArr[s][0];
    for (std::vector<Vertex>::size_type i = 0; i < sSize; ++i) {
      Vertex t = sVertexArr[i];
      if (tStatusArr[t] == eStatusProcessed) {
        continue;
      }
      Weight stWeight = sWeightArr[i];
      Weight slack = (sDualArr[s] + tDualArr[t]) - stWeight;
      if (slack > 0) {
        continue;
      }
      tLengthArr[t] = sLengthArr[s] + 1;
      tProcessedQue.Push(t);
      tStatusArr[t] = eStatusProcessed;
      ++(*numVisitedVertices);
      Vertex ss = tMateArr[t];
      if (ss == cNullVertex) {
        *shortestAugPathLength = tLengthArr[t];
      } else {
        sLengthArr[ss] = tLengthArr[t] + 1;
        sProcessableQue.Push(ss);
        sStatusArr[ss] = eStatusProcessable;
        ++(*numVisitedVertices);
      }
    }
  }
}



template<class Queue, class Stack, class IndexedQueue>
inline void BipartiteGraph::FindMaximalSetShortestAugPathsEdgWght(
    Vertex* sMateArr, Vertex* tMateArr, Weight* sDualArr, Weight* tDualArr,
    Vertex* sPtrArr, Size* sLengthArr, Size* tLengthArr, Status* sStatusArr,
    Status* tStatusArr, Stack& sProcessableStk, IndexedQueue& sExposedQue,
    Size shortestAugPathLength, Queue& sLastQue, Queue& tLastQue,
    Count* numVisitedVertices, bool reverse) const {
  *numVisitedVertices = 0;
  const std::vector<Vertex>* sVertexVecArr = (reverse == false) ?
    ((mSNumVertices == 0) ? NULL : &mSVertexVecVec[0]) :
    ((mTNumVertices == 0) ? NULL : &mTVertexVecVec[0]);
  const std::vector<Weight>* sWeightVecArr = (reverse == false) ?
    ((mSNumVertices == 0) ? NULL : &mSWeightVecVec[0]) :
    ((mTNumVertices == 0) ? NULL : &mTWeightVecVec[0]);
  Vertex sFirst = sExposedQue.First();
  while (sFirst != cNullVertex) {
    sProcessableStk.Push(sFirst);
    ++(*numVisitedVertices);
    while (sProcessableStk.Empty() == false) {
      Vertex s = sProcessableStk.Top();
      sProcessableStk.Pop();
      sStatusArr[s] = eStatusProcessed2;
      std::vector<Vertex>::size_type sSize = sVertexVecArr[s].size();
      const Vertex* sVertexArr = (sSize == 0) ? NULL : &sVertexVecArr[s][0];
      const Weight* sWeightArr = (sSize == 0) ? NULL : &sWeightVecArr[s][0];
      for (std::vector<Vertex>::size_type i = 0; i < sSize; ++i) {
        Vertex t = sVertexArr[i];
        Vertex ss = tMateArr[t];
        if ((ss != cNullVertex) && (sStatusArr[ss] == eStatusProcessed2)) {
          continue;
        }
        Weight stWeight = sWeightArr[i];
        Weight slack = (sDualArr[s] + tDualArr[t]) - stWeight;
        if (slack > 0) {
          continue;
        }
        if (tLengthArr[t] != sLengthArr[s] + 1) {
          continue;
        }
        if (tStatusArr[t] == eStatusLast) {
          continue;
        }
        if (tLengthArr[t] == shortestAugPathLength) {
          if (ss == cNullVertex) {
            sProcessableStk.Clear();
            sLastQue.Push(s);
            tLastQue.Push(t);
            tStatusArr[t] = eStatusLast;
            ++(*numVisitedVertices);
            break;
          }
          continue;
        }
        ++(*numVisitedVertices);
        sPtrArr[ss] = s;
        sProcessableStk.Push(ss);
        ++(*numVisitedVertices);
      }
    }
    sFirst = sExposedQue.Next(sFirst);
  }
}

template<class Queue, class IndexedPriorityQueue>
inline void BipartiteGraph::ProcessEdgWghtPerfSinglePath(Vertex* sMateArr,
    Vertex* tMateArr, Weight* sDualArr, Weight* tDualArr, Vertex* sPtrArr,
    Vertex* tPtrArr, Weight* sDistanceArr, Weight* tDistanceArr,
    Status* sStatusArr, Status* tStatusArr, Queue& sProcessableQue,
    Queue& sProcessedQue, IndexedPriorityQueue& tProcessableQue,
    Queue& tProcessedQue, Vertex* sLast, Vertex* tLast, Weight* dualShift,
    Count* numVisitedVertices, bool reverse) const {
  const std::vector<Vertex>* sVertexVecArr = (reverse == false) ?
    ((mSNumVertices == 0) ? NULL : &mSVertexVecVec[0]) :
    ((mTNumVertices == 0) ? NULL : &mTVertexVecVec[0]);
  const std::vector<Weight>* sWeightVecArr = (reverse == false) ?
    ((mSNumVertices == 0) ? NULL : &mSWeightVecVec[0]) :
    ((mTNumVertices == 0) ? NULL : &mTWeightVecVec[0]);
  Weight currentDistance = 0.0;
  while (true) {
    while (sProcessableQue.Empty() == false) {
      Vertex s = sProcessableQue.Front();
      sProcessableQue.Pop();
      sProcessedQue.Push(s);
      sStatusArr[s] = eStatusProcessed;
      std::vector<Vertex>::size_type sSize = sVertexVecArr[s].size();
      const Vertex* sVertexArr = (sSize == 0) ? NULL : &sVertexVecArr[s][0];
      const Weight* sWeightArr = (sSize == 0) ? NULL : &sWeightVecArr[s][0];
      for (std::vector<Vertex>::size_type i = 0; i < sSize; ++i) {
        Vertex t = sVertexArr[i];
        if (tStatusArr[t] == eStatusProcessed) {
          continue;
        }
        Weight stWeight = sWeightArr[i];
        Weight slack = (sDualArr[s] + tDualArr[t]) - stWeight;
        if (slack <= 0) {
          if (tStatusArr[t] == eStatusProcessable) {
            tProcessableQue.Erase(t);
          }
          tDistanceArr[t] = currentDistance;
          tProcessedQue.Push(t);
          tStatusArr[t] = eStatusProcessed;
          ++(*numVisitedVertices);
          Vertex ss = tMateArr[t];
          if (ss == cNullVertex) {
            *sLast = s;
            *tLast = t;
            *dualShift = currentDistance;
            return;
          }
          sDistanceArr[ss] = currentDistance;
          sPtrArr[ss] = s;
          sProcessableQue.Push(ss);
          sStatusArr[ss] = eStatusProcessable;
          ++(*numVisitedVertices);
        } else {
          Weight newDistance = currentDistance + slack;
          if (tDistanceArr[t] > newDistance) {
            tDistanceArr[t] = newDistance;
            tPtrArr[t] = s;
            if (tStatusArr[t] == eStatusIdle) {
              tProcessableQue.Push(t, tDistanceArr[t]);
              tStatusArr[t] = eStatusProcessable;
            } else {
              tProcessableQue.IncreasePriority(t, tDistanceArr[t]);
            }
          }
        }
      }
    }
    {
      Vertex tMin = cNullVertex;
      Weight tMinDistance = cPositiveInfinity;
      if (tProcessableQue.Empty() == false) {
        tMin = tProcessableQue.Top();
        tMinDistance = tDistanceArr[tMin];
      }
      if (tMin == cNullVertex) {
        return;
      } else {
        tProcessableQue.Erase(tMin);
        tProcessedQue.Push(tMin);
        tStatusArr[tMin] = eStatusProcessed;
        ++(*numVisitedVertices);
        Vertex ss = tMateArr[tMin];
        if (ss == cNullVertex) {
          *sLast = tPtrArr[tMin];
          *tLast = tMin;
          *dualShift = tMinDistance;
          return;
        }
        currentDistance = tMinDistance;
        sDistanceArr[ss] = currentDistance;
        sPtrArr[ss] = tPtrArr[tMin];
        sProcessableQue.Push(ss);
        sStatusArr[ss] = eStatusProcessable;
        ++(*numVisitedVertices);
      }
    }
  }
}



template<class Queue, class IndexedPriorityQueue>
inline void BipartiteGraph::FindAugPathEdgWghtPerfSingleSource(
    Vertex* sMateArr, Vertex* tMateArr, Weight* sDualArr, Weight* tDualArr,
    Vertex* sPtrArr, Vertex* tPtrArr, Weight* sDistanceArr,
    Weight* tDistanceArr, Status* sStatusArr, Status* tStatusArr,
    Queue& sProcessableQue, Queue& sProcessedQue,
    IndexedPriorityQueue& tProcessableQue, Queue& tProcessedQue, Vertex sFirst,
    Vertex* sLast, Vertex* tLast, Weight* dualShift, Count* numVisitedVertices,
    bool reverse) const {
  *sLast = cNullVertex;
  *tLast = cNullVertex;
  *dualShift = cPositiveInfinity;
  sDistanceArr[sFirst] = 0.0;
  sProcessableQue.Push(sFirst);
  sStatusArr[sFirst] = eStatusProcessable;
  *numVisitedVertices = 1;
  ProcessEdgWghtPerfSinglePath<Queue, IndexedPriorityQueue>
    (sMateArr, tMateArr, sDualArr, tDualArr, sPtrArr, tPtrArr, sDistanceArr,
     tDistanceArr, sStatusArr, tStatusArr, sProcessableQue, sProcessedQue,
     tProcessableQue, tProcessedQue, sLast, tLast, dualShift,
     numVisitedVertices, reverse);
}

#endif // BIPARTITE_GRAPH_H
