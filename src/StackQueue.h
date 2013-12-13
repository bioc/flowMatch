// List implementation of an indexed priority queue.
//
// Item is assumed to be an integral type and is supposed to be instantiated as
// either Graph::Vertex or BipartiteGraph::Vertex. Value is assumed to be a
// floating point type and it is supposed to be instantiated as either
// Graph::Weight of BipartiteGraph::Weight. This class is a template rather
// than a nontemplate only in order to accomodate Vertex and Weight types for
// both Graph and BipartiteGraph objects.

#ifndef LST_INDEXED_MAX_PRIORITY_QUEUE_H
#define LST_INDEXED_MAX_PRIORITY_QUEUE_H

#include <vector>
#include <list>
#include <algorithm>
#include <utility>
#include <cassert>

template<class Item>
class VecStack {
public:
  typedef int Size;

  explicit VecStack(Size maxSize = 0): mMaxSize(maxSize), mSize(0),
      mItemArr(NULL) {
    std::vector<Item> tmpItemVec;
    tmpItemVec.reserve(mMaxSize);
    tmpItemVec.resize(mMaxSize);
    mItemVec.swap(tmpItemVec);
    if (mMaxSize > 0) {
      mItemArr = &mItemVec[0];
    }
  }

  ~VecStack() {}

  Size GetMaxSize(void) const { return mMaxSize; };

  Size GetSize(void) const { return mSize; };

  bool Empty(void) const { return (mSize == 0); };

  bool Full(void) const { return (mSize == mMaxSize); };

  void Push(const Item& item) {
    assert((item > -1) && (item < mMaxSize));
    assert(mSize < mMaxSize);
    mItemArr[mSize] = item;
    ++mSize;
  }

  void Pop(void) {
    assert(mSize > 0);
    --mSize;
  }

  const Item& Top(void) const {
    assert(mSize > 0);
    return mItemArr[mSize - 1];
  }

  void Clear(void) {
    mSize = 0;
  }

private:
  Size mMaxSize;
  Size mSize;
  std::vector<Item> mItemVec;
  Item* mItemArr;
};






/*
 * VecQueueVecQueueVecQueueVecQueueVecQueueVecQueueVecQueueVecQueueVecQueueVecQueueVecQueueVecQueue
 */

template<class Item>
class VecQueue {
public:
  typedef int Size;

  explicit VecQueue(Size maxSize = 0): mMaxSize(maxSize), mSize(0), mFront(0),
      mBack(0), mItemArr(NULL) {
    std::vector<Item> tmpItemVec;
    tmpItemVec.reserve(mMaxSize);
    tmpItemVec.resize(mMaxSize);
    mItemVec.swap(tmpItemVec);
    if (mMaxSize > 0) {
      mItemArr = &mItemVec[0];
    }
  }

  ~VecQueue() {}

  Size GetMaxSize(void) const { return mMaxSize; };

  Size GetSize(void) const { return mSize; };

  bool Empty(void) const { return (mSize == 0); };

  bool Full(void) const { return (mSize == mMaxSize); };

  void Push(const Item& item) {
    assert((item > -1) && (item < mMaxSize));
    assert(mSize < mMaxSize);
    mItemArr[mBack] = item;
    if (mBack < mMaxSize - 1) {
      ++mBack;
    } else {
      mBack = 0;
    }
    ++mSize;
  }

  void Pop(void) {
    assert(mSize > 0);
    if (mFront < mMaxSize - 1) {
      ++mFront;
    } else {
      mFront = 0;
    }
    --mSize;
  }

  const Item& Front(void) const {
    assert(mSize > 0);
    return mItemArr[mFront];
  }

  void Clear(void) {
    mFront = 0;
    mBack = 0;
    mSize = 0;
  }

private:
  Size mMaxSize;
  Size mSize;
  Size mFront;
  Size mBack;
  std::vector<Item> mItemVec;
  Item* mItemArr;
};



template<class Item>
class LstIndexedQueue {
public:
  typedef int Size;

  explicit LstIndexedQueue(Size maxSize = 0): mMaxSize(maxSize), mSize(0),
      mItArr(NULL) {
    std::vector<typename std::list<Item>::iterator> tmpItVec;
    tmpItVec.reserve(mMaxSize);
    tmpItVec.resize(mMaxSize);
    mItVec.swap(tmpItVec);
    if (mMaxSize > 0) {
      mItArr = &mItVec[0];
    }
    if (mItArr != NULL) {
      std::fill(&mItArr[0], &mItArr[mMaxSize], mItemLst.end());
    }
  }

  ~LstIndexedQueue() {}

  Size GetMaxSize(void) const { return mMaxSize; };

  Size GetSize(void) const { return mSize; };

  bool Empty(void) const { return (mSize == 0); };

  bool Full(void) const { return (mSize == mMaxSize); };

  void Push(const Item& item) {
    assert((item > -1) && (item < mMaxSize));
    assert(mItArr[item] == mItemLst.end());
    mItArr[item] = mItemLst.insert(mItemLst.end(), item);
    ++mSize;
  }

  void Pop(void) {
    assert(mSize > 0);
    mItArr[mItemLst.front()] = mItemLst.end();
    mItemLst.pop_front();
    --mSize;
  }

  void Erase(const Item& item) {
    assert((item > -1) && (item < mMaxSize));
    assert(mItArr[item] != mItemLst.end());
    mItemLst.erase(mItArr[item]);
    mItArr[item] = mItemLst.end();
    --mSize;
  }

  const Item& Front(void) const {
    assert(mSize > 0);
    return mItemLst.front();
  }

  void Clear(void) {
    for (typename std::list<Item>::iterator it = mItemLst.begin();
         it != mItemLst.end(); ++it) {
      mItArr[*it] = mItemLst.end();
    }
    mItemLst.clear();
    mSize = 0;
  }

  Item First(void) const {
    if (mSize == 0) {
      return -1;
    }
    return mItemLst.front();
  }

  Item Next(Item item) const {
    assert((item > -1) && (item < mMaxSize));
    assert(mItArr[item] != mItemLst.end());
    typename std::list<Item>::const_iterator it = mItArr[item];
    if (++it == mItemLst.end()) {
      return -1;
    }
    return *it;
  }

  Item Last(void) const {
    if (mSize == 0) {
      return -1;
    }
    return mItemLst.back();
  }

  Item Previous(Item item) const {
    assert((item > -1) && (item < mMaxSize));
    assert(mItArr[item] != mItemLst.end());
    typename std::list<Item>::const_reverse_iterator it(mItArr[item]);
    if (it == mItemLst.rend()) {
      return -1;
    }
    return *it;
  }

private:
  Size mMaxSize;
  Size mSize;
  std::vector<typename std::list<Item>::iterator> mItVec;
  typename std::list<Item>::iterator* mItArr;
  std::list<Item> mItemLst;
};







template<class Item, class Value>
class LstIndexedMinPriorityQueue {
public:
  typedef int Size;

  explicit LstIndexedMinPriorityQueue(Size maxSize = 0): mMaxSize(maxSize),
      mSize(0), mItArr(NULL) {
    std::vector<typename std::list<std::pair<Item, Value> >::iterator>
      tmpItVec;
    tmpItVec.reserve(mMaxSize);
    tmpItVec.resize(mMaxSize);
    mItVec.swap(tmpItVec);
    if (mMaxSize > 0) {
      mItArr = &mItVec[0];
    }
    if (mItArr != NULL) {
      std::fill(&mItArr[0], &mItArr[mMaxSize], mItemValueLst.end());
    }
  }

  ~LstIndexedMinPriorityQueue() {}

  Size GetMaxSize(void) const { return mMaxSize; };

  Size GetSize(void) const { return mSize; };

  bool Empty(void) const { return (mSize == 0); };

  bool Full(void) const { return (mSize == mMaxSize); };

  void Push(const Item& item, const Value& value) {
    assert((item > -1) && (item < mMaxSize));
    assert(mItArr[item] == mItemValueLst.end());
    mItArr[item] = mItemValueLst.insert(mItemValueLst.end(),
                                        std::pair<Item, Value>(item, value));
    ++mSize;
  }

  void Pop(void) {
    assert(mSize > 0);
    typename std::list<std::pair<Item, Value> >::iterator
      topIt = mItemValueLst.begin();
    typename std::list<std::pair<Item, Value> >::iterator it = topIt;
    for (++it; it != mItemValueLst.end(); ++it) {
      if (it->second < topIt->second) {
        topIt = it;
      }
    }
    mItArr[topIt->first] = mItemValueLst.end();
    mItemValueLst.erase(topIt);
    --mSize;
  }

  void Erase(const Item& item) {
    assert((item > -1) && (item < mMaxSize));
    assert(mItArr[item] != mItemValueLst.end());
    mItemValueLst.erase(mItArr[item]);
    mItArr[item] = mItemValueLst.end();
    --mSize;
  }

  void IncreasePriority(const Item& item, const Value& value) {
    assert((item > -1) && (item < mMaxSize));
    assert(mItArr[item] != mItemValueLst.end());
    typename std::list<std::pair<Item, Value> >::iterator it = mItArr[item];
    assert(value < it->second);
    it->second = value;
  }

  const Item& Top(void) const {
    assert(mSize > 0);
    typename std::list<std::pair<Item, Value> >::const_iterator
      topIt = mItemValueLst.begin();
    typename std::list<std::pair<Item, Value> >::const_iterator it = topIt;
    for (++it; it != mItemValueLst.end(); ++it) {
      if (it->second < topIt->second) {
        topIt = it;
      }
    }
    return topIt->first;
  }

  const Item& Top(Value* value) const {
    assert(mSize > 0);
    typename std::list<std::pair<Item, Value> >::const_iterator
      topIt = mItemValueLst.begin();
    typename std::list<std::pair<Item, Value> >::const_iterator it = topIt;
    for (++it; it != mItemValueLst.end(); ++it) {
      if (it->second < topIt->second) {
        topIt = it;
      }
    }
    *value = topIt->second;
    return topIt->first;
  }

  void Clear(void) {
    for (typename std::list<std::pair<Item, Value> >::iterator
         it = mItemValueLst.begin(); it != mItemValueLst.end(); ++it) {
      mItArr[it->first] = mItemValueLst.end();
    }
    mItemValueLst.clear();
    mSize = 0;
  }

private:
  Size mMaxSize;
  Size mSize;
  std::vector<typename std::list<std::pair<Item, Value> >::iterator> mItVec;
  typename std::list<std::pair<Item, Value> >::iterator* mItArr;
  std::list<std::pair<Item, Value> > mItemValueLst;
};





#endif // LST_INDEXED_MAX_PRIORITY_QUEUE_H
