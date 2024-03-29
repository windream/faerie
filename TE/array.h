/*
  $Id: array.h 5151 2010-03-24 23:58:00Z abehm $

  Copyright (C) 2010 by The Regents of the University of California

  Redistribution of this file is permitted under the terms of the
  BSD license.

  An array that dynamically grows.  Reasons to write our own: (1) STL
  Vector doubles its array size every time we need to do a realloc, which
  causes too much storage overhead; (2) BOOST Array has a static size.

  WARNING: Since internally we use "malloc" and "realloc" to allocate
  memory, it only works for basic types such as int, unsigned, float,
  bool, and address pointer. DO NOT USE THIS CLASS FOR CLASSES SUCH AS
  STRING. FOR THESE CLASSES, YOU SHOULD USE STL VECTOR.

  Date: 05/16/2007
  Author: Chen Li <chenli (at) ics.uci.edu>
*/

#ifndef _array_h_
#define _array_h_

#include <cstdlib>
#include <cassert>
#include <iostream>
#include <vector>
#ifndef _WIN32
#include <tr1/unordered_set>
#else
#include <unordered_set>
#endif
#include <cmath>
#include <cstdio>

using namespace std;
using namespace tr1;

typedef enum {
  REALLOC_ADD,
  REALLOC_MULT
} ReallocPolicy;

template <typename T>
class Array {
private:
  T*        data;
  unsigned  capacity;
  unsigned  elementNum;

  int  te;

  float  incrementalSize;
  ReallocPolicy reallocPolicy;
  unordered_set<T>* elementSet;

  template <typename W> friend ostream& operator<<(ostream& out, const Array<W> &array);

public:
  typedef T* iterator;
  typedef T elementType;

  Array() : capacity(10), elementNum(0), incrementalSize(10.f), reallocPolicy(REALLOC_ADD), elementSet(NULL), te(0){
    init();
  };

  Array(unsigned initCapacity, float incrementalSize = 10.0f, ReallocPolicy reallocPolicy = REALLOC_ADD)
    : capacity(initCapacity), elementNum(0), incrementalSize(incrementalSize),
      reallocPolicy(reallocPolicy), elementSet(NULL) {
    init();
  };

  Array(bool useHashSet,
	unsigned initCapacity,
	float incrementalSize,
	ReallocPolicy reallocPolicy = REALLOC_ADD)
    : capacity(initCapacity), elementNum(0), incrementalSize(incrementalSize),
      reallocPolicy(reallocPolicy) {

    init();
    if(useHashSet) elementSet = new unordered_set<T>;
    else elementSet = NULL;
  };

  // helper for the constructor
  void init() {
    if(capacity == 0) capacity = 10;
    data = (T*)malloc(sizeof(T) * capacity);
    assert(data != NULL);
  };

  // destructor
  ~Array() {
    if(data != NULL)
      free(data);

    if(elementSet != NULL)
      delete elementSet;
  };

  iterator begin() {
    return data;
  }

  iterator end() {
    return &data[elementNum];
  }

  unsigned size() const {
    return elementNum;
  };

  static size_t elementSize() {
    return sizeof(T);
  };

  unsigned compressedSize() const {
    return 0;
  };

  // added by shengyue to manipulate size of array directly to boost speed
  void setSize(unsigned size) {

    capacity = size;
    data = (T *)realloc(data, sizeof(T) * capacity);
//    assert(data != NULL);
    elementNum = size;
  };

  void sette(int t)
  {
      te = t;
  };

  int gette()
  {
      return te;
  };

  unsigned getCapacity() const {
    return capacity;
  };

  void extendCapacity() {
    // calculate new capacity
    switch(reallocPolicy) {
    case REALLOC_ADD: capacity = capacity + (unsigned)incrementalSize;
      break;
    case REALLOC_MULT: capacity = ceil(capacity * incrementalSize);
      break;
    default: capacity = capacity + (unsigned)incrementalSize;
      break;
    }

    data = (T *)realloc(data, sizeof(T) * capacity);
    assert(data != NULL);
  }

  void append(const T &element) {
    if (elementNum >= capacity)
      extendCapacity();

    // store this element
    data[elementNum] = element;
    elementNum++;

    // if there are too many elements, we are willing to build a hashset
    // for this array in order to improve the performance of lookups
    //if (elementNum > 10000) {
    //
    //    if (elementNum == 50000 && elementSet == NULL &&  && rand() % 10
    //    <= 2 ) {
    /*if (elementNum == 10000 && elementSet == NULL) {
      elementSet = new unordered_set<T>;
      for (unsigned i = 0; i < elementNum - 1; i ++)
	elementSet->insert(data[i]);
	}*/

    // store the latest one in the hashset
    if (elementSet != NULL){
      //cout<<"insert to hastable " <<endl;
      elementSet->insert(element);
    }//end if
  };

  // delete an element at a given position (starting from 0)
  void erase(const unsigned position) {

    // check corner cases
    if (elementNum == 0 || position >= elementNum)
      return;

    // shift elements to the left
    // be careful about the usage of the unsigned type, and make sure
    // "elementNum -1" is never negative.
    for (register unsigned i = position + 1; i <= elementNum - 1; i ++) {
      data[i-1] = data[i];
    }

    elementNum --;
  }

  // insert an element at a given position (starting from 0)
  void insert(const unsigned position, const T &element) {

    // check corner cases
    if (position > elementNum) // OK to insert at the last position + 1
      return;

    // allocate space
    if (elementNum >= capacity)
      extendCapacity();

    // shift elements to the right
    // changed from unsigned to int, to fix bug when position == 0 by shengyue
    for (register int i = elementNum - 1; i >= (int)position; i --) {
      data[i+1] = data[i];
    }

    // insert the lement
    data[position] = element;

    elementNum ++;
  }

  // added by ALEX to save memory when final size of array is known
  void finalize() {
    data = (T *)realloc(data, sizeof(T) * elementNum);
    assert(data != NULL);
  }

  //added by jiaheng for delete the last element
  void removeLastElement(){
    if(elementNum>0)
      elementNum--;
  };

  T front() const {
#ifdef ARRAY_COMPRESSION
    return at(0);
#else
    return data[0];
#endif
  };

  T back() const {
#ifdef ARRAY_COMPRESSION
    return at(elementNum-1);
#else
    return data[elementNum - 1];
#endif
  };

  bool empty() const {
    return elementNum == 0;
  };

  T& at(unsigned pos) const {
    return data[pos];
  };

  const T& operator[] (unsigned pos) const {
    return data[pos];
  }

  T& operator[] (unsigned pos) {
    return data[pos];
  }

  //
  // Returns the position of the element between [first, last] that is the
  // smallest among those that are >= the given element.  If the element
  // exists in the array, returns that position.  If the key is greater
  // than all the elements in the array, returns the current element
  // number + 1.
  //
  // WARNING: We assume the array is already sorted in an ascending
  // order. it's the caller's responsibility to keep this order.
  //
  // http://www.fredosaurus.com/notes-cpp/algorithms/searching/binarysearch.html

  unsigned binarySearch(T key,  unsigned first) const {
    //assert(capacity > 0);
    unsigned last  = elementNum - 1;

    while (first <= last) {
      unsigned mid = (first + last) / 2;  // compute mid point.
      if (key > data[mid])
    first = mid + 1;  // repeat search in top half.
      else if (key < data[mid])
       {
        if (mid==0) return first; // avoid last == -1
    last = mid - 1; // repeat search in bottom half.
       }
      else
    return mid;     // found it. return the position
    }

    return first;    // failed to find key

  };

  unsigned binarySearch(T key) const {

    return binarySearch(key,0); // default: starting from 0
  }

  /*
   * Returns the position of the first element equal or smaller than the key.
   * Assumes that the array is sorted in ascending order.
   * The search is done form end to start.
   * The algorithm used is a variant of the Jump Search
   * http://www.nist.gov/dads/HTML/jumpsearch.html
   * where the jump step is initially 1 and it is incremented by 1 after each jump
   * Starts the search at element pos - 1, so the maximum value for pos is
   * the length of the array.
   */
  unsigned jumpIncRevSearch(T key, unsigned pos) const {
    unsigned t = 1;
    // t designates the jump step
    // the jump step is incremented by one every time
    while (pos > t && data[pos - t] > key)
      pos -= t++;

    if (pos > 0)
      pos--;
    // t designates the end position
    // end position is the end of current block or the end of the array
    t = pos > t? pos - t:0;
    while (pos > t && data[pos] > key)
      pos--;
    return pos;

  }

  bool has(const T key) const {

    if (elementSet != NULL)
	return elementSet->find(key) != elementSet->end();

    // if no hash_set is available, we have to do a binary search
    unsigned pos = binarySearch(key, 0);
    if (pos < size() && at(pos) == key)
      return true;

    return false;
  }
};

template <typename T>
ostream& operator<<(ostream& out, const Array<T> &array) {
  out << '[';
  for (unsigned i = 0; i < array.elementNum; i++) {
    if (i != 0)
      out << ", ";
    out << array.data[i];
  }
  out << ']';
  return out;
}

template <typename T>
bool operator== (const Array<T> &a, const Array<T> &b) {
  if (&a == &b)
    return true;
  if (a.size() != b.size())
    return false;
  for (unsigned i = 0; i < a.size(); i++)
    if (a[i] != b[i])
      return false;
  return true;
}

//
// the following functions is for Array operations as intersection, union, difference
// add by Wang Bin on Dec 13, 2007
//

// Array union operation input sorted array
void arrayUnion(Array<unsigned> * srcArray, Array<unsigned> * otherArray, Array<unsigned> *& resultArray);

// Array intersection operation is between two sorted array
void arrayIntersection(Array<unsigned> *& array1, Array<unsigned> *& array2, Array<unsigned> *& res);

// Array difference operation is from 'srcArray' minus 'subArray'
void arrayDifference(Array<unsigned> * srcArray, Array<unsigned> *& subArray);

// Array difference oepration outputs stored in 'resultArray' is from 'srcArray' minus 'subArray'
void arrayDifference(Array<unsigned> * srcArray, Array<unsigned> * subArray, Array<unsigned> *& resultArray);
void intersect(Array<unsigned> * srcArray, Array<unsigned> *& otherArray);


#endif
