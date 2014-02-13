//
//  Simulator.h
//  Visualizer
//
//  Created by eschweickart on 2/10/14.
//
//

#ifndef __Visualizer__Simulator__
#define __Visualizer__Simulator__

#include <iostream>
#include <boost/thread/thread.hpp>
#include "Yarn.h"

//#define PARALLEL
#define NUM_THREADS 4
#define MAYBE_PARALLEL(f, threads, num, limit) if(num<limit){f(0,num);}else{for(int i=0;i<NUM_THREADS;i++){threads[i]=boost::thread(f,(i*num)/NUM_THREADS,((i+1)*num)/NUM_THREADS);}for(int i=0;i<NUM_THREADS;i++){threads[i].join();}}

typedef Eigen::Vector2f Vec2f;

template <class T>
class offVec {
private:
  std::vector<T> vec;
  int offset;
public:
  offVec<T>() { offset = -1; }
  offVec<T>(int i) : offset(i) {}
  void setOffset(int i) { offset = i; }
  T& operator[]( const size_t index ) {
    assert(index+offset>=0 && index+offset<vec.size());
    return vec[index + offset];
  }
  const T& operator[]( const size_t index ) const { return vec[index + offset]; }
  void push_back( T elt ) { vec.push_back(elt); }
  const size_t size() const { return vec.size(); }
  const size_t end() const { return vec.size() + offset; }
  void clear() { vec.clear(); }
};

struct Workspace
{
  /// "Thread pool" for parallel procedures. In quotes because, as far as I know, threads are
  /// re-initialized each time boost::thread(f) is called.
  boost::thread threads[NUM_THREADS];
  /// Curvature binormals for the centerline. Notice this is undefined for the first and final
  /// control points.
  
  offVec<offVec<Vec2f>> restMatCurvature;
  
};

class Simulator
{
  /// The yarn at time 0.
  Yarn restYarn;
  
  /// The yarn at time t.
  Yarn* curYarn;
  /// The yarn at time t+1.
  Yarn* nextYarn;
  
  /// Space for calculations to avoid repetitive, expensive allocation.
  Workspace ws;
  
  /// Compute F(t+1).
  void computeForces();
public:
  /// Entry point for simulator; called after the app has finished initializing
  void run();
};

#endif /* defined(__Visualizer__Simulator__) */
