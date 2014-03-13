//
//  Util.h
//  Visualizer
//
//  Created by eschweickart on 2/24/14.
//
//

#ifndef Visualizer_Util_h
#define Visualizer_Util_h

#include <boost/thread/thread.hpp>

//#define PARALLEL
#define NUM_THREADS 4
/// A macro for parallelizing a for loop. f is a closure that takes (startIndex, endIndex) as
/// arguments, threads is an array of `boost::thread`s of size NUM_THREADS, num is the number
/// of elements in the for loop, limit is a threshold such that if num<threshold, then the
/// for loop is executed on a single thread to reduce the overhead of spawning threads.
#define MAYBE_PARALLEL(f, threads, num, limit) \
 if (num<limit) { f(0,num); } \
 else { \
   for(int i=0; i<NUM_THREADS; i++) { \
      threads[i]=boost::thread(f,(i*num)/NUM_THREADS,((i+1)*num)/NUM_THREADS);\
   } \
   for(int i=0; i<NUM_THREADS ;i++) { threads[i].join(); } \
 }

// #define ENABLE_CHECK_NAN

#ifdef ENABLE_CHECK_NAN
#define CHECK_NAN(f) if(isnan(f)) assert(false && "NaN Detected.")
#else
#define CHECK_NAN(f)
#endif
/// A vector-like class that allows an offset to indices. That is, if the only elements stored are
/// in indices 5, 6, and 7, this class allows you to specify an offse (-5 in this case) such that
/// only 3 elements are stored in memory, but they are accessed via indices 5, 6, and 7.
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
  const T& operator[]( const size_t index ) const {
    assert(index+offset>=0 && index+offset<vec.size());
    return vec[index + offset];
  }
  void push_back( T elt ) { vec.push_back(elt); }
  const size_t size() const { return vec.size(); }
  const size_t end() const { return vec.size() + offset; }
  void clear() { vec.clear(); }
};

typedef Eigen::Vector2f Vec2f;

/// A space to keep common data. Currently in flux.
struct Workspace
{
  /// "Thread pool" for parallel procedures. In quotes because, as far as I know, threads are
  /// re-initialized each time boost::thread(f) is called.
  boost::thread threads[NUM_THREADS];
  
  /// Storage space for the material curvature at rest. This is defined at each internal control
  /// point i for edges i-1 and i.
  offVec<offVec<Vec2f>> restMatCurvature;
  
};

#endif
