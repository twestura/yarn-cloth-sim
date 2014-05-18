//
//  Util.h
//  Visualizer
//
//  Created by eschweickart on 2/24/14.
//
//

#ifndef Visualizer_Util_h
#define Visualizer_Util_h

#include "Eigen/Dense"
#include <boost/asio/io_service.hpp>
#include <boost/thread/thread.hpp>
#include <boost/timer.hpp>
#include <boost/bind.hpp>

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

/// A space to keep common data. Currently in flux.
struct Workspace
{
  boost::asio::io_service ioService;
  /// "Thread pool" for parallel procedures. In quotes because, as far as I know, threads are
  /// re-initialized each time boost::thread(f) is called.
  boost::thread_group threads;
  
};

static inline bool is_approx(float f1, float f2) {
  return std::abs(f1 - f2) < 1e-5;
}

class Profiler {
  struct Stopwatch {
    boost::timer t;
    double e = 0;
    bool running = true;
  };
  
public:
  typedef std::map<std::string, Stopwatch> TimerMap;
  static TimerMap map;
  
  void static start(std::string name) {
    TimerMap::iterator it = map.find(name);
    if (it == map.end()) {
      map.insert(std::pair<std::string, Stopwatch>(name, Stopwatch{}));
    } else {
      assert(!it->second.running && "Timer was already running!");
      it->second.t.restart();
      it->second.running = true;
    }
  }
  
  void static stop(std::string name) {
    TimerMap::iterator it = map.find(name);
    assert(it != map.end() && "Timer doesn't exist!");
    assert(it->second.running && "Timer was already stopped!");
    it->second.e += it->second.t.elapsed();
    it->second.running = false;
  }
  
  void static reset(std::string name) {
    TimerMap::iterator it = map.find(name);
    assert(it != map.end() && "Timer doesn't exist!");
    it->second.t.restart();
    it->second.e = 0;
  }
  
  double static elapsed(std::string name) {
    TimerMap::iterator it = map.find(name);
    assert(it != map.end() && "Timer doesn't exist!");
    if (it->second.running) {
      return it->second.e + it->second.t.elapsed();
    }
    return it->second.e;
  }
  
  void static resetAll() {
    TimerMap::iterator it = map.begin();
    while (it != map.end()) {
      reset(it->first);
      ++it;
    }
  }
  
  void static printElapsed() {
    TimerMap::iterator it = map.begin();
    while (it != map.end()) {
      std::cout << it->first << ": " << elapsed(it->first) << "\n";
      ++it;
    }
  }
};

#define DECLARE_PROFILER() Profiler::TimerMap Profiler::map{}

const ci::Vec3f static toCi(const Eigen::Vector3f& v) {
  return ci::Vec3f(v.x(), v.y(), v.z());
}

#endif
