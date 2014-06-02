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
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>
#include <boost/timer.hpp>

// #define ENABLE_CHECK_NAN

#ifdef ENABLE_CHECK_NAN
#define CHECK_NAN(f) assert(!isnan(f) && "NaN Detected.")
#define CHECK_NAN_VEC(v) assert(v.allFinite() && "NaN or Inf Detected.")
#else
#define CHECK_NAN(f) // NOP
#define CHECK_NAN_VEC(v) // NOP
#endif

/// A space for parallelizing operations. Currently in flux.
/// See http://stackoverflow.com/questions/19500404/how-to-create-a-thread-pool-using-boost-in-c
/// for an example of how this should be used.
struct Threadpool
{
  // TODO: make these static
  boost::asio::io_service ioService;
  boost::thread_group threads;
};

/// Convenience method for checking if two floats are within some epsilon of each other.
static inline bool is_approx(float f1, float f2, float eps = 1e-5) {
  return std::abs(f1 - f2) < eps;
}

/// Convenience method for converting an Eigen Vec3f to a Cinder Vec3f. Useful for draw-time ops.
const ci::Vec3f static inline toCi(const Eigen::Vector3f& v) {
  return ci::Vec3f(v.x(), v.y(), v.z());
}

/// A class for profiling operations over several iterations.
class Profiler {
  struct Stopwatch {
    boost::timer t;
    double e = 0;
    bool running = true;
  };
  
public:
  typedef std::map<std::string, Stopwatch> TimerMap;
  // This shouldn't really be public, but needs to be so it can be declared in a separate .cpp file.
  static TimerMap map;
  
  /// Start a timer with the given name. If the timer does not already exist, it is created. If the
  /// timer is already running, this method fails.
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
  
  /// Stop a timer with a given name. The timer must already exist and be running, or this method
  /// will fail.
  void static stop(std::string name) {
    TimerMap::iterator it = map.find(name);
    assert(it != map.end() && "Timer doesn't exist!");
    assert(it->second.running && "Timer was already stopped!");
    it->second.e += it->second.t.elapsed();
    it->second.running = false;
  }
  
  /// Reset the timer with the given name (i.e. sets the elapsed time to 0). The timer must exist.
  /// Note that the timer is not stopped if it is currently running.
  void static reset(std::string name) {
    TimerMap::iterator it = map.find(name);
    assert(it != map.end() && "Timer doesn't exist!");
    it->second.t.restart();
    it->second.e = 0;
  }
  
  /// Returns the amount of time (in seconds) that the given timer has recorded. The timer must
  /// exist.
  double static elapsed(std::string name) {
    TimerMap::iterator it = map.find(name);
    assert(it != map.end() && "Timer doesn't exist!");
    if (it->second.running) {
      return it->second.e + it->second.t.elapsed();
    }
    return it->second.e;
  }
  
  /// Resets all timers. Note that timers are not stopped if they are currently running.
  void static resetAll() {
    TimerMap::iterator it = map.begin();
    while (it != map.end()) {
      reset(it->first);
      ++it;
    }
  }
  
  /// Prints the amount of time recorded by each timer. Not that timers are not stopped if they are
  /// currently running.
  void static printElapsed() {
    TimerMap::iterator it = map.begin();
    while (it != map.end()) {
      std::cout << it->first << ": " << elapsed(it->first) << "\n";
      ++it;
    }
  }
};

/// Call this in a .cpp file to initialize the static profiler struct.
#define DECLARE_PROFILER() Profiler::TimerMap Profiler::map{}

#endif
