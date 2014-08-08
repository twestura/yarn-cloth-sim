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
#define CHECK_NAN_VEC(v) assert(!v.hasNaN() && "NaN Detected.")
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

/// Convenience method for checking if two reals are within some epsilon of each other.
static inline bool is_approx(real f1, real f2, real eps = 1.0e-5) {
  return std::abs(f1 - f2) < eps;
}

/// Convenience method for converting a Vec3e to a Vec3c. Useful for draw-time ops.
const Vec3c static inline EtoC(const Vec3e& v) {
  return Vec3c(v.x(), v.y(), v.z());
}

/// Convenience method for converting a Vec3c to a Vec3e.
const Vec3e static inline CtoE(const Vec3c& v) {
  return Vec3e(v.x, v.y, v.z);
}

/// Returns a parallel transported vector given a previous vector and its orthogonal u component.
Vec3e static parallelTransport(const Vec3e vecPrev, const Vec3e vecCur, const Vec3e uPrev) {
  Vec3e cross = vecPrev.cross(vecCur).normalized();
  real cosT = vecCur.dot(vecPrev)/(vecCur.norm() * vecPrev.norm());
  if (!cross.hasNaN() && cosT < 1.0 && cosT >= -1.0) {
    // Form rotation matrix
    real oneMinusCosT = 1.0 - cosT;
    real sinT = sqrt(1.0 - (cosT * cosT));
    real xyc = cross.x() * cross.y() * oneMinusCosT;
    real xzc = cross.x() * cross.z() * oneMinusCosT;
    real yzc = cross.y() * cross.z() * oneMinusCosT;
    real xs = cross.x() * sinT;
    real ys = cross.y() * sinT;
    real zs = cross.z() * sinT;
    Mat3e rotMat;
    rotMat << cosT + cross.x() * cross.x() * oneMinusCosT, xyc - zs, xzc + ys,
    xyc + zs, cosT + cross.y() * cross.y() * oneMinusCosT, yzc - xs,
    xzc - ys, yzc + xs, cosT + cross.z() * cross.z() * oneMinusCosT;
    return rotMat * uPrev;
  }
  return uPrev;
}

/// A class for profiling operations over several iterations.
class Profiler {
  struct Stopwatch {
    boost::timer t;
    double e = 0.0;
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
    it->second.e = 0.0;
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
  
  /// Prints the amount of time recorded by each timer. Note that timers are not stopped if they are
  /// currently running.
  void static printElapsed() {
    TimerMap::iterator it = map.find("Total");
    bool hasTotal = it != map.end();
    double total = hasTotal ? it->second.e + it->second.t.elapsed() : 0.0;
    it = map.begin();
    while (it != map.end()) {
      double e = elapsed(it->first);
      std::cout << it->first << ": " << e;
      if (hasTotal) {
        std::cout << " (" << 100.0 * e / total << "%)";
      }
      std::cout << "\n";
      ++it;
    }
  }
};

#ifdef ENABLE_PROFILER
#define PROFILER_START(x) Profiler::start(x)
#define PROFILER_STOP(x) Profiler::stop(x)
#define PROFILER_PRINT_ELAPSED() Profiler::printElapsed()
#define PROFILER_RESET_ALL() Profiler::resetAll()
#else
#define PROFILER_START(x) // NOP
#define PROFILER_STOP(x) // NOP
#define PROFILER_PRINT_ELAPSED() // NOP
#define PROFILER_RESET_ALL() // NOP
#endif

/// Call this in a .cpp file to initialize the static profiler struct.
#define DECLARE_PROFILER() Profiler::TimerMap Profiler::map{}

template <typename T>
struct NewTypeStruct {
protected:
  T t;
public:
  NewTypeStruct(T t) : t(t) {}
  inline T& operator*() { return t; }
  inline const T& operator*() const { return t; }
  inline T& operator->() { return t; }
  inline const T& operator->() const { return t; }
};

struct Millimeter;

struct Meter : public NewTypeStruct<real> {
public:
  Meter(real r) : NewTypeStruct<real>(r) {}
  Millimeter toMillimeters(); // { return Millimeter(t * 1000.0); }
};

struct Millimeter : public NewTypeStruct<real> {
  Millimeter(real r) : NewTypeStruct<real>(r) {}
  Meter toMeters() { return Meter(t / 1000.0); }
};

#endif
