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
#include "Eigen/Dense"
#include <boost/thread/thread.hpp>

#define NUM_THREADS 4

typedef Eigen::Vector3f Vec3f;

struct CtrlPoint
{
  /// Position of the control point in space.
  Vec3f pos;
  /// The accumulated force at this control point.
  Vec3f force;
};


// TODO: move to separate class?
class Segment
{
  /// Reference to the first control point that defines this segment.
  const CtrlPoint& first;
  /// Reference to the second control point that defines this segment.
  const CtrlPoint& second;
  /// The u vector of the twist-free reference (Bishop) frame.
  Vec3f u;
  /// The v vector of the twist-free reference (Bishop) frame.
  Vec3f v;
  /// Rotation in radians of the material frame.
  float rot = 0;

public:
  /// Default Segment constructor
  Segment(const CtrlPoint&, const CtrlPoint&);
  /// Get the vector that represents this segment. Not necessarily unit length.
  const Vec3f inline vec() const;
  /// Length of the segment at rest.
  float restLength;
  
  // TODO: Need Catmull-Rom interpolation for contact resolution.
};


// TODO: move to separate class?
class Yarn
{
public:
  /// Control points (n+1 total) that define the yarn's position.
  std::vector<CtrlPoint> points;
  /// Segments (n total) that define the yarn.
  std::vector<Segment> segments;
  /// Get the number of control points on the yarn.
  const size_t inline numCPs() const;
  /// Get the number of control points associated with 2 edges.
  const size_t inline numIntCPs() const;
  /// Get the number of segments in the yarn.
  const size_t inline numSegs() const;
  
  // TODO: enforce vector size invariants
};

struct Workspace
{
  /// "Thread pool" for parallel procedures. In quotes because, as far as I know, threads are
  /// re-initialized each time boost::thread(f) is called.
  boost::thread threads[NUM_THREADS];
  /// Curvature binormals for the centerline. Notice this is undefined for the first and final
  /// control points.
  std::vector<Vec3f> curvatureBinormals;
};

class Simulator
{
  Yarn y;
  Workspace ws;
  
  void computeForces();
public:
  /// Entry point for simulator; called after the app has finished initializing
  void run();
};

#endif /* defined(__Visualizer__Simulator__) */
