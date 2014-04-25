//
//  Yarn.h
//  Visualizer
//
//  Created by eschweickart on 2/12/14.
//
//

#ifndef __Visualizer__Yarn__
#define __Visualizer__Yarn__

#include "Segment.h"
#include <vector>

/// A Yarn at a specific point in time.
struct YarnStr
{
public:
  /// Control points (n+1 total) that define the yarn's position.
  std::vector<CtrlPoint> points;
  /// Segments (n total) that define the yarn.
  std::vector<Segment> segments;
  
  // TODO: Need Catmull-Rom interpolation for contact resolution.
};

/// A Yarn stepped though points in time.
class Yarn
{
private:
  /// The yarn at rest.
  YarnStr  restYS;
  /// The yarn at the current point in time.
  YarnStr* curYS;
  /// The yarn at the next point in time.
  YarnStr* nextYS;
public:
  /// Constructs a yarn that is initially at rest given a vector of initial positions.
  /// The U vector for the first segment is then propagated down the yarn by parallel transport
  /// through space.
  Yarn(std::vector<Vec3f>& points, Vec3f u0) {
    for (Vec3f p : points) {
      CtrlPoint cp;
      cp.pos = p;
      cp.vel.setZero();
      cp.accel.setZero();
      restYS.points.push_back(cp);
    }
    curYS  = new YarnStr(restYS);
    nextYS = new YarnStr(restYS);
    for (int i=0; i<numCPs()-1; i++) {
      Segment restSeg(restYS.points[i], restYS.points[i+1], u0);
      restYS.segments.push_back(restSeg);
      
      Segment curSeg(curYS->points[i], curYS->points[i+1], u0);
      curYS->segments.push_back(curSeg);
      
      Segment nextSeg(nextYS->points[i], nextYS->points[i+1], u0);
      nextYS->segments.push_back(nextSeg);
    }
    for (int i=1; i<numSegs(); i++) {
      restYS.segments[i].parallelTransport(restYS.segments[i-1]);
      curYS->segments[i].setU(restYS.segments[i].getU());
      nextYS->segments[i].setU(restYS.segments[i].getU());
    }
  }
  
  /// Get a reference to the yarn at the current point in time.
  YarnStr& cur()  { return *curYS;  }
  /// Get a const reference to the yarn at the current point in time.
  const YarnStr& cur() const  { return *curYS;  }
  /// Get a reference to the yarn at the next point in time.
  YarnStr& next() { return *nextYS; }
  /// Get a const reference to the yarn at the next point in time.
  const YarnStr& next() const { return *nextYS; }
  /// Get a const reference to the yarn at rest.
  const YarnStr& rest() const { return restYS; }
  
  /// Get the number of control points on the yarn.
  const size_t inline numCPs() const { return restYS.points.size(); }
  /// Get the number of control points associated with 2 edges.
  const size_t inline numIntCPs() const { return std::max((int)restYS.points.size()-2, 0); }
  /// Get the number of segments in the yarn.
  const size_t inline numSegs() const { return restYS.segments.size(); }
  
  /// Swaps the current and next yarns, e.g. at the end of a timestep.
  void inline swapYarns() {
    YarnStr* temp = curYS;
    curYS = nextYS;
    nextYS = temp;
  }
  
  ~Yarn() {
    delete curYS;
    delete nextYS;
  }
};

class Spline {
  typedef Eigen::Vector4f Vec4f;
  
  // TODO: make this const?
  Vec3f p[4];
  float ti[4] = {0, 0, 0, 0};
    
  
public:
  
  Spline(const Vec3f p0, const Vec3f p1, const Vec3f p2, const Vec3f p3) {
    p[0] = p0; p[1] = p1; p[2] = p2; p[3] = p3;
    ti[1] = powf((p1-p0).norm(), 0.5);
    ti[2] = powf((p2-p1).norm(), 0.5) + ti[1];
    ti[3] = powf((p3-p2).norm(), 0.5) + ti[2];
  }
  
  Spline(const CtrlPoint& cp0, const CtrlPoint& cp1, const CtrlPoint& cp2, const  CtrlPoint& cp3) :
  Spline(cp0.pos, cp1.pos, cp2.pos, cp3.pos) { }
  
  Spline(const CtrlPoint& cp0, const CtrlPoint& cp1, const CtrlPoint&cp2, bool start) :
  Spline(start ? cp0.pos : (2*cp0.pos) - cp1.pos, start ? cp1.pos : cp0.pos,
         start ? cp2.pos : cp1.pos, start ? (2*cp2.pos) - cp1.pos : cp2.pos) { }

  
  Vec3f eval(float t, bool test) const {
    if (test) {
      t = (1-t)*ti[1] + t*ti[2];
      
      Vec3f l01  = ((ti[1] - t)*p[0] + (t - ti[0])*p[1]) / (ti[1] - ti[0]);
      Vec3f l12  = ((ti[2] - t)*p[1] + (t - ti[1])*p[2]) / (ti[2] - ti[1]);
      Vec3f l23  = ((ti[3] - t)*p[2] + (t - ti[2])*p[3]) / (ti[3] - ti[2]);
      Vec3f l012 = ((ti[2] - t)*l01  + (t - ti[0])*l12)  / (ti[2] - ti[0]);
      Vec3f l123 = ((ti[3] - t)*l12  + (t - ti[1])*l23)  / (ti[3] - ti[1]);
      Vec3f c12  = ((ti[2] - t)*l012 + (t - ti[1])*l123) / (ti[2] - ti[1]);
      return c12;
    }
    
    Vec4f u(t*t*t, t*t, t, 1);
    Vec4f uBasis(u.dot(constants::basis[0]),
                 u.dot(constants::basis[1]),
                 u.dot(constants::basis[2]),
                 u.dot(constants::basis[3]));
    Vec4f x(p[0].x(), p[1].x(), p[2].x(), p[3].x());
    Vec4f y(p[0].y(), p[1].y(), p[2].y(), p[3].y());
    Vec4f z(p[0].z(), p[1].z(), p[2].z(), p[3].z());
    return Vec3f(uBasis.dot(x), uBasis.dot(y), uBasis.dot(z));
  }
  
};

#endif /* defined(__Visualizer__Yarn__) */
