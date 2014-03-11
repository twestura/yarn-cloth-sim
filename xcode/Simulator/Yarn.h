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
      cp.force.setZero();
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

#endif /* defined(__Visualizer__Yarn__) */
