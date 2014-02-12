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

class Yarn
{
public:
  /// Control points (n+1 total) that define the yarn's position.
  std::vector<CtrlPoint> points;
  /// Segments (n total) that define the yarn.
  std::vector<Segment> segments;
  /// Get the number of control points on the yarn.
  const size_t inline numCPs() const { return points.size(); }
  /// Get the number of control points associated with 2 edges.
  const size_t inline numIntCPs() const { return fmax(points.size()-2, 0); }
  /// Get the number of segments in the yarn.
  const size_t inline numSegs() const { return segments.size(); }
  
  // TODO: enforce vector size invariants
  // TODO: Need Catmull-Rom interpolation for contact resolution.
};

#endif /* defined(__Visualizer__Yarn__) */
