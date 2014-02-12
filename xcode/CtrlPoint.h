//
//  CtrlPoint.h
//  Visualizer
//
//  Created by eschweickart on 2/12/14.
//
//

#ifndef Visualizer_CtrlPoint_h
#define Visualizer_CtrlPoint_h

#include "Eigen/Dense"

typedef Eigen::Vector3f Vec3f;

struct CtrlPoint
{
  /// Position of the control point in space.
  Vec3f pos;
  /// Velocity of the control point.
  Vec3f vel;
  /// The accumulated force at this control point.
  Vec3f force;
};

#endif