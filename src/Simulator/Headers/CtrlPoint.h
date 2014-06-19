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

struct CtrlPoint
{
  /// Position of the control point in space.
  Vec3e pos;
  /// Velocity of the control point.
  Vec3e vel;
  /// Acceleration of the control point.
  Vec3e accel;
};

#endif