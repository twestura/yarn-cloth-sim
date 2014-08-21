//
//  Constraint.cpp
//  Visualizer
//
//  Created by eschweickart on 6/6/14.
//
//

#include "Constraint.h"

// C(p1, p2) = |p1 - p2| - |\bar{p1} - \bar{p2}|
/// A constraint enforcing that the distance between control points does not change.
bool Length::eval(VecXe& xStar, real omega) {
  VecXe xDelta = VecXe::Zero(xStar.rows());
  for (int i=0; i<r.numEdges(); i++) {
    Vec3e p1 = xStar.segment<3>(3*i);
    Vec3e p2 = xStar.segment<3>(3*(i+1));
    real c = (p1 - p2).norm() - r.rest().edgeLength(i);
    Vec3e p1gradC = (p1 - p2).normalized();
    Vec3e p2gradC = -p1gradC;
    // TODO: Assumes mass is identity
    real p1lambda = - c / p1gradC.dot(p1gradC) / 2.0;
    real p2lambda = - c / p2gradC.dot(p2gradC) / 2.0;
    
    xDelta.segment<3>(3*i) += omega * p1gradC * p1lambda; // / (i == 0 ? 1.0 : 2.0);
    xDelta.segment<3>(3*(i+1)) += omega * p2gradC * p2lambda; // / (i == r.numEdges()-1 ? 1.0 : 2.0);
  }
  xStar += xDelta;
  return true;
}