//
//  Constraint.cpp
//  Visualizer
//
//  Created by eschweickart on 6/6/14.
//
//

#include "Constraint.h"

// C(p1, p2) = |p1 - p2| - |\bar{p1} - \bar{p2}|
bool Length::eval(VecXf& xStar, float omega) {
  VecXf xDelta = VecXf::Zero(xStar.rows());
  for (int i=0; i<y.numSegs(); i++) {
    Vec3f p1 = xStar.block<3,1>(3*i, 0);
    Vec3f p2 = xStar.block<3,1>(3*(i+1), 0);
    float c = (p1 - p2).norm() - y.rest().segments[i].length();
    Vec3f p1gradC = (p1 - p2).normalized();
    Vec3f p2gradC = -p1gradC;
    // TODO: Assumes mass is identity
    float p1lambda = - c / p1gradC.dot(p1gradC) / 2.0f;
    float p2lambda = - c / p2gradC.dot(p2gradC) / 2.0f;
    
    xDelta.block<3,1>(3*i, 0) += omega * p1gradC * p1lambda; // / (i == 0 ? 1.0f : 2.0f);
    xDelta.block<3,1>(3*(i+1), 0) += omega * p2gradC * p2lambda; // / (i == y.numSegs()-1 ? 1.0f : 2.0f);
  }
  xStar += xDelta;
  return true;
}