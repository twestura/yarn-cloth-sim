//
//  ConstraintIntegrator.cpp
//  Visualizer
//
//  Created by eschweickart on 6/6/14.
//
//

#include "ConstraintIntegrator.h"


ConstraintIntegrator::ConstraintIntegrator(Yarn& y, std::vector<YarnEnergy*> energies,
                                           std::vector<YarnConstraint*> constraints) :
Integrator(y, energies), constraints(constraints) {
  
}

bool ConstraintIntegrator::integrate(Clock& c) {
  // Evaluate all (explicit) forces
  size_t numEqs = 3*y.numCPs();
  VecXf dqdot = VecXf::Zero(numEqs);
  
  for (YarnEnergy* e : energies) {
    e->eval(VecXf::Zero(numEqs), c, dqdot);
  }
  
  // TODO: mass matrix may not be the Identity
  // Find candidate positions
  VecXf xStar = VecXf(numEqs);
  for (int i=0; i<y.numCPs(); i++) {
    Vec3f velStar = y.cur().points[i].vel - dqdot.block<3,1>(3*i, 0);
    xStar.block<3,1>(3*i, 0) = y.cur().points[i].pos + c.timestep() * velStar;
  }
  
  // TODO: apply mass scaling (??)
  
  // Pre-stabilization
  
  // Solve constraints
  const int gaussJacobiMaxIter = 4;
  int gaussJacobiIter = 0;
  float omega = 1;
  while (gaussJacobiIter < gaussJacobiMaxIter) {
    for (YarnConstraint* c : constraints) {
      c->eval(xStar, omega); // need number of constraints per point?
    }
    gaussJacobiIter++;
  }
  
  // Update positions/velocities
  for (int i=0; i<y.numCPs(); i++) {
    Vec3f delta = xStar.block<3,1>(3*i, 0) - y.cur().points[i].pos;
    if (delta.norm() < 1e-4) { // FIXME
      y.next().points[i].vel = Vec3f::Zero();
      y.next().points[i].pos = y.cur().points[i].pos;
    } else {
      y.next().points[i].vel = delta / c.timestep();
      y.next().points[i].pos = xStar.block<3,1>(3*i, 0);
    }
  }
  
  return true;
}