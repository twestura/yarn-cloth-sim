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
  size_t dof = 3*y.numCPs();
  VecXe forces = VecXe::Zero(dof);
  
  for (YarnEnergy* e : energies) {
    e->eval(&forces);
  }
  
  // Find candidate positions
  VecXe xStar = VecXe(dof);
  for (int i=0; i<y.numCPs(); i++) {
    Vec3e velStar = y.cur().points[i].vel + y.getInvMass().diag(i)*c.timestep()*forces.segment<3>(3*i);
    xStar.segment<3>(3*i) = y.cur().points[i].pos + c.timestep()*velStar;
  }
  
  // TODO: apply mass scaling (??)
  
  // Pre-stabilization
  
  // Solve constraints
  const int gaussJacobiMaxIter = 4;
  int gaussJacobiIter = 0;
  real omega = 1.0;
  while (gaussJacobiIter < gaussJacobiMaxIter) {
    for (YarnConstraint* c : constraints) {
      c->eval(xStar, omega); // need number of constraints per point?
    }
    gaussJacobiIter++;
  }
  
  // Update positions/velocities
  for (int i=0; i<y.numCPs(); i++) {
    Vec3e delta = xStar.segment<3>(3*i) - y.cur().points[i].pos;
    if (delta.norm() < 1.0e-4) { // FIXME
      y.next().points[i].vel = Vec3e::Zero();
      y.next().points[i].pos = y.cur().points[i].pos;
    } else {
      y.next().points[i].vel = delta / c.timestep();
      y.next().points[i].pos = xStar.segment<3>(3*i);
    }
  }
  
  return true;
}