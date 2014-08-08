//
//  ConstraintIntegrator.cpp
//  Visualizer
//
//  Created by eschweickart on 6/6/14.
//
//

#include "ConstraintIntegrator.h"


ConstraintIntegrator::ConstraintIntegrator(Yarn& y, std::vector<YarnEnergy*>& energies,
                                           std::vector<YarnConstraint*>& constraints) :
Integrator(y, energies), constraints(constraints) {
  
}

bool ConstraintIntegrator::integrate(Clock& c) {
  // Evaluate all (explicit) forces
  VecXe forces = VecXe::Zero(y.numDOF());
  
  for (YarnEnergy* e : energies) {
    e->eval(&forces);
  }
  
  VecXe xStar = y.cur().pos + c.timestep() * (y.cur().vel +
                                              c.timestep() * (y.getInvMass().sparse * forces));
  
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
  y.next().vel = (xStar - y.cur().pos) / c.timestep();
  // TODO: filter velocities to prevent jitter
  y.next().pos = xStar;
  
  return true;
}