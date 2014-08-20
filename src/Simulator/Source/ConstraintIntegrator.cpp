//
//  ConstraintIntegrator.cpp
//  Visualizer
//
//  Created by eschweickart on 6/6/14.
//
//

#include "ConstraintIntegrator.h"
#include "Simulator_Prefix.pch"


ConstraintIntegrator::ConstraintIntegrator(Rod& r, std::vector<RodEnergy*>& energies,
                                           std::vector<RodConstraint*>& constraints) :
Integrator(r, energies), constraints(constraints) {
  
}

bool ConstraintIntegrator::integrate(Clock& c) {
  // Evaluate all (explicit) forces
  VecXe forces = VecXe::Zero(r.numDOF());
  
  for (RodEnergy* e : energies) {
    e->eval(&forces);
  }
  
  VecXe xStar = r.cur().pos + c.timestep() * (r.cur().vel +
                                              c.timestep() * (r.getInvMass().sparse * forces));
  
  // TODO: apply mass scaling (??)
  
  // Pre-stabilization
  
  // Solve constraints
  const int gaussJacobiMaxIter = 4;
  int gaussJacobiIter = 0;
  real omega = 1.0;
  while (gaussJacobiIter < gaussJacobiMaxIter) {
    for (RodConstraint* c : constraints) {
      c->eval(xStar, omega); // need number of constraints per point?
    }
    gaussJacobiIter++;
  }
  
  // Update positions/velocities
  r.next().vel = (xStar - r.cur().pos) / c.timestep();
  // TODO: filter velocities to prevent jitter
  r.next().pos = xStar;
  
  return true;
}