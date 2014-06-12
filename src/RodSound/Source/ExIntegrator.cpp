//
//  ExIntegrator.cpp
//  Visualizer
//
//  Created by eschweickart on 6/10/14.
//
//

#include "ExIntegrator.h"

typedef Eigen::VectorXf VecXf;

ExIntegrator::ExIntegrator(Yarn& y, std::vector<YarnEnergy*>& energies) : Integrator(y, energies) {}

bool ExIntegrator::integrate(Clock& c) {
  
  size_t NumEqs = y.numCPs() * 3;
  VecXf dqdot = VecXf::Zero(NumEqs);
  
  for (YarnEnergy* e : energies) {
    if (!e->eval(VecXf::Zero(NumEqs), c, dqdot)) return false;
  }
  
  for (int i=0; i<y.numCPs(); i++) {
    // WARNING: assumes the mass matrix is the identity
    Vec3f curdqdot = dqdot.block<3,1>(3*i, 0);
    y.next().points[i].accel = -curdqdot;
    y.next().points[i].vel = y.cur().points[i].vel - curdqdot;
    y.next().points[i].pos = y.cur().points[i].pos + c.timestep() * y.next().points[i].vel;
  }
  
  return true;
}