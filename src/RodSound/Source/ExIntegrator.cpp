//
//  ExIntegrator.cpp
//  Visualizer
//
//  Created by eschweickart on 6/10/14.
//
//

#include "ExIntegrator.h"

typedef Eigen::VectorXf VecXf;

ExIntegrator::ExIntegrator(Yarn& y, std::vector<YarnEnergy*>& energies) : Integrator(y, energies) {
  const float alpha1 = 1.5e-11;
  const float alpha2 = 2e-5;
  
  std::vector<Triplet> triplets;
  size_t NumEqs = y.numCPs() * 3;
  VecXf unused(NumEqs); // Hacky, but we need a non-const vector for eval()
  Clock c(1.0f); // Also hacky, but we need a clock for eval()
  for (YarnEnergy* e : energies) {
    // Assumes external forces have 0 space derivative, which is blatantly false in most cases.
    // Also assumes that everything was intended to be evaluated explicitly
    e->setEvalType(Implicit);
    e->eval(VecXf::Zero(NumEqs), c, unused, &triplets);
    e->setEvalType(Explicit);
  }
  
  Eigen::SparseMatrix<float> stiffness(NumEqs, NumEqs);
  stiffness.setFromTriplets(triplets.begin(), triplets.end());
  Eigen::SparseMatrix<float> mass(NumEqs, NumEqs);
  // WARNING: assumes the mass matrix is the identity
  mass.setIdentity();
  
  // Rayleigh damping matrix
  damping = alpha1 * stiffness + alpha2 * mass;
}

bool ExIntegrator::integrate(Clock& c) {
  
  size_t NumEqs = y.numCPs() * 3;
  VecXf dqdot = VecXf::Zero(NumEqs);
  
  for (YarnEnergy* e : energies) {
    if (!e->eval(VecXf::Zero(NumEqs), c, dqdot)) return false;
  }
  
  // Damping calculations
  VecXf vel(NumEqs);
  for (int i=0; i<y.numCPs(); i++) {
    vel.block<3,1>(3*i, 0) = y.cur().points[i].vel;
  }
  dqdot += damping * vel;
  
  for (int i=0; i<y.numCPs(); i++) {
    // WARNING: assumes the mass matrix is the identity
    Vec3f curdqdot = dqdot.block<3,1>(3*i, 0);
    y.next().points[i].accel = -curdqdot;
    y.next().points[i].vel = y.cur().points[i].vel - curdqdot;
    y.next().points[i].pos = y.cur().points[i].pos + c.timestep() * y.next().points[i].vel;
  }
  
  return true;
}