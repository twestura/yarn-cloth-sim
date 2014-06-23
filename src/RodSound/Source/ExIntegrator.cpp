//
//  ExIntegrator.cpp
//  Visualizer
//
//  Created by eschweickart on 6/10/14.
//
//

#include "ExIntegrator.h"

ExIntegrator::ExIntegrator(Yarn& y, std::vector<YarnEnergy*>& energies) : Integrator(y, energies) {
  const real alpha1 = 6.615e-7;
  const real alpha2 = 0.882;
  
  std::vector<Triplet> triplets;
  size_t NumEqs = y.numCPs() * 3;
  for (YarnEnergy* e : energies) {
    if (e->energySource() == Internal) {
      e->eval(nullptr, &triplets);
    }
  }
  
  Eigen::SparseMatrix<real> stiffness(NumEqs, NumEqs);
  stiffness.setFromTriplets(triplets.begin(), triplets.end()); // Remember to negate this!
  
  // Eigen::SelfAdjointEigenSolver<Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>> saes(-stiffness.toDense());
  // std::cout << saes.eigenvalues() << "\n\n";
  
  // Rayleigh damping matrix
  damping = -alpha1 * stiffness + alpha2 * y.getMass().sparse;
}

bool ExIntegrator::integrate(Clock& c) {
  
  size_t NumEqs = y.numCPs() * 3;
  VecXe forces = VecXe::Zero(NumEqs);
  
  for (YarnEnergy* e : energies) {
    if (!e->eval(&forces)) return false;
  }
  
  // Damping calculations
  VecXe vel(NumEqs);
  for (int i=0; i<y.numCPs(); i++) {
    vel.block<3,1>(3*i, 0) = y.cur().points[i].vel;
  }
  forces -= damping * vel;
  forces = y.getInvMass().sparse * forces;
  
  for (int i=0; i<y.numCPs(); i++) {
    Vec3e dqdot = forces.block<3,1>(3*i, 0) * c.timestep();
    y.next().points[i].accel = dqdot;
    y.next().points[i].vel = y.cur().points[i].vel + dqdot;
    y.next().points[i].pos = y.cur().points[i].pos + y.next().points[i].vel * c.timestep();
  }
  
  return true;
}