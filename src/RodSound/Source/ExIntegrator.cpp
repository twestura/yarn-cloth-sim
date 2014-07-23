//
//  ExIntegrator.cpp
//  Visualizer
//
//  Created by eschweickart on 6/10/14.
//
//

#include "ExIntegrator.h"

ExIntegrator::ExIntegrator(Yarn& y, std::vector<YarnEnergy*>& energies) : Integrator(y, energies) {
  alpha1 = 1.0e-8;
  alpha2 = 0.0;
  
  size_t dof = y.numCPs() * 3;
  stiffness.resize(dof, dof);
  setDamping();
  
  Eigen::SparseMatrix<real> lklt = y.getInvMass().sparse * -stiffness;
  
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>> saes(lklt.toDense());
  std::cout << saes.eigenvalues() << "\n\n";
  
  for (int i=0; i<4; i++) {
    real eigval = saes.eigenvalues()(dof-4+i);
    std::cout << "eigval " << dof-4+i << " (" << eigval << "): ";
    if (eigval < 0) {
      std::cout << "eigval is negative, skipping.\n";
      continue;
    }
    real temp = alpha1 * eigval + alpha2;
    real disc = temp * temp - 4.0 * eigval;
    if (disc >= 0) {
      std::cout << "disc is non-negative: " << disc << "\n";
      continue;
    }
    typedef std::complex<real> comp;
    comp freqi = std::sqrt(comp(disc)) / 2.0; // only care about imaginary part
    // comp freq2i = (-temp - std::sqrt(comp(disc))) / 2.0;
    real freq = fabs(freqi.imag()) / (2.0 * constants::pi);
    // real freq2 = fabs(freq2i.imag()) / (2.0 * constants::pi);
    std::cout << "Freq: " << freq << "\n";
  }
  std::cout << "\n";
  
  // std::cout << saes.eigenvectors().topRightCorner(dof, 4) << "\n\n";
}

bool ExIntegrator::integrate(Clock& c) {
  PROFILER_START("Integrate");
  
  size_t dof = y.numCPs() * 3;
  VecXe forces = VecXe::Zero(dof);
  
  for (YarnEnergy* e : energies) {
    if (!e->eval(&forces)) return false;
  }
  
  // Damping calculations
  VecXe vel(dof);
  for (int i=0; i<y.numCPs(); i++) {
    vel.segment<3>(3*i) = y.cur().points[i].vel;
  }
  if (c.getTicks() % 1000 == 0 && c.getTicks() != 0) { // Periodically update stiffness matrix
    setDamping();
  }
  forces -= damping * vel;
  
  forces = y.getInvMass().sparse * forces;
  
  // Symplectic Euler
  for (int i=0; i<y.numCPs(); i++) {
    Vec3e dqdot = forces.segment<3>(3*i) * c.timestep();
    y.next().points[i].accel = dqdot;
    y.next().points[i].vel = y.cur().points[i].vel + dqdot;
    y.next().points[i].pos = y.cur().points[i].pos + y.next().points[i].vel * c.timestep();
  }
  
  // Explicit Newmark -- NOT STABLE
  /*
  real lambda = 0.5;
  for (int i=0; i<y.numCPs(); i++) {
    Vec3e dqdot = forces.segment<3>(3*i) * c.timestep();
    y.next().points[i].accel = dqdot;
    y.next().points[i].vel = y.cur().points[i].vel + (1.0-lambda) * y.cur().points[i].accel + lambda * dqdot;
    y.next().points[i].pos = y.cur().points[i].pos + c.timestep() * y.cur().points[i].vel + 0.5 * c.timestep() * y.cur().points[i].accel;
  }
   */
  
  
  PROFILER_STOP("Integrate");
  return true;
}

void ExIntegrator::setDamping() {
  std::vector<Triplet> triplets;
  for (YarnEnergy* e : energies) {
    if (e->energySource() == Internal) {
      e->eval(nullptr, &triplets);
    }
  }
  stiffness.setFromTriplets(triplets.begin(), triplets.end());
  stiffness *= -1;
  
  // Rayleigh damping matrix
  damping = alpha1 * stiffness + alpha2 * y.getMass().sparse;
}