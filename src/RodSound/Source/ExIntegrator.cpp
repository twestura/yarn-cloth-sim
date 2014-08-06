//
//  ExIntegrator.cpp
//  Visualizer
//
//  Created by eschweickart on 6/10/14.
//
//

#include "ExIntegrator.h"

ExIntegrator::ExIntegrator(Yarn& y, std::vector<YarnEnergy*>& energies,
                           std::vector<YarnConstraint*>* c) : Integrator(y, energies) {
  // Stiffness coefficient
  alpha1 = 1.0e-9;
  // Mass coefficient
  alpha2 = 2.0;
  
  size_t dof = y.numCPs() * 3;
  stiffness.resize(dof, dof);
  setDamping();
  
#ifdef DRAW_EIGENMODE
  // Test eigenvalues to extract frequencies.
  Eigen::SparseMatrix<real> lklt = y.getInvMass().sparse * stiffness;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>> saes(lklt.toDense());
  for (int i=0; i<dof; i++) {
    real eigval = saes.eigenvalues()(i);
    if (eigval < 0) continue; // Only positive eigenvalues contribute.
    real temp = alpha1 * eigval + alpha2;
    real disc = temp * temp - 4.0 * eigval;
    if (disc >= 0) continue; // No imaginary part, which specifies the frequency.
    typedef std::complex<real> comp;
    comp freqi = std::sqrt(comp(disc)) / 2.0;
    real freq = fabs(freqi.imag()) / (2.0 * constants::pi);
    if (freq > 20.0 && freq < 20000.0) {
      std::cout << "eigval " << i << " (" << eigval << "): " << "Freq: " << freq << "\n";
    }
  }
  std::cout << "\n";
  modes = saes.eigenvectors();
#endif // ifdef DRAW_EIGENMODE
  
  if (c) constraints = *c;
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
  
  // Constraint-based integration
  /*
  // Find candidate positions
  VecXe xStar = VecXe(dof);
  for (int i=0; i<y.numCPs(); i++) {
    Vec3e velStar = y.cur().points[i].vel + c.timestep()*forces.segment<3>(3*i);
    xStar.segment<3>(3*i) = y.cur().points[i].pos + c.timestep()*velStar;
  }
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
    y.next().points[i].vel = delta / c.timestep();
    y.next().points[i].pos = xStar.segment<3>(3*i);
    y.next().points[i].accel = y.next().points[i].vel - y.cur().points[i].vel;
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
  stiffness *= -1.0;
  
  // Rayleigh damping matrix
  damping = alpha1 * stiffness + alpha2 * y.getMass().sparse;
}


void ExIntegrator::draw() {
#ifdef DRAW_EIGENMODE
  size_t eigval = 56;
  VecXe mode = modes.col(eigval);
  
  ci::gl::color(0.7, 0.1, 0.6);
  for (int i=0; i<y.numCPs(); i++) {
    Vec3e flux = mode.segment<3>(i*3);
    ci::gl::drawLine(EtoC(y.cur().points[i].pos), EtoC(y.cur().points[i].pos + flux));
  }
#endif DRAW_EIGENMODE
}