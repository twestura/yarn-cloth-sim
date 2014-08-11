//
//  ExIntegrator.cpp
//  Visualizer
//
//  Created by eschweickart on 6/10/14.
//
//

#include "ExIntegrator.h"

ExIntegrator::ExIntegrator(Rod& r, std::vector<RodEnergy*>& energies,
                           std::vector<RodConstraint*>* c) : Integrator(r, energies) {
  // Stiffness coefficient
  alpha1 = 1.0e-9;
  // Mass coefficient
  alpha2 = 2.0;
  
  stiffness.resize(r.numDOF(), r.numDOF());
  setDamping();
  
#ifdef DRAW_EIGENMODE
  // Test eigenvalues to extract frequencies.
  Eigen::SparseMatrix<real> lklt = r.getInvMass().sparse * stiffness;
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
  
  size_t dof = r.numCPs() * 3;
  VecXe forces = VecXe::Zero(dof);
  
  for (RodEnergy* e : energies) {
    if (!e->eval(&forces)) return false;
  }
  
  // Damping calculations
  if (c.getTicks() % 1000 == 0 && c.getTicks() != 0) { // Periodically update stiffness matrix
    setDamping();
  }
  forces -= damping * r.cur().vel;
  
  forces = r.getInvMass().sparse * forces;
  
  // Symplectic Euler
  r.next().dVel = forces * c.timestep();
  r.next().vel = r.cur().vel + r.next().dVel;
  r.next().pos = r.cur().pos + r.next().vel * c.timestep();
  
  // Explicit Newmark -- NOT STABLE
  /*
  real lambda = 0.5;
   
  r.next().dVel = forces * c.timestep();
  r.next().vel = r.cur().vel + (1.0-lambda) * r.cur().dVel + lambda * r.next().dVel;
  r.next().pos = r.cur().pos + c.timestep() * r.cur().vel + 0.5 * c.timestep() * r.cur().dVel;
   */
  
  // Constraint-based integration
  /*
  // Find candidate positions
  VecXe xStar = r.cur().pos + c.timestep() * (r.cur().vel + c.timestep() * forces)

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
  r.next().pos = xStar;
  r.next().dVel = r.next().vel - r.cur().vel
   */
  
  
  PROFILER_STOP("Integrate");
  return true;
}

void ExIntegrator::setDamping() {
  std::vector<Triplet> triplets;
  for (RodEnergy* e : energies) {
    if (e->energySource() == Internal) {
      e->eval(nullptr, &triplets);
    }
  }
  stiffness.setFromTriplets(triplets.begin(), triplets.end());
  stiffness *= -1.0;
  
  // Rayleigh damping matrix
  damping = alpha1 * stiffness + alpha2 * r.getMass().sparse;
}


void ExIntegrator::draw() {
#ifdef DRAW_EIGENMODE
  size_t eigval = 56;
  VecXe mode = modes.col(eigval);
  
  ci::gl::color(0.7, 0.1, 0.6);
  for (int i=0; i<r.numCPs(); i++) {
    Vec3e flux = mode.segment<3>(i*3);
    ci::gl::drawLine(EtoC(r.cur().points[i].pos), EtoC(r.cur().points[i].pos + flux));
  }
#endif // ifdef DRAW_EIGENMODE
}