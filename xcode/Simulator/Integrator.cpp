//
//  Integrator.cpp
//  Visualizer
//
//  Created by eschweickart on 2/24/14.
//
//

#include "Integrator.h"

const float ConvergenceThreshold = 0.001;

DECLARE_DIFFSCALAR_BASE(); // Initialization of static struct

void Integrator::integrate(Yarn& y, Clock& c) {
  
  const size_t NumEqs = y.numCPs() * 3;
  
  Eigen::VectorXf Fx    = Eigen::VectorXf::Zero(NumEqs);
  Eigen::VectorXf dqdot = Eigen::VectorXf::Zero(NumEqs);
  Eigen::VectorXf sol   = Eigen::VectorXf::Zero(NumEqs);
  Eigen::SparseMatrix<float> GradFx(NumEqs, NumEqs);
  std::vector<Triplet> triplets;
  
  // WARNING: push_back() on a vector is thread-safe if no allocation is performed. Make sure
  // the correct amount of space is reserved here!
  size_t numTriplets = 9*9*y.numIntCPs();
  
  bool converge = false;
  int iterations = 0;
  
  for (int i=0; i<y.numCPs(); i++) {
#ifdef ENABLE_CHECK_NAN
    assert(y.cur().points[i].pos.allFinite());
    assert(y.cur().points[i].vel.allFinite());
#endif
    y.next().points[i].pos = y.cur().points[i].pos;
    y.next().points[i].vel = y.cur().points[i].vel;
  }
  
  // Perform Newton iteration to solve IMEX equations
  while (!converge) {
    if (iterations != 0) {
      triplets.clear();
      Fx.setZero();
    }
    
    triplets.reserve(numTriplets);
    iterations++;
  
    // Add up energies
    for (YarnEnergy* e : energies) {
      e->eval(Fx, triplets, dqdot, c);
    }
    
    // TODO: Mass matrix may not be I
    for (int i=0; i<NumEqs; i++) {
      Fx(i) += dqdot(i);
      triplets.push_back(Triplet(i, i, 1));
    }
    
    // Solve equations for updates to changes in position and velocity using Conjugate Gradient
    GradFx.setFromTriplets(triplets.begin(), triplets.end()); // sums up duplicates automagically
    
#ifdef ENABLE_CHECK_NANS
    assert(Fx.allFinite());
#endif
    
    Eigen::ConjugateGradient<Eigen::SparseMatrix<float>> cg;
    cg.compute(GradFx);
    Eigen::VectorXf temp = sol;
    sol = cg.solveWithGuess(Fx, temp);
    
    if (cg.info() == Eigen::NoConvergence) {
      std::cerr << "No convergence!! Newton iterate: " << iterations << "\n";
      std::cerr << "Fx max coeff: " << Fx.maxCoeff() << "\n";
      std::cerr << "GradFx max coeff: " << GradFx.toDense().maxCoeff() << "\n";
      assert(false);
    } else {
//      std::cout << "CG iters: " << cg.iterations() << "\n";
    }
    
    dqdot -= sol;
    
    if (sol.maxCoeff() < ConvergenceThreshold) {
      converge = true;
    } else if (iterations > 5) {
      std::cerr << "Too many newton iterations, breaking.\n";
      break;
    }
  }
  
  // Update changes to position and velocity
  for (int i=0; i<y.numCPs(); i++) {
    Vec3f curdqdot = dqdot.block<3, 1>(3*i, 0);
    y.next().points[i].vel += curdqdot;
    y.next().points[i].pos += c.timestep()*y.next().points[i].vel;
  }
}