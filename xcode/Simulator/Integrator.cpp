//
//  Integrator.cpp
//  Visualizer
//
//  Created by eschweickart on 2/24/14.
//
//

#include "Integrator.h"
#include <boost/timer.hpp>

const float ConvergenceThreshold = 0.001;

DECLARE_DIFFSCALAR_BASE(); // Initialization of static struct
DECLARE_PROFILER();

void Integrator::integrate(Yarn& y, Clock& c) {
  
  Profiler::start("Total");
  
  // Compute timestep
  for (YarnEnergy* e : energies) {
    e->suggestTimestep(c);
  }
  
  const size_t NumEqs = y.numCPs() * 3;
  bool success = false;
  int iter = 0;

  while (!success) {
    iter++;
//    assert(iter < 10 && "Could not find a stable timestep!!");
    
    Eigen::VectorXf Fx    = Eigen::VectorXf::Zero(NumEqs);
    Eigen::VectorXf FxEx  = Eigen::VectorXf::Zero(NumEqs);
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

    
    // Calculate Fx contribution from Explicit energies
    // (these do not change between Newton iterations)
    for (YarnEnergy* e : energies) {
      if (e->evalType() == Explicit) {
        success = e->eval(FxEx, triplets, dqdot, c); // Ignores triplets and dqdot
        if (!success) break;
      }
    }
    if (!success) continue;
    
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
        if (e->evalType() == Implicit) {
          success = e->eval(Fx, triplets, dqdot, c);
          if (!success) break;
        }
      }
      if (!success) break;
      
      // TODO: Mass matrix may not be I
      for (int i=0; i<NumEqs; i++) {
        Fx(i) += dqdot(i) + FxEx(i);
        triplets.push_back(Triplet(i, i, 1));
      }
      
      // Solve equations for updates to changes in position and velocity using Conjugate Gradient
      GradFx.setFromTriplets(triplets.begin(), triplets.end()); // sums up duplicates automagically
      
#ifdef ENABLE_CHECK_NANS
      assert(Fx.allFinite());
#endif
      
      Profiler::start("CG Solver");
      Eigen::ConjugateGradient<Eigen::SparseMatrix<float>, Eigen::Lower, Eigen::IncompleteLUT<float>> cg;
      cg.compute(GradFx);
      Eigen::VectorXf temp = sol;
      sol = cg.solveWithGuess(Fx, temp);
      
      if (cg.info() == Eigen::NoConvergence) {
        if (c.canDecreaseTimestep()) {
          Profiler::stop("CG Solver");
          c.suggestTimestep(c.timestep()/2);
          std::cout << "Warning: No convergence in CG solver. New timestep: " << c.timestep() << "\n";
          std::cout << "GradFx max coeff: " << GradFx.toDense().maxCoeff() << "\n";
          success = false;
          break;
        }
        std::cerr << "No convergence!! Newton iterate: " << iterations << "\n";
        std::cerr << "Fx all finite: " << Fx.allFinite() << "\n";
        std::cerr << "GradFx all finite: " << GradFx.toDense().allFinite() << "\n";
        std::cerr << "Fx max coeff: " << Fx.maxCoeff() << "\n";
        std::cerr << "GradFx max coeff: " << GradFx.toDense().maxCoeff() << "\n";
        assert(false);
      } else {
        //      std::cout << "CG iters: " << cg.iterations() << "\n";
      }
      
      dqdot -= sol;
      Profiler::stop("CG Solver");
      
      
      if (sol.maxCoeff() < ConvergenceThreshold) {
        converge = true;
      } else if (iterations > 3) {
//        std::cerr << "Too many newton iterations, breaking.\n";
        break;
      }
    }
    if (!success) continue;
    
    //#define NEWMARK_BETA
#ifdef NEWMARK_BETA
    
    // Newmark-Beta update
    const float gamma = 0.5;
    const float beta = 0.25;
    for (int i=0; i<y.numCPs(); i++) {
      Vec3f curdqdot = dqdot.block<3, 1>(3*i, 0);
      y.next().points[i].vel = y.cur().points[i].vel + (1-gamma)*y.cur().points[i].accel + gamma*curdqdot;
      y.next().points[i].pos = y.cur().points[i].pos + c.timestep()*(y.cur().points[i].vel +
                                                                     (1-2*beta)/2*y.cur().points[i].accel +
                                                                     beta*curdqdot);
      y.next().points[i].accel = curdqdot;
    }
    
#else // ifdef NEWMARK_BETA
    
    // Update changes to position and velocity
    for (int i=0; i<y.numCPs(); i++) {
      Vec3f curdqdot = dqdot.block<3, 1>(3*i, 0);
      y.next().points[i].vel += curdqdot;
      y.next().points[i].pos += c.timestep()*y.next().points[i].vel;
    }
    
#endif // ifdef NEWMARK_BETA
  }
  
  Profiler::stop("Total");
  //  Profiler::printElapsed();
  Profiler::resetAll();
  //  std::cout << "\n";
}