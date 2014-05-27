//
//  Integrator.cpp
//  Visualizer
//
//  Created by eschweickart on 2/24/14.
//
//

#include "Integrator.h"
#include <boost/timer.hpp>

const float ConvergenceThreshold = 0.5;
const float ConvergenceTolerance = 50;

DECLARE_DIFFSCALAR_BASE(); // Initialization of static struct
DECLARE_PROFILER();

bool Integrator::integrate(Yarn& y, Clock& c) {
  frames.clear();
  
  Profiler::start("Total");
  
  // Compute timestep
  for (YarnEnergy* e : energies) {
    e->suggestTimestep(c);
  }
  
  const size_t NumEqs = y.numCPs() * 3;
  bool evalSuccess = false;
  bool newtonConverge = false;
  int newtonIterations = 0;

  while (!evalSuccess) {
    Eigen::VectorXf Fx    = Eigen::VectorXf::Zero(NumEqs);
    Eigen::VectorXf FxEx  = Eigen::VectorXf::Zero(NumEqs);
    Eigen::VectorXf dqdot = Eigen::VectorXf::Zero(NumEqs);
    Eigen::VectorXf sol   = Eigen::VectorXf::Zero(NumEqs);
    Eigen::SparseMatrix<float> GradFx(NumEqs, NumEqs);
    Eigen::SparseMatrix<float> id(NumEqs, NumEqs);
    id.setIdentity();
    std::vector<Triplet> triplets;
    
    // TODO: Figure out a thread-safe way to fill this
    // TODO: Query energies to figure out a good estimate
    size_t numTriplets = 9*9*y.numIntCPs();
    
    // Calculate Fx contribution from Explicit energies
    // (these do not change between Newton iterations)
    for (YarnEnergy* e : energies) {
      if (e->evalType() == Explicit) {
        evalSuccess = e->eval(dqdot, c, FxEx); // Ignores dqdot
        CHECK_NAN_VEC(FxEx);
        if (!evalSuccess) break;
      }
    }
    if (!evalSuccess) continue;
    
    
    // Perform Newton iteration to solve IMEX equations
    while (!newtonConverge) {
      if (newtonIterations != 0) {
        triplets.clear();
        Fx.setZero();
      }
      
      triplets.reserve(numTriplets);
      newtonIterations++;
      
      // Add up energies
      for (YarnEnergy* e : energies) {
        if (e->evalType() == Implicit) {
          evalSuccess = e->eval(dqdot, c, Fx, &triplets);
          if (!evalSuccess) break;
        }
      }
      if (!evalSuccess) break;
      
      // TODO: Mass matrix may not be I
      for (int i=0; i<NumEqs; i++) {
        Fx(i) += dqdot(i) + FxEx(i);
        triplets.push_back(Triplet(i, i, 1));
      }
      
      // Solve equations for updates to changes in position and velocity using Conjugate Gradient
      GradFx.setFromTriplets(triplets.begin(), triplets.end()); // sums up duplicates automagically
//      GradFx += id;
      
      CHECK_NAN_VEC(Fx);
      CHECK_NAN_VEC(GradFx.toDense());
      
      Profiler::start("CG Solver");
      Eigen::ConjugateGradient<Eigen::SparseMatrix<float>, Eigen::Upper, Eigen::IncompleteLUT<float>> cg;
      cg.compute(GradFx);
      Eigen::VectorXf temp = sol;
      sol = cg.solveWithGuess(Fx, temp);
      
      if (cg.info() == Eigen::NoConvergence) {
        if (c.canDecreaseTimestep()) {
          Profiler::stop("CG Solver");
          c.suggestTimestep(c.timestep()/2);
          std::cout << "Warning: No convergence in CG solver. New timestep: " << c.timestep() << "\n";
          std::cout << "GradFx max coeff: " << GradFx.toDense().maxCoeff() << "\n";
          evalSuccess = false;
          break;
        }
        std::cerr << "No convergence!! Newton iterate: " << newtonIterations << "\n";
        std::cerr << "Fx all finite: " << Fx.allFinite() << "\n";
        std::cerr << "GradFx all finite: " << GradFx.toDense().allFinite() << "\n";
        std::cerr << "Fx max coeff: " << Fx.maxCoeff() << "\n";
        std::cerr << "GradFx max coeff: " << GradFx.toDense().maxCoeff() << "\n";
        assert(false);
      }
      
      dqdot -= sol;
      Profiler::stop("CG Solver");
      
      VecXf FxIm = VecXf::Zero(NumEqs);
      for (YarnEnergy* e : energies) {
        if (e->evalType() == Implicit) {
          e->eval(dqdot, c, FxIm);
        }
      }
      
      
      VecXf error = dqdot + FxIm + FxEx;
      float residual = error.norm();
      if (residual < ConvergenceThreshold) {
        newtonConverge = true;
      } else if (newtonIterations > 4) {
        std::cerr << "resid: " << residual << "\n";
        
        if (residual < ConvergenceTolerance) newtonConverge = true;
        
        float maxcoeff = error.maxCoeff();
        Yarn* yp = &y;
        for (int i=0; i<y.numCPs(); i++) {
          Vec3f curerror = error.block<3,1>(3*i, 0);
          frames.push_back([yp, i, maxcoeff, curerror] () {
            ci::gl::color(curerror[0]/maxcoeff, curerror[1]/maxcoeff, curerror[2]/maxcoeff, 0.7);
            ci::gl::drawSphere(toCi(yp->cur().points[i].pos), constants::radius*2);
          });
        }
        
        break;
      }
    }
    if (!evalSuccess) continue;
    
    // Update yarn positions
    if (newtonConverge) {
// #define NEWMARK_BETA
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
      y.next().points[i].vel = y.cur().points[i].vel + curdqdot;
      y.next().points[i].pos = c.timestep()*y.next().points[i].vel + y.cur().points[i].pos;
    }
    
#endif // ifdef NEWMARK_BETA
    }
  }
  
  Profiler::stop("Total");
  //  Profiler::printElapsed();
  Profiler::resetAll();
  //  std::cout << "\n";
  
  return newtonConverge;
}

void const Integrator::draw() {
  for (std::function<void(void)> f : frames) {
    f();
  }
}