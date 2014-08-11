//
//  IMEXIntegrator.cpp
//  Visualizer
//
//  Created by eschweickart on 2/24/14.
//
//

#include "IMEXIntegrator.h"

DECLARE_DIFFSCALAR_BASE(); // Initialization of static struct
DECLARE_PROFILER();

IMEXIntegrator::IMEXIntegrator(std::vector<RodEnergy*>& energies, Rod& r) :
Integrator(r, energies) {
  // Fill hess base
  std::vector<Triplet> triplets;
  for (int i=1; i < r.numEdges()-1; i++) {
    if (i > 1) {
      triplets.push_back(Triplet(i-1, i-2, -2.0*r.getCS()[i].twistCoeff()/r.restVoronoiLength(i)));
    }
    if (i < r.numEdges()-2) {
      triplets.push_back(Triplet(i-1, i, -2.0*r.getCS()[i+1].twistCoeff()/r.restVoronoiLength(i+1)));
    }
  }
  hessBase = Eigen::SparseMatrix<real>(r.numEdges()-2, r.numEdges()-2);
  hessBase.setFromTriplets(triplets.begin(), triplets.end());
}

bool IMEXIntegrator::integrate(Clock& c) {
  PROFILER_START("Total");
  
  const real ConvergenceThreshold = 0.0125 * r.numCPs();
  const real ConvergenceTolerance = 1.25 * r.numCPs();
  
  bool evalSuccess = false;
  bool newtonConverge = false;
  int newtonIterations = 0;

  while (!evalSuccess) {
    VecXe Fx    = VecXe::Zero(r.numDOF());
    VecXe FxEx  = VecXe::Zero(r.numDOF());
    VecXe dqdot = VecXe::Zero(r.numDOF());
    VecXe sol   = VecXe::Zero(r.numDOF());
    Eigen::SparseMatrix<real> GradFx(r.numDOF(), r.numDOF());
    std::vector<Triplet> triplets;
    
    // TODO: Figure out a thread-safe way to fill this
    // TODO: Query energies to figure out a good estimate
    size_t numTriplets = 9*9*r.numIntCPs();
    
    // Calculate Fx contribution from Explicit energies
    // (these do not change between Newton iterations)
    for (RodEnergy* e : energies) {
      if (e->evalType() == Explicit) {
        evalSuccess = e->eval(&FxEx);
        CHECK_NAN_VEC(FxEx);
        if (!evalSuccess) break;
      }
    }
    if (!evalSuccess) continue;
    
    
    // Perform Newton iteration to solve IMEX equations
    while (!newtonConverge) {
      triplets.clear();
      triplets.reserve(numTriplets);
      Fx = FxEx;
      
      // Find offset for implicit evaluation
      VecXe offset = c.timestep() * (r.cur().vel + dqdot);
      
      // Add up energies
      for (RodEnergy* e : energies) {
        if (e->evalType() == Implicit) {
          evalSuccess = e->eval(&Fx, &triplets, &offset);
          if (!evalSuccess) break;
        }
      }
      if (!evalSuccess) break;
      
      Fx *= -c.timestep();
      Fx += r.getMass().sparse * dqdot;
      
      CHECK_NAN_VEC(Fx);
      
      // Test for convergence
      real residual = Fx.norm();
      if (residual < ConvergenceThreshold) {
        newtonConverge = true;
        break;
      } else if (newtonIterations > 4) {
        std::cout << "resid: " << residual << "\n";
        if (residual < ConvergenceTolerance) {
          newtonConverge = true;
        }
        break;
      }
      
      // Error too high; perform one Newton iteration to update dqdot
      newtonIterations++;
      
      // Create sparse Jacobian of Fx
      GradFx.setFromTriplets(triplets.begin(), triplets.end()); // sums up duplicates automagically
      GradFx *= -c.timestep() * c.timestep();
      GradFx += r.getMass().sparse;
      
      CHECK_NAN_VEC(GradFx.toDense());
      
      PROFILER_START("CG Solver");
      Eigen::ConjugateGradient<Eigen::SparseMatrix<real>, Eigen::Upper,
                               Eigen::IncompleteLUT<real>> cg;
      cg.compute(GradFx);
      assert(cg.preconditioner().info() == Eigen::ComputationInfo::Success);
      VecXe guess = sol;
      sol = cg.solveWithGuess(-Fx, guess); // H(x_n) (x_n+1 - x_n) = -F(x_n)
      
      if (cg.info() == Eigen::NoConvergence) {
        if (c.canDecreaseTimestep()) {
          PROFILER_STOP("CG Solver");
          c.suggestTimestep(c.timestep()/2.0);
          std::cout << "No convergence in CG solver. New timestep: " << c.timestep() << "\n";
          evalSuccess = false;
          break;
        }
        std::cerr << "No convergence!! Newton iterate: " << newtonIterations << "\n";
        std::cerr << "Fx has NaN: " << Fx.hasNaN() << "\n";
        std::cerr << "GradFx has NaN: " << GradFx.toDense().hasNaN() << "\n";
        std::cerr << "Fx max coeff: " << Fx.maxCoeff() << "\n";
        std::cerr << "GradFx max coeff: " << GradFx.toDense().maxCoeff() << "\n";
        assert(false);
      }
      
      dqdot += sol;
      PROFILER_STOP("CG Solver");
    }
    
    // Update rod positions
    if (newtonConverge) {
#ifdef NEWMARK_BETA
    // Newmark-Beta update
      const real gamma = 0.5;
      const real beta = 0.25;
      r.next().dVel = dqdot;
      r.next().vel = r.cur().vel + (1.0-gamma) * r.cur().dVel + gamma * dqdot;
      r.next().pos = r.cur().pos + c.timestep() * (r.cur().vel +
                                                   ((1.0-2.0*beta) / 2.0) * r.cur().dVel +
                                                   beta * dqdot);
    
#else // ifdef NEWMARK_BETA
      
      // Update changes to position and velocity
      r.next().vel = r.cur().vel + dqdot;
      r.next().pos = r.cur().pos + c.timestep() * r.next().vel;
      
#endif // ifdef NEWMARK_BETA
    }
  }
  
  PROFILER_STOP("Total");
  PROFILER_PRINT_ELAPSED();
  PROFILER_RESET_ALL();
  
  return newtonConverge;
}

void static calcRotEqs(const Rod& r, const VecXe& rot, const std::vector<Vec3e>& curveBinorm,
                      VecXe& grad, std::vector<Triplet>& triplets) {
  Eigen::Matrix<real, 2, 2> J;
  J << 0.0, -1.0, 1.0, 0.0;
  for (int i=1; i<r.numEdges()-1; i++) {
    Vec3e m1 = cos(rot(i)) * r.next().u[i] + sin(rot(i)) * r.next().v(i);
    Vec3e m2 = -sin(rot(i)) * r.next().u[i] + cos(rot(i)) * r.next().v(i);
    Vec2e curvePrev(curveBinorm[i-1].dot(m2), -curveBinorm[i-1].dot(m1)); // omega ^i _i
    Vec2e curveNext(curveBinorm[i].dot(m2), -curveBinorm[i].dot(m1)); // omega ^i _i+1
    real dWprev = 1.0 / r.restVoronoiLength(i) *
      curvePrev.dot(J * r.getCS()[i].bendMat() * (curvePrev - r.restCurveNext(i)));
    real dWnext = 1.0 / r.restVoronoiLength(i+1) *
      curveNext.dot(J * r.getCS()[i+1].bendMat() * (curveNext - r.restCurvePrev(i+1)));
    real twistPrev = rot(i) - rot(i-1) + r.next().refTwist(i);
    real twistNext = rot(i+1) - rot(i) + r.next().refTwist(i+1);
    grad(i-1) = -(dWprev + dWnext + 2.0 * r.getCS()[i].twistCoeff() *
                  (twistPrev/r.restVoronoiLength(i) - twistNext/r.restVoronoiLength(i+1)));
    
    real hess = 2.0*(r.getCS()[i].twistCoeff()/r.restVoronoiLength(i) +
                     r.getCS()[i+1].twistCoeff()/r.restVoronoiLength(i+1));
    hess += 1.0 / r.restVoronoiLength(i) *
      (curvePrev.dot(J.transpose() * r.getCS()[i].bendMat() * J * curvePrev)
       - curvePrev.dot(r.getCS()[i].bendMat() * (curvePrev - r.restCurveNext(i))));
    hess += 1.0 /r.restVoronoiLength(i+1) *
      (curveNext.dot(J.transpose() * r.getCS()[i+1].bendMat() * J * curveNext)
       - curveNext.dot(r.getCS()[i+1].bendMat() * (curveNext - r.restCurvePrev(i+1))));
    triplets.push_back(Triplet(i-1, i-1, hess));
  }
}

bool IMEXIntegrator::setRotations() const {
  const real newtonThreshold = 1.0e-5; //should be able to get this almost exact
  std::vector<Triplet> triplets;
  Eigen::SparseMatrix<real> hess(r.numEdges()-2, r.numEdges()-2);
  VecXe rot = r.next().rot;
  VecXe grad = VecXe::Zero(r.numEdges()-2); // Assumes edges are clamped
  bool newtonConverge = false;

  std::vector<Vec3e> curveBinorm;
  for (int i=1; i<r.numCPs()-1; i++) {
    Vec3e tPrev = r.next().edge(i-1).normalized();
    Vec3e tNext = r.next().edge(i).normalized();
    real chi = 1.0 + (tPrev.dot(tNext));
    curveBinorm.push_back(2.0*tPrev.cross(tNext)/chi);
  }
  int newtonIterations = 0;
  
  do {
    triplets.clear();
    calcRotEqs(r, rot, curveBinorm, grad, triplets);
    real resid = grad.norm();
    if (resid < newtonThreshold || newtonIterations > 4) {
      if (resid > 1.0e-5 * r.numEdges()) { return false; }
      newtonConverge = true;
      break;
    }
    newtonIterations++;
    hess.setFromTriplets(triplets.begin(), triplets.end());
    hess += hessBase;
    
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<real>> sLDLT;
    sLDLT.compute(hess);
    VecXe sol = sLDLT.solve(grad);
    assert(sLDLT.info() == Eigen::Success);
    rot.block(1, 0, r.numEdges()-2, 1) += sol;
  } while (!newtonConverge);
  
  if (newtonConverge) {
    r.next().rot = rot;
  }
  
  return newtonConverge;
}