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

IMEXIntegrator::IMEXIntegrator(std::vector<YarnEnergy*>& energies, Yarn& y) :
Integrator(y, energies) {
  // Fill hess base
  std::vector<Triplet> triplets;
  for (int i=1; i < y.numSegs()-1; i++) {
    if (i > 1) {
      triplets.push_back(Triplet(i-1, i-2, -2.0f*y.twistCoeff()/y.restVoronoiLength(i)));
    }
    if (i < y.numSegs()-2) {
      triplets.push_back(Triplet(i-1, i, -2.0f*y.twistCoeff()/y.restVoronoiLength(i+1)));
    }
  }
  hessBase = Eigen::SparseMatrix<float>(y.numSegs()-2, y.numSegs()-2);
  hessBase.setFromTriplets(triplets.begin(), triplets.end());
}

bool IMEXIntegrator::integrate(Clock& c) {
  Profiler::start("Total");
  
  const float ConvergenceThreshold = 0.0125f * y.numCPs();
  const float ConvergenceTolerance = 1.25f * y.numCPs();
  
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
    std::vector<Triplet> triplets;
    
    // TODO: Figure out a thread-safe way to fill this
    // TODO: Query energies to figure out a good estimate
    size_t numTriplets = 9*9*y.numIntCPs();
    
    // Calculate Fx contribution from Explicit energies
    // (these do not change between Newton iterations)
    for (YarnEnergy* e : energies) {
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
      VecXf offset(NumEqs);
      for (int i=0; i<y.numCPs(); i++) {
        offset.block<3,1>(3*i, 0) = c.timestep()*(y.cur().points[i].vel + dqdot.block<3,1>(3*i, 0));
      }
      
      // Add up energies
      for (YarnEnergy* e : energies) {
        if (e->evalType() == Implicit) {
          evalSuccess = e->eval(&Fx, &triplets, &offset);
          if (!evalSuccess) break;
        }
      }
      if (!evalSuccess) break;
      
      // WARNING: assumes the mass matrix is the identity
      for (int i=0; i<NumEqs; i++) {
        Fx(i) = dqdot(i) + (-c.timestep()*Fx(i));
//        triplets.push_back(Triplet(i, i, -1.0f / c.timestep() / c.timestep()));
      }
      
      CHECK_NAN_VEC(Fx);
      
      // Test for convergence
      float residual = Fx.norm();
      if (residual < ConvergenceThreshold) {
        newtonConverge = true;
        break;
      } else if (newtonIterations > 4) {
        std::cout << "resid: " << residual << "\n";
        if (residual < ConvergenceTolerance) {
          newtonConverge = true;
        }
#ifdef DRAW_IMEX_INTEGRATOR
        float maxcoeff = error.maxCoeff();
        Yarn* yp = &y;
        for (int i=0; i<y.numCPs(); i++) {
          Vec3f curerror = error.block<3,1>(3*i, 0);
          frames.push_back([yp, i, maxcoeff, curerror] () {
            ci::gl::color(curerror[0]/maxcoeff, curerror[1]/maxcoeff, curerror[2]/maxcoeff, 0.7f);
            ci::gl::drawSphere(toCi(yp->cur().points[i].pos), constants::radius*2.0f);
          });
        }
#endif // ifdef DRAW_IMEX_INTEGRATOR
        break;
      }
      
      // Error too high; perform one Newton iteration to update dqdot
      newtonIterations++;
      
      // Create sparse Jacobian of Fx
      GradFx.setFromTriplets(triplets.begin(), triplets.end()); // sums up duplicates automagically
      GradFx *= -c.timestep() * c.timestep();
      // WARNING: assumes the mass matrix is the identity
      Eigen::SparseMatrix<float> id(NumEqs, NumEqs);
      id.setIdentity();
      GradFx += id;
      
      CHECK_NAN_VEC(GradFx.toDense());
      
      Profiler::start("CG Solver");
      Eigen::ConjugateGradient<Eigen::SparseMatrix<float>, Eigen::Upper,
                               Eigen::IncompleteLUT<float>> cg;
      cg.compute(GradFx);
      Eigen::VectorXf guess = sol;
      sol = cg.solveWithGuess(-Fx, guess); // H(x_n) (x_n+1 - x_n) = -F(x_n)
      
      if (cg.info() == Eigen::NoConvergence) {
        if (c.canDecreaseTimestep()) {
          Profiler::stop("CG Solver");
          c.suggestTimestep(c.timestep()/2.0f);
          std::cout << "No convergence in CG solver. New timestep: " << c.timestep() << "\n";
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
      
      dqdot += sol;
      Profiler::stop("CG Solver");
    }
    
    // Update yarn positions
    if (newtonConverge) {
#ifdef NEWMARK_BETA
    // Newmark-Beta update
      const float gamma = 0.5f;
      const float beta = 0.25f;
      for (int i=0; i<y.numCPs(); i++) {
        Vec3f curdqdot = dqdot.block<3, 1>(3*i, 0);
        y.next().points[i].vel = y.cur().points[i].vel +
                                 (1.0f-gamma)*y.cur().points[i].accel + gamma*curdqdot;
        y.next().points[i].pos = y.cur().points[i].pos +
                                 c.timestep()*(y.cur().points[i].vel +
                                               (1.0f-2.0f*beta)/2.0f*y.cur().points[i].accel +
                                               beta*curdqdot);
        y.next().points[i].accel = curdqdot;
      }
    
#else // ifdef NEWMARK_BETA
    
      // Update changes to position and velocity
      for (int i=0; i<y.numCPs(); i++) {
        Vec3f curdqdot = dqdot.block<3, 1>(3*i, 0);
        y.next().points[i].vel = y.cur().points[i].vel + curdqdot;
        y.next().points[i].pos = y.cur().points[i].pos + c.timestep()*y.next().points[i].vel;
      }
      
#endif // ifdef NEWMARK_BETA
    }
  }
  
  Profiler::stop("Total");
//  Profiler::printElapsed();
//  std::cout << "\n";
  Profiler::resetAll();
  
  return newtonConverge;
}

void static calcRotEqs(const Yarn& y, const VecXf& rot, const std::vector<Vec3f>& curveBinorm,
                      VecXf& grad, std::vector<Triplet>& triplets) {
  Eigen::Matrix2f J;
  J << 0.0f, -1.0f, 1.0f, 0.0f;
  // This assumes the yarn is isotropic
  for (int i=1; i<y.numSegs()-1; i++) {
    const Segment& s = y.next().segments[i];
    Vec3f m1 = cos(rot(i)) * s.getU() + sin(rot(i)) * s.v();
    Vec3f m2 = -sin(rot(i)) * s.getU() + cos(rot(i)) * s.v();
    Vec2f curvePrev(curveBinorm[i-1].dot(m2), -curveBinorm[i-1].dot(m1)); // omega ^i _i
    Vec2f curveNext(curveBinorm[i].dot(m2), -curveBinorm[i].dot(m1)); // omega ^i _i+1
    float dWprev = y.bendCoeff() / y.restVoronoiLength(i) *
      curvePrev.dot(J * (curvePrev - y.restCurveNext(i)));
    float dWnext = y.bendCoeff() / y.restVoronoiLength(i+1) *
      curveNext.dot(J * (curveNext - y.restCurvePrev(i+1)));
    float twistPrev = rot(i) - rot(i-1) + y.next().segments[i].getRefTwist();
    float twistNext = rot(i+1) - rot(i) + y.next().segments[i+1].getRefTwist();
    grad(i-1) = -(dWprev + dWnext + 2.0f*y.twistCoeff()*
                  (twistPrev/y.restVoronoiLength(i) - twistNext/y.restVoronoiLength(i+1)));
    
    float hess = 2.0f*(y.twistCoeff()/y.restVoronoiLength(i) +
                       y.twistCoeff()/y.restVoronoiLength(i+1));
    hess += y.bendCoeff()/y.restVoronoiLength(i) *
      (curvePrev.dot(curvePrev) - curvePrev.dot(curvePrev - y.restCurveNext(i)));
    hess += y.bendCoeff()/y.restVoronoiLength(i+1) *
      (curveNext.dot(curveNext) - curveNext.dot(curveNext - y.restCurvePrev(i+1)));
    triplets.push_back(Triplet(i-1, i-1, hess));
  }
}

bool IMEXIntegrator::setRotations() const {
  const float newtonThreshold = 1e-5; //should be able to get this almost exact
  std::vector<Triplet> triplets;
  Eigen::SparseMatrix<float> hess(y.numSegs()-2, y.numSegs()-2);
  VecXf rot(y.numSegs());
  VecXf grad = VecXf::Zero(y.numSegs()-2); // Assumes edges are clamped
  bool newtonConverge = false;
  for(int i=0; i<y.numSegs(); i++) {
    rot(i) = y.next().segments[i].getRot();
  }
  std::vector<Vec3f> curveBinorm;
  for (int i=1; i<y.numCPs()-1; i++) {
    Vec3f tPrev = y.next().segments[i-1].vec().normalized();
    Vec3f tNext = y.next().segments[i].vec().normalized();
    float chi = 1.0f + (tPrev.dot(tNext));
    curveBinorm.push_back(2.0f*tPrev.cross(tNext)/chi);
  }
  int newtonIterations = 0;
  
  // TESTING
  // VecXf guess = VecXf::Zero(y.numSegs()-2);
  // float totalTwist = y.next().segments[y.numSegs()-1].getRot() - y.next().segments[0].getRot();
  
  
  do {
    triplets.clear();
    calcRotEqs(y, rot, curveBinorm, grad, triplets);
    float resid = grad.norm();
    if (resid < newtonThreshold || newtonIterations > 4) {
      if (resid > 1e-5 * y.numSegs()) { return false; }
      newtonConverge = true;
      break;
    }
    newtonIterations++;
    hess.setFromTriplets(triplets.begin(), triplets.end());
    hess += hessBase;
    
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>> sLDLT;
    sLDLT.compute(hess);
    VecXf sol = sLDLT.solve(grad);
    assert(sLDLT.info() == Eigen::Success);
    rot.block(1, 0, y.numSegs()-2, 1) += sol;
  } while (!newtonConverge);
  
  if (newtonConverge) {
    for (int i=1; i<y.numSegs()-1; i++) {
      y.next().segments[i].setRot(rot(i));
    }
  }
  
  return newtonConverge;
}