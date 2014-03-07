//
//  Integrator.cpp
//  Visualizer
//
//  Created by eschweickart on 2/24/14.
//
//

#include "Integrator.h"
#include <stdexcept>

#define NUM_VARS 9
const float ConvergenceThreshold = 0.001;

typedef Eigen::Matrix<float, NUM_VARS, 1> Gradient;
typedef Eigen::Matrix<float, NUM_VARS, NUM_VARS> Hessian;
typedef DScalar2<float, Gradient, Hessian> DScalar;
typedef Eigen::Triplet<float> Triplet;

DECLARE_DIFFSCALAR_BASE(); // Initialization of static struct

void Integrator::integrate(Yarn& curYarn, Yarn& nextYarn, Yarn& restYarn, Workspace& ws, Clock& c, Vec3f mousePos) {
  
  float h = c.timestep();
  const size_t NumEqs = curYarn.numCPs() * 3;
  DiffScalarBase::setVariableCount(NUM_VARS);
  
  Eigen::VectorXf Fx    = Eigen::VectorXf::Zero(NumEqs);
  Eigen::VectorXf dqdot = Eigen::VectorXf::Zero(NumEqs);
  Eigen::VectorXf sol   = Eigen::VectorXf::Zero(NumEqs);
  Eigen::SparseMatrix<float> GradFx(NumEqs, NumEqs);
  std::vector<Triplet> triplets;
  
  // WARNING: push_back() on a vector is thread-safe if no allocation is performed. Make sure
  // the correct amount of space is reserved here!
  size_t numTriplets = NUM_VARS*NUM_VARS*curYarn.numIntCPs();
  
  bool converge = false;
  int iterations = 0;
  
  for (int i=0; i<curYarn.numCPs(); i++) {
#ifdef ENABLE_CHECK_NAN
    assert(curYarn.points[i].pos.allFinite());
    assert(curYarn.points[i].vel.allFinite());
#endif
    nextYarn.points[i].pos = curYarn.points[i].pos;
    nextYarn.points[i].vel = curYarn.points[i].vel;
  }
  
  // Perform Newton iteration to solve IMEX equations
  while (!converge) {
    if (iterations != 0) {
      triplets.clear();
      Fx.setZero();
    }
    
    triplets.reserve(numTriplets);
    iterations++;
  
#ifndef PARALLEL
    // NON-PARALLEL:
    
    // Get equations for each point and each velocity
    for (int i=0; i<curYarn.numCPs(); i++) {
      CtrlPoint& curPoint = curYarn.points[i];
      
      if (i > 0 && i<curYarn.numCPs()-1) {
        CtrlPoint& prevPoint = curYarn.points[i-1];
        CtrlPoint& nextPoint = curYarn.points[i+1];
        Segment&   prevSeg   = curYarn.segments[i-1];
        Segment&   nextSeg   = curYarn.segments[i];
        
        typedef DScalar::DVector3 DVector3;
        typedef DScalar::DVector2 DVector2;
        
        // Redefine NUM_VARS if you change these
        DVector3 dPrevPoint(
                      DScalar(0, prevPoint.pos.x() + h*dqdot(3*(i-1))   + h*prevPoint.vel.x()),
                      DScalar(1, prevPoint.pos.y() + h*dqdot(3*(i-1)+1) + h*prevPoint.vel.y()),
                      DScalar(2, prevPoint.pos.z() + h*dqdot(3*(i-1)+2) + h*prevPoint.vel.z()));
        
        DVector3 dCurPoint(
                      DScalar(3, curPoint.pos.x() + h*dqdot(3*i)   + h*curPoint.vel.x()),
                      DScalar(4, curPoint.pos.y() + h*dqdot(3*i+1) + h*curPoint.vel.y()),
                      DScalar(5, curPoint.pos.z() + h*dqdot(3*i+2) + h*curPoint.vel.z()));
        
        DVector3 dNextPoint(
                      DScalar(6, nextPoint.pos.x() + h*dqdot(3*(i+1))   + h*nextPoint.vel.x()),
                      DScalar(7, nextPoint.pos.y() + h*dqdot(3*(i+1)+1) + h*nextPoint.vel.y()),
                      DScalar(8, nextPoint.pos.z() + h*dqdot(3*(i+1)+2) + h*nextPoint.vel.z()));
        
        DVector3 dPrevSeg = dCurPoint - dPrevPoint;
        DVector3 dNextSeg = dNextPoint - dCurPoint;
        assert(dPrevSeg.norm() != 0 && dNextSeg.norm() != 0 && "Edge length is 0");
        DVector3 dPrevSegN = dPrevSeg.normalized();
        DVector3 dNextSegN = dNextSeg.normalized();
        DScalar dotProd = dPrevSegN.dot(dNextSegN);
        assert(dotProd != -1 && "Segments are pointing in exactly opposite directions");
        
        DVector3 curveBinorm = (DScalar(2)*dPrevSegN.cross(dNextSegN))/(1+dotProd);
        
        DVector3 d1prev, d1next, d2prev, d2next;
        
        d1prev << DScalar(prevSeg.m1().x()),
                  DScalar(prevSeg.m1().y()),
                  DScalar(prevSeg.m1().z());
        
        d2prev << DScalar(prevSeg.m2().x()),
                  DScalar(prevSeg.m2().y()),
                  DScalar(prevSeg.m2().z());
        
        d1next << DScalar(nextSeg.m1().x()),
                  DScalar(nextSeg.m1().y()),
                  DScalar(nextSeg.m1().z());
        
        d2next << DScalar(nextSeg.m2().x()),
                  DScalar(nextSeg.m2().y()),
                  DScalar(nextSeg.m2().z());
        
        DVector2 matCurvePrev(curveBinorm.dot(d2prev), -curveBinorm.dot(d1prev));
        DVector2 matCurveNext(curveBinorm.dot(d2next), -curveBinorm.dot(d1next));
        DVector2 matCurve = DScalar(0.5)*(matCurvePrev + matCurveNext);
        DVector2 restMatCurvePrev, restMatCurveNext;
        
        restMatCurvePrev << DScalar(ws.restMatCurvature[i][i-1].x()),
                            DScalar(ws.restMatCurvature[i][i-1].y());
        
        restMatCurveNext << DScalar(ws.restMatCurvature[i][i].x()),
                            DScalar(ws.restMatCurvature[i][i].y());
        DVector2 restMatCurve = DScalar(0.5)*(restMatCurvePrev + restMatCurveNext);
        
        // TODO: bending matrix may not be I
        
        DScalar voronoiCell = 0.5*(dPrevSeg.norm()+dNextSeg.norm());
        
        DScalar bendEnergy = 0.5*(1/voronoiCell)*
        (matCurve - restMatCurve).dot(matCurve - restMatCurve);
        
        Gradient grad = bendEnergy.getGradient();
        Hessian hess = bendEnergy.getHessian();
        
        for (int j=0; j<NUM_VARS; j++) {
          // TODO: mass matrix may not be I
          Fx(3*(i-1)+j) -= h*grad(j);
          for (int k=0; k<NUM_VARS; k++) {
            float val = -h*h*hess(j,k);
            CHECK_NAN(val);
            if (val != 0) {
              triplets.push_back(Triplet(3*(i-1)+j, 3*(i-1)+k, val));
            }
          }
        }
      }
      
      if (i > 1) {
        Fx(3*i+1) += 10*h; // gravity hack
      }
      
      if (i+1 == curYarn.numCPs()) {
        Fx(3*(i-1)+1) -= h * (mousePos.y() - curYarn.points[i].pos.y());
        Fx(3*(i-1)+2) -= h * (mousePos.z() - curYarn.points[i].pos.z());
      }
    }
#else
    // PARALLEL

    // Check for data race issues
    assert(triplets.size() == numTriplets);
#endif // ifndef PARALLEL
    
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
      assert(false);
    } else {
      std::cout << "CG iters: " << cg.iterations() << "\n";
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
  for (int i=0; i<curYarn.numCPs(); i++) {
    Vec3f curdqdot = dqdot.block<3, 1>(3*i, 0);
    nextYarn.points[i].pos += h*nextYarn.points[i].vel + h*curdqdot;
    nextYarn.points[i].vel += curdqdot;
  }
}