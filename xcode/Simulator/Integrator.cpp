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

void Integrator::integrate(Yarn& curYarn, Yarn& nextYarn, Yarn& restYarn, Workspace& ws, Clock& c) {
  float h = c.timestep();
  const size_t NumEqs = curYarn.numCPs() * 3;
  DiffScalarBase::setVariableCount(NUM_VARS);
  
  Eigen::VectorXf Fx(NumEqs);
  Eigen::VectorXf dqdot(NumEqs);
  Eigen::VectorXf sol(NumEqs);
  Eigen::SparseMatrix<float> GradFx(NumEqs, NumEqs);
  std::vector<Triplet> triplets;
  
  // WARNING: push_back() on a vector is thread-safe if no allocation is performed. Make sure
  // the correct amount of space is reserved here!
  size_t numTriplets = NUM_VARS*NUM_VARS*curYarn.numIntCPs();
  
  bool converge = false;
  int iterations = 0;
  
  for (int i=0; i<curYarn.numCPs(); i++) {
    nextYarn.points[i].pos = curYarn.points[i].pos;
    nextYarn.points[i].vel = curYarn.points[i].vel;
  }
  
  // Perform Newton iteration to solve IMEX equations
  while (!converge) {
    if (iterations != 0) {
      triplets.clear();
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
        
        typedef Eigen::Matrix<DScalar, 3, 1> DVector3;
        typedef Eigen::Matrix<DScalar, 2, 1> DVector2;
        
        // Redefine NUM_VARS if you change these
        DVector3 dPrevPoint, dCurPoint, dNextPoint;
        dPrevPoint << DScalar(0, prevPoint.pos.x() + h*dqdot(3*(i-1))   + h*prevPoint.vel.x()),
                      DScalar(1, prevPoint.pos.y() + h*dqdot(3*(i-1)+1) + h*prevPoint.vel.y()),
                      DScalar(2, prevPoint.pos.z() + h*dqdot(3*(i-1)+2) + h*prevPoint.vel.z());
        
        dCurPoint << DScalar(3, curPoint.pos.x() + h*dqdot(3*i)   + h*curPoint.vel.x()),
                     DScalar(4, curPoint.pos.y() + h*dqdot(3*i+1) + h*curPoint.vel.y()),
                     DScalar(5, curPoint.pos.z() + h*dqdot(3*i+2) + h*curPoint.vel.z());
        
        dNextPoint << DScalar(6, nextPoint.pos.x() + h*dqdot(3*(i+1))   + h*nextPoint.vel.x()),
                      DScalar(7, nextPoint.pos.y() + h*dqdot(3*(i+1)+1) + h*nextPoint.vel.y()),
                      DScalar(8, nextPoint.pos.z() + h*dqdot(3*(i+1)+2) + h*nextPoint.vel.z());
        
        DVector3 dPrevSeg = dCurPoint - dPrevPoint;
        DVector3 dNextSeg = dNextPoint - dCurPoint;
        
        DVector3 curveBinorm = (DScalar(2)*dPrevSeg.normalized().cross(dNextSeg.normalized())) /
           (1+dPrevSeg.normalized().dot(dNextSeg.normalized()));
        
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
          Fx(3*(i-1)+j) += dqdot(3*(i-1)+j) -h*grad(j); // Add external forces here
          for (int k=0; k<NUM_VARS; k++) {
            triplets.push_back(Triplet(3*(i-1)+j, 3*(i-1)+k, (j==k ? 1 : 0)-h*h*hess(j, k)));
          }
        }
      }
    }
#else
    // PARALLEL

    // Check for data race issues
    assert(triplets.size() == numTriplets);
#endif // ifndef PARALLEL
    
    // Solve equations for updates to changes in position and velocity using Conjugate Gradient
    GradFx.setFromTriplets(triplets.begin(), triplets.end()); // sums up duplicates automagically
    Eigen::ConjugateGradient<Eigen::SparseMatrix<float>> cg;
    cg.compute(GradFx);
    sol = cg.solveWithGuess(Fx, sol); // Aliasing problems???
    
    if (cg.info() == Eigen::NoConvergence) {
      throw std::runtime_error("Error: Conjugate Gradient solver in integrator did not converge.");
    }
    
    // Update changes to position and velocity
    dqdot -= sol;
    for (int i=0; i<curYarn.numCPs(); i++) {
      Vec3f curdqdot = dqdot.block<3, 1>(3*i, 1);
      nextYarn.points[i].pos += h*nextYarn.points[i].vel + h*curdqdot;
      nextYarn.points[i].vel += curdqdot;
    }
    
    
    if (sol.maxCoeff() < ConvergenceThreshold) {
      converge = true;
    } else if (iterations > 50) {
      std::cerr << "Too many newton iterations; we likely diverged.";
      break;
    }
  }
}