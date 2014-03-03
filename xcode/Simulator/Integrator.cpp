//
//  Integrator.cpp
//  Visualizer
//
//  Created by eschweickart on 2/24/14.
//
//

#include "Integrator.h"
#include <stdexcept>

#define NUM_VARS 6
const float ConvergenceThreshold = 0.001;

typedef Eigen::Matrix<float, NUM_VARS, 1> Gradient;
typedef Eigen::Matrix<float, NUM_VARS, NUM_VARS> Hessian;
typedef DScalar2<float, Gradient, Hessian> DScalar;
typedef Eigen::Triplet<float> Triplet;

DECLARE_DIFFSCALAR_BASE(); // Initialization of static struct

void Integrator::integrate(Yarn& curYarn, Yarn& nextYarn, Yarn& restYarn, Workspace& ws, Clock& c) {
  
  DiffScalarBase::setVariableCount(NUM_VARS);
  
  // TODO: remove hard-coded numbers
  Eigen::VectorXf Fx(curYarn.numCPs()*6);
  Eigen::VectorXf dq(curYarn.numCPs()*3);
  Eigen::VectorXf dqdot(curYarn.numCPs()*3);
  Eigen::SparseMatrix<float> GradFx(curYarn.numCPs()*6, curYarn.numCPs()*6);
  std::vector<Triplet> triplets;
  
  // WARNING: push_back() on a vector is thread-safe if no allocation is performed. Make sure
  // the correct amount of space is reserved here!
  triplets.reserve(NUM_VARS*NUM_VARS*curYarn.numCPs());
  
  bool converge = false;
  
  // TODO: set nextYarn appropriately
  
  // Perform Newton iteration to solve IMEX equations
  while (!converge) {
  
#ifndef PARALLEL
    // NON-PARALLEL:
    
    // Get equations for each point and each velocity
    for (int i=0; i<curYarn.numCPs(); i++) {
      
      if (i > 0 && i<curYarn.numCPs()-1) {
        
        typedef Eigen::Matrix<DScalar, 3, 1> DVector3;
        typedef Eigen::Matrix<DScalar, 2, 1> DVector2;
        
        // Redefine NUM_VARS if you change these
        DVector3 pos, vel;
        pos << DScalar(0, curYarn.points[i].pos.x() + dq(NUM_VARS*i)),
               DScalar(1, curYarn.points[i].pos.y() + dq(NUM_VARS*i+1)),
               DScalar(2, curYarn.points[i].pos.z() + dq(NUM_VARS*i+2));
        
        vel << DScalar(3, curYarn.points[i].vel.x() + dqdot(NUM_VARS*i)),
               DScalar(4, curYarn.points[i].vel.y() + dqdot(NUM_VARS*i+1)),
               DScalar(5, curYarn.points[i].vel.z() + dqdot(NUM_VARS*i+2));

        DVector3 prevSeg, nextSeg;
        prevSeg << pos.x() - curYarn.points[i-1].pos.x(),
                   pos.y() - curYarn.points[i-1].pos.y(),
                   pos.z() - curYarn.points[i-1].pos.z();
        
        nextSeg << curYarn.points[i+1].pos.x() - pos.x(),
                   curYarn.points[i+1].pos.y() - pos.y(),
                   curYarn.points[i+1].pos.z() - pos.z();
        
        DVector3 curveBinorm = (DScalar(2)*prevSeg.normalized().cross(nextSeg.normalized())) /
           (1+prevSeg.normalized().dot(nextSeg.normalized()));
        
        DVector3 d1prev, d1next, d2prev, d2next;
        
        d1prev << DScalar(curYarn.segments[i-1].m1().x()),
                  DScalar(curYarn.segments[i-1].m1().y()),
                  DScalar(curYarn.segments[i-1].m1().z());
        
        d2prev << DScalar(curYarn.segments[i-1].m2().x()),
                  DScalar(curYarn.segments[i-1].m2().y()),
                  DScalar(curYarn.segments[i-1].m2().z());
        
        d1next << DScalar(curYarn.segments[i].m1().x()),
                  DScalar(curYarn.segments[i].m1().y()),
                  DScalar(curYarn.segments[i].m1().z());
        
        d2next << DScalar(curYarn.segments[i].m2().x()),
                  DScalar(curYarn.segments[i].m2().y()),
                  DScalar(curYarn.segments[i].m2().z());
        
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
        
        DScalar voronoiCell = DScalar(0.5)*(prevSeg.norm()+nextSeg.norm());
        
        DScalar bendEnergy = DScalar(0.5)*(DScalar(1)/voronoiCell)*
        (matCurve - restMatCurve).dot(matCurve - restMatCurve);
        
        Gradient grad = bendEnergy.getGradient();
        Hessian hess = bendEnergy.getHessian();
        
        // TODO: remove hard-coded 3
        for (int j=0; j<3; j++) {
          // TODO: mass matrix may not be I
          Fx(3*i+j) = dqdot(3*i+j) -c.timestep()*grad(j);
          Fx(3*curYarn.numCPs()+3*i+j) = dq(3*i+j) - c.timestep()*dqdot(3*i+j)
           - c.timestep()*curYarn.points[i].vel(j);
          for (int k=0; k<3; k++) {
            triplets.push_back(Triplet(3*i+j, 3*i+k, -c.timestep()*grad(j, k)));
            if (j==k) {
              triplets.push_back(Triplet(3*i+j, 3*curYarn.numCPs()*3*i+k, 1));
              triplets.push_back(Triplet(3*curYarn.numCPs()*3*i+j, 3*i+k, 1));
              triplets.push_back(Triplet(3*curYarn.numCPs()*3*i+j, 3*curYarn.numCPs()*3*i+k, c.timestep()));
            }
          }
        }
      }
    }
#else
    // PARALLEL

#endif // ifndef PARALLEL
    
    // Check for data race issues
    assert(triplets.size() == NUM_VARS*NUM_VARS*curYarn.numCPs());
    // Solve equations for updates to changes in position and velocity using Conjugate Gradient
    GradFx.setFromTriplets(triplets.begin(), triplets.end());
    Eigen::ConjugateGradient<Eigen::SparseMatrix<float>> cg;
    cg.compute(GradFx);
    Eigen::VectorXf sol(6*curYarn.numCPs());
    // TODO: solve with guess?
    sol = cg.solve(Fx);
    
    if (cg.info() == Eigen::NoConvergence) {
      throw std::runtime_error("Error: Conjugate Gradient solver in integrator did not converge.");
    }
    
    // Update changes to position and velocity
    dq += sol.block(0, 0, 3*curYarn.numCPs(), 1);
    dqdot += sol.block(0, 1, 3*curYarn.numCPs(), 1);
    
    if (sol.maxCoeff() < ConvergenceThreshold) converge = true;
    
  }
}