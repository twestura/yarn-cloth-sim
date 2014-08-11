//
//  ExIntegrator.h
//  Visualizer
//
//  Created by eschweickart on 6/10/14.
//
//

#ifndef __Visualizer__ExIntegrator__
#define __Visualizer__ExIntegrator__

#include "Integrator.h"
#include "Constraint.h"

class ExIntegrator : public Integrator {
  Eigen::SparseMatrix<real> damping;
  Eigen::SparseMatrix<real> stiffness;
  real alpha1;
  real alpha2;
  void setDamping();
  
  std::vector<RodConstraint*> constraints;
  Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> modes;
public:
  ExIntegrator(Rod&, std::vector<RodEnergy*>&, std::vector<RodConstraint*>* = nullptr);
  bool integrate(Clock&);
  void draw();
};

#endif /* defined(__Visualizer__ExIntegrator__) */
