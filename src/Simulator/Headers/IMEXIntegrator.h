//
//  Integrator.h
//  Visualizer
//
//  Gathers energies (either implicit or explicit) and integrates the given
//  object by a given amount using an IMEX scheme.
//
//  Created by eschweickart on 2/24/14.
//
//

#ifndef Visualizer_IMEXIntegrator_h
#define Visualizer_IMEXIntegrator_h

#include "Integrator.h"
#include "Util.h"
#include "autodiff.h"
#include "Eigen/Sparse"


class IMEXIntegrator : public Integrator {
private:
  Eigen::SparseMatrix<float> hessBase;
public:
  IMEXIntegrator(std::vector<YarnEnergy*>& energies, Yarn& y);
  
  bool integrate(Clock& c);
  bool setRotations() const;  
};


#endif
