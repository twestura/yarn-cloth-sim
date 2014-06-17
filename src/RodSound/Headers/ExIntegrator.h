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

class ExIntegrator : public Integrator {
  Eigen::SparseMatrix<float> damping;
public:
  ExIntegrator(Yarn&, std::vector<YarnEnergy*>&);
  bool integrate(Clock&);
};

#endif /* defined(__Visualizer__ExIntegrator__) */
