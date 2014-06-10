//
//  ConstraintIntegrator.h
//  Visualizer
//
//  Created by eschweickart on 6/6/14.
//
//

#ifndef __Visualizer__ConstraintIntegrator__
#define __Visualizer__ConstraintIntegrator__

#include "Integrator.h"
#include "Constraint.h"

class ConstraintIntegrator : public Integrator {
private:
  std::vector<YarnConstraint*> constraints;
  
public:
  
  ConstraintIntegrator(Yarn&, std::vector<YarnEnergy*>, std::vector<YarnConstraint*>);
  bool integrate(Clock& c);
};


#endif /* defined(__Visualizer__ConstraintIntegrator__) */
