//
//  ConstraintIntegrator.h
//  Visualizer
//
//  Created by eschweickart on 6/6/14.
//
//

#ifndef __Visualizer__ConstraintIntegrator__
#define __Visualizer__ConstraintIntegrator__

#include <iostream>
#include "Energy.h"
#include "Constraint.h"

class ConstraintIntegrator {
  Yarn& y;
  std::vector<YarnEnergy*> energies;
  std::vector<YarnConstraint*> constraints;
  
public:
  
  ConstraintIntegrator(Yarn&, std::vector<YarnEnergy*>, std::vector<YarnConstraint*>);
  bool integrate(Clock& c);
};


#endif /* defined(__Visualizer__ConstraintIntegrator__) */
