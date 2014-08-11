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

/// An integrator that combines explicitly evaluated forces and position-based constraints.
class ConstraintIntegrator : public Integrator {
private:
  /// The position-based constraints imposed on the rod.
  std::vector<RodConstraint*> constraints;
  
public:
  
  ConstraintIntegrator(Rod&, std::vector<RodEnergy*>&, std::vector<RodConstraint*>&);
  bool integrate(Clock& c);
};


#endif /* defined(__Visualizer__ConstraintIntegrator__) */
