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

#ifndef Visualizer_Integrator_h
#define Visualizer_Integrator_h

#include <stdexcept>
#include "Yarn.h"
#include "Util.h"
#include "autodiff.h"
#include "Eigen/Sparse"
#include "Clock.h"
#include "Energy.h"


class Integrator {
private:
  std::vector<YarnEnergy*>& energies;
public:
  Integrator(std::vector<YarnEnergy*>& energies) : energies(energies) {}
  
  void integrate(Yarn& y, Clock& c);
  
};


#endif
