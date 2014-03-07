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

#include "Yarn.h"
#include "Util.h"
#include "autodiff.h"
#include "Eigen/Sparse"
#include "Clock.h"


class Integrator
{
  public:
  static void integrate(Yarn& curYarn, Yarn& nextYarn, Yarn& restYarn, Workspace& ws, Clock& c, Vec3f mousePos);
  
};


#endif
