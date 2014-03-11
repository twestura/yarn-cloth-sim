//
//  Simulator.h
//  Visualizer
//
//  Created by eschweickart on 2/10/14.
//
//

#ifndef __Visualizer__Simulator__
#define __Visualizer__Simulator__

#include <iostream>
#include <boost/thread/thread.hpp>
#include "Yarn.h"
#include "Integrator.h"
#include "Util.h"

class Simulator
{
  Yarn* y;
  Integrator* integrator;

public:
  /// Entry point for simulator; called after the app has finished initializing
  void run();
  
  void computeForces();
};

#endif /* defined(__Visualizer__Simulator__) */
