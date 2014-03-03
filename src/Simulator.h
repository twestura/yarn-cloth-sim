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
  /// The yarn at time 0.
  Yarn restYarn;
  /// The yarn at time t.
  Yarn* curYarn;
  /// The yarn at time t+1.
  Yarn* nextYarn;
  
  /// Space for calculations to avoid repetitive, expensive allocation.
  Workspace ws;
  
  /// Compute F(t+1).
  void computeForces();
public:
  /// Entry point for simulator; called after the app has finished initializing
  void run();
};

#endif /* defined(__Visualizer__Simulator__) */
