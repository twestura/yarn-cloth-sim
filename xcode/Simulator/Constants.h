//
//  Constants.h
//  Visualizer
//
//  Created by eschweickart on 2/13/14.
//
//

#ifndef Visualizer_Constants_h
#define Visualizer_Constants_h

namespace constants {

  /// Timestep of the simulation in seconds.
  const float INITIAL_TIMESTEP = (1.0/30.0);

  /// Determines centerline plasticity due to bending forces on the yarn.
  const float pPlastic = 0.01;
  /// Determines the maximum amount of plasticity the centerline of the yarn exhibits.
  const float pPlasticMax = 2.5;

  /// Determines the radius of the yarn.
  const float radius = 0.05;

  /// Pi. You know the one.
  const float pi = 3.1415926535;
  
}

#endif
