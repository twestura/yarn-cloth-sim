//
//  Constants.h
//  Visualizer
//
//  Created by eschweickart on 2/13/14.
//
//

#ifndef Visualizer_Constants_h
#define Visualizer_Constants_h

typedef Eigen::Vector4f Vec4f;

namespace constants {

  /// Timestep of the simulation in seconds.
  const float INITIAL_TIMESTEP = (1.0/60.0);

  /// Determines centerline plasticity due to bending forces on the yarn.
  const float pPlastic = 0.01;
  /// Determines the maximum amount of plasticity the centerline of the yarn exhibits.
  const float pPlasticMax = 2.5;

  /// Determines the radius of the yarn.
  const float radius = 0.15;

  /// Pi. You know the one.
  const float pi = 3.1415926535;
  
  const Vec4f basis[4] = { Vec4f(-0.5,  1,  -0.5, 0),
                           Vec4f( 1.5, -2.5, 0,   1),
                           Vec4f(-1.5,  2,   0.5, 0),
                           Vec4f( 0.5, -0.5, 0,   0) };
  
}

#endif
