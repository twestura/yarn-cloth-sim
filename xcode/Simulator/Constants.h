//
//  Constants.h
//  Visualizer
//
//  Created by eschweickart on 2/13/14.
//
//

#ifndef Visualizer_Constants_h
#define Visualizer_Constants_h

/// Timestep of the simulation in seconds.
const float h = (1.0/24000.0);

/// Determines centerline plasticity due to bending forces on the yarn.
const float pPlastic = 0.01;
/// Determines the maximum amount of plasticity the centerline of the yarn exhibits.
const float pPlasticMax = 2.5;


#endif
