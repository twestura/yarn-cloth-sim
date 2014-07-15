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
  const real INITIAL_TIMESTEP = 1.0/60.0;

  // TODO: Use these.
  /// Determines centerline plasticity due to bending forces on the yarn.
  const real pPlastic = 0.01;
  /// Determines the maximum amount of plasticity the centerline of the yarn exhibits.
  const real pPlasticMax = 2.5;

  /// The default radius of the rod.
  const real radius = 0.15;
  
  /// The default Young's modulus of the rod.
  constexpr real youngsModulus = 2.0e9;
  
  /// The default Poisson's ratio of the rod.
  constexpr real poissonRatio = 0.3;
  
  /// The default shear modulus of the rod.
  constexpr real shearModulus = youngsModulus / (2.0 * (1.0 + poissonRatio));
  
  /// Pi. You know the one.
  const real pi = 3.1415926535;
  
  /// The speed of sound in air, measured in m/s.
  const real cAir = 340.0;
  
  /// The density of air, measured in kg/m^3.
  const real rhoAir = 1.23;
  
  /// The spline basis for the Catmull-Rom spline, separated into 4 vectors for convenience.
  const Vec4e basis[4] = { Vec4e(-0.5,  1.0, -0.5, 0.0),
                           Vec4e( 1.5, -2.5,  0.0, 1.0),
                           Vec4e(-1.5,  2.0,  0.5, 0.0),
                           Vec4e( 0.5, -0.5,  0.0, 0.0) };
  
  /// The number of quadrature points per spline segment.
  const int numQuadPoints = 6;
  
  /// The path for outputting media files
  const std::string ResultPath = "../../result/";
}

#endif
