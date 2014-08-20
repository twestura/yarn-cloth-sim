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

  /// Target timestep of the simulation in seconds.
  const real INITIAL_TIMESTEP = 1.0/60.0;

  /// The default radius of the rod in meters.
  const real radius = 0.00635;
  // Mascarenhas test: 0.00955
  // 1/2 inch diameter = 1/4 inch radius: 0.00635
  // 3/4 inch diameter = 3/8 inch radius: 0.009525
  
  /// The default Young's modulus of the rod in N/m^2.
  const real youngsModulus = 6.89e10;
  // Steel (1018): 2.05e11
  // Aluminum (6061) : 6.89e10
  
  /// The default Poisson's ratio of the rod (dimensionless).
  const real poissonRatio = 0.33;
  // Steel (1018): 0.29
  // Aluminum (6061) : 0.33
  
  /// The default shear modulus of the rod in N/m^2.
  const real shearModulus = 2.6e10; //  youngsModulus / (2.0 * (1.0 + poissonRatio));
  // Steel (1018): 8.0e10
  // Aluminum (6061): 2.6e10
  
  /// The default density of the rod in kg/m^3.
  const real rhoRod = 2700;
  // Steel (1018): 7870
  // Aluminum (6061): 2700
  
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
  
  /// The path for outputting media files.
  const std::string ResultPath = "../../result/";
}

#endif
