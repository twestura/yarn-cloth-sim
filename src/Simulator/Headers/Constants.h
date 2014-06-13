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
  const float INITIAL_TIMESTEP = 1.0f/60.0f;

  /// Determines centerline plasticity due to bending forces on the yarn.
  const float pPlastic = 0.01f;
  /// Determines the maximum amount of plasticity the centerline of the yarn exhibits.
  const float pPlasticMax = 2.5f;

  /// The default radius of the yarn.
  const float radius = 0.15f;
  
  // The default shear modulus of the yarn.
  const float shearModulus = 8e2f;

  // the default Young's modulus of the yarn.
  const float youngsModulus = 2e7f;
  
  /// Pi. You know the one.
  const float pi = 3.1415926535f;
  
  /// The spline basis for the Catmull-Rom spline, separated into 4 vectors for convenience.
  const Vec4f basis[4] = { Vec4f(-0.5f,  1.0f, -0.5f, 0.0f),
                           Vec4f( 1.5f, -2.5f,  0.0f, 1.0f),
                           Vec4f(-1.5f,  2.0f,  0.5f, 0.0f),
                           Vec4f( 0.5f, -0.5f,  0.0f, 0.0f) };
  
  /// The number of quadrature points per spline segment.
  const int numQuadPoints = 6;
  
  /// The path for outputting media files
  const std::string ResultPath = "../../result/";

  
}

#endif
