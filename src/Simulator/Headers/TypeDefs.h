//
//  Simulator_Prefix.pch
//  Visualizer
//
//  Created by eschweickart on 3/5/14.
//
//

#ifndef Visualizer_Simulator_Prefix_pch
#define Visualizer_Simulator_Prefix_pch

// Include any system framework and library headers here that should be included in all compilation units.
// You will also need to set the Prefix Header build setting of one or more of your targets to reference this file.

#include "Macros.h"

#ifdef __OBJC__
#import <Cocoa/Cocoa.h>
#endif

#if defined( __cplusplus )
#include "cinder/Cinder.h"

#include "cinder/app/AppBasic.h"
#include "cinder/gl/gl.h"

#include "cinder/CinderMath.h"
#include "cinder/Matrix.h"
#include "cinder/Vector.h"
#include "cinder/Quaternion.h"

#include "Eigen/dense"
#include "Eigen/sparse"

// Global typedefs
typedef double real;

typedef Eigen::Matrix<real, 2, 2> Mat2e;
typedef Eigen::Matrix<real, 3, 3> Mat3e;

typedef Eigen::Matrix<real, 2, 1> Vec2e;
typedef Eigen::Matrix<real, 3, 1> Vec3e;
typedef Eigen::Matrix<real, 4, 1> Vec4e;
typedef Eigen::Matrix<real, Eigen::Dynamic, 1> VecXe;

typedef ci::Vec2<real> Vec2c;
typedef ci::Vec3<real> Vec3c;

typedef Eigen::Triplet<real> Triplet;
#endif

#endif
