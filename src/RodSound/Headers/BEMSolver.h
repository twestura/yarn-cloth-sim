//
//  BEMSolver.h
//  Visualizer
//
//  Created by eschweickart on 7/23/14.
//
//

#ifndef __Visualizer__BEMSolver__
#define __Visualizer__BEMSolver__

#include "Constants.h"
#include "Util.h"

static inline Vec2e circle_pos(real t, real r) {
  t *= 2.0 * constants::pi;
  return Vec2e(r * cos(t), r * sin(t));
}

static inline Vec2e circle_normal(real t) {
  t *= 2.0 * constants::pi;
  return Vec2e(cos(t), sin(t));
}


static Mat2e solveBEM(real r) {
  /// This is the number of samples to use on the cross section we're solving for.
  uint32 n_N = 100;
  
  /// This is a real factor that determines the radius (and sampling rate) of the outer circle.
  /// At infinity, this is most accurate: The Krichoff vector has a 1/alpha falloff, which we
  /// approximate as 0 at the outer boundary.
  uint32 alpha = 40;
  
  uint32 n_D = n_N * alpha;
  
  uint32 n = n_N + n_D;
  
  Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> C(n, n);
  Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> A_N(n, n_N);
  Eigen::Matrix<real, Eigen::Dynamic, 2> n_i(n_N, 2);
  Eigen::Matrix<real, Eigen::Dynamic, 2> u(n, 2);
  
  for (int i=0; i<n; i++) {
    if (i < n_N) {
      n_i.row(i) = circle_normal(((real) i) / n_N).transpose();
    }
    
    Vec2e xi, xi1;
    if (i < n_N) {
      xi = circle_pos(((real) i) / n_N, r);
    } else {
      xi = circle_pos(((real) i-n_N) / n_D, r*alpha);
    }
    if (i < n_N-1) {
      xi1 = circle_pos(((real) i+1) / n_N, r);
    } else if (i == n_N-1) {
      xi1 = circle_pos(0, r);
    } else if (i < n-1) {
      xi1 = circle_pos(((real) i+1-n_N) / n_D, r*alpha);
    } else {
      xi1 = circle_pos(0, r*alpha);
    }
    Vec2e yi = (xi - xi1) / 2.0;
    
    for (int j=0; j<n; j++) {
      if (i == j) {
        real h = (xi1 - xi).norm();
        if (i < n_N) {
          C(i, i) = 0.5;
          A_N(i, i) = -0.5 / constants::pi * h * (log(h / 2.0) - 2.0);
        } else {
          C(i, i) = 0.5 / constants::pi * h * (log(h / 2.0) - 2.0);
        }
      } else {
        Vec2e xj, xj1;
        if (j < n_N) {
          xj = circle_pos(((real) j) / n_N, r);
        } else {
          xj = circle_pos(((real) j-n_N) / n_D, r*alpha);
        }
        if (j < n_N-1) {
          xj1 = circle_pos(((real) j+1) / n_N, r);
        } else if (j == n_N-1) {
          xj1 = circle_pos(0, r);
        } else if (j < n-1) {
          xj1 = circle_pos(((real) j+1-n_N) / n_D, r*alpha);
        } else {
          xj1 = circle_pos(0, r*alpha);
        }
        Vec2e x = xj1 - xj;
        
        real a = ((-x.y())*yi.x() + (x.x())*yi.y() + xj.x()*xj1.y() - xj1.x()*xj.y()) / x.norm();
        
        if (fabs(a) < 1.0e-7) { // Colinear
          real upperLimit = (xj1 - yi).norm();
          real lowerLimit = (xj - yi).norm();
          if (j < n_N) {
            C(i, j) = 0;
            A_N(i, j) = -0.5 / constants::pi * ((upperLimit * (log(upperLimit) - 1.0)) -
                                                (lowerLimit * (log(lowerLimit) - 1.0)));
          } else {
            C(i, j) = 0.5 / constants::pi * ((upperLimit * (log(upperLimit) - 1.0)) -
                                             (lowerLimit * (log(lowerLimit) - 1.0)));
          }
        } else { // Not Colinear
          Vec2e xNorm = x.normalized();
          Vec2e z = xj + (yi - xj).dot(xNorm) * xNorm;
          real cosTheta1 = (xj - yi).dot(z - yi) / (xj - yi).norm() / (z - yi).norm();
          real cosTheta2 = (xj1 - yi).dot(z - yi) / (xj1 - yi).norm() / (z - yi).norm();
          
          real theta1, theta2;
          if (cosTheta1 >= 1.0) { // cosTheta will never be <= 0
            theta1 = 0.0;
          } else {
            theta1 = acos(cosTheta1);
          }
          if (cosTheta2 >= 1.0) {
            theta2 = 0.0;
          } else {
            theta2 = acos(cosTheta2);
          }
          
          // Correct the signs
          if (a < 0.0) {
            theta1 = -theta1;
            theta2 = -theta2;
          }
          if ((z - xj).dot(x) > 0.0) theta1 = -theta1;
          if ((z - xj1).dot(x) > 0.0) theta2 = -theta2;
          
          a = fabs(a);
          
          if (j < n_N) {
            C(i, j) = -0.5 / constants::pi * (theta2 - theta1);
            A_N(i, j) = -0.5 / constants::pi * ((tan(theta2)*(log(a / cosTheta2) - 1.0) + theta2) -
                                                (tan(theta1)*(log(a / cosTheta1) - 1.0) + theta1));
          } else {
            C(i, j) = 0.5 / constants::pi * ((tan(theta2)*(log(a / cosTheta2) - 1.0) + theta2) -
                                             (tan(theta1)*(log(a / cosTheta1) - 1.0) + theta1));
          }
        }
      }
    }
  }
  
  Eigen::Matrix<real, Eigen::Dynamic, 2> f;
  f.noalias() = A_N * n_i;
  Eigen::Matrix<real, Eigen::Dynamic, 2> sol = C.colPivHouseholderQr().solve(f).topRows(n_N);
  
  Mat2e ret;
  ret.noalias() = n_i.transpose() * sol;
  return ret;
}



#endif /* defined(__Visualizer__BEMSolver__) */
