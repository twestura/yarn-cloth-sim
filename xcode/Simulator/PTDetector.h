//
//  PTDetector.h
//  Visualizer
//
//  Created by eschweickart on 4/25/14.
//
//

#ifndef __Visualizer__PTDetector__
#define __Visualizer__PTDetector__

#include <iostream>
#include "Eigen/Sparse"

/*
class PTDetector {
  Eigen::SparseMatrix<float> mat;
  
public:
  
  PTDetector(size_t size) : mat(Eigen::SparseMatrix<float>(size, size)) {}
  
  void contact(size_t i, size_t j, float priority) {
    if (i > j) { // store upper half
      size_t temp = i;
      i = j;
      j = temp;
    }
    assert(mat.coeff(i, j) == 0 && "You forgot to reset!");
    mat.insert(i, j) = priority;
  }
  
  void iter() {
    
  }
  
  float topPriority() {
    return mat.toDense().maxCoeff(); // not great...
  }
  
  void reset() {
    mat.setZero();
  }
  
};
*/
 
class PTDetector {
  int i;
  int j;
  
  float s1dot;
  float s2dot;
public:
  
  PTDetector(int i, int j, float s1dot, float s2dot) :
    i(i), j(j), s1dot(s1dot), s2dot(s2dot) { }
  
  std::pair<int, int> id() { return std::pair<int, int>(i, j); }
  
  bool pass(float s1dotnew, float s2dotnew) {
    bool ret = ((s1dotnew > 0 && s1dot < 0) || (s1dotnew < 0 && s1dot > 0)) &&
               ((s2dotnew > 0 && s2dot < 0) || (s2dotnew < 0 && s2dot > 0));
    s1dot = s1dotnew;
    s2dot = s2dotnew;
    return ret;
  }
  
};

#endif /* defined(__Visualizer__PTDetector__) */
