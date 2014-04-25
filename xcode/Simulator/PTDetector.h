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
