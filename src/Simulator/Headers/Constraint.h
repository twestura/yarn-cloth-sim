//
//  Constraint.h
//  Visualizer
//
//  Created by eschweickart on 6/6/14.
//
//

#ifndef __Visualizer__Constraint__
#define __Visualizer__Constraint__

#include <iostream>
#include "Yarn.h"

typedef Eigen::VectorXf VecXf;

class YarnConstraint {
protected:
  Yarn& y;
public:
  YarnConstraint(Yarn& y) : y(y) { }
  virtual bool eval(VecXf&, float)=0;
};

class Length : public YarnConstraint {
public:
  Length(Yarn& y) : YarnConstraint(y) { }
  bool eval(VecXf&, float);
};

#endif /* defined(__Visualizer__Constraint__) */
