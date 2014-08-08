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

class YarnConstraint {
protected:
  Yarn& y;
public:
  YarnConstraint(Yarn& y) : y(y) { }
  virtual bool eval(VecXe&, real)=0;
  virtual ~YarnConstraint() { }
};

class Length : public YarnConstraint {
public:
  Length(Yarn& y) : YarnConstraint(y) { }
  bool eval(VecXe&, real);
};

#endif /* defined(__Visualizer__Constraint__) */
