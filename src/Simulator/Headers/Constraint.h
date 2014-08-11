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
#include "Rod.h"

/// An abstract class representing a position-based constraint on a rod.
class RodConstraint {
protected:
  Rod& r;
public:
  RodConstraint(Rod& r) : r(r) { }
  virtual bool eval(VecXe&, real)=0;
  virtual ~RodConstraint() { }
};

class Length : public RodConstraint {
public:
  Length(Rod& r) : RodConstraint(r) { }
  bool eval(VecXe&, real);
};

#endif /* defined(__Visualizer__Constraint__) */
