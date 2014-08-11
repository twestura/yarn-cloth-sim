//
//  Integrator.h
//  Visualizer
//
//  Created by eschweickart on 6/10/14.
//
//

#ifndef Visualizer_Integrator_h
#define Visualizer_Integrator_h

#include "Rod.h"
#include "Clock.h"
#include "Energy.h"

/// An abstract class that moves a rod forward in time.
class Integrator {
protected:
  /// A reference to the rod affected.
  Rod& r;
  /// The energies affecting the rod.
  std::vector<RodEnergy*>& energies;
  /// Lambda functions that draw elements related to the integrator.
  std::vector<std::function<void(void)>> drawFuncs;
public:
  Integrator(Rod&, std::vector<RodEnergy*>&);
  virtual ~Integrator() { }
  /// Attempt to move the rod forward in time one tick. Returns false on failure.
  virtual bool integrate(Clock&)=0;
  /// Draw elements related to the integrator.
  virtual void draw();
  
};

#endif
