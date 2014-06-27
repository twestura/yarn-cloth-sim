//
//  Integrator.h
//  Visualizer
//
//  Created by eschweickart on 6/10/14.
//
//

#ifndef Visualizer_Integrator_h
#define Visualizer_Integrator_h

#include "Yarn.h"
#include "Clock.h"
#include "Energy.h"

class Integrator {
protected:
  Yarn& y;
  std::vector<YarnEnergy*>& energies;
  std::vector<std::function<void(void)>> frames;
public:
  Integrator(Yarn&, std::vector<YarnEnergy*>&);
  virtual ~Integrator() { }
  virtual bool integrate(Clock&)=0;
  virtual void draw();
  
};

#endif
