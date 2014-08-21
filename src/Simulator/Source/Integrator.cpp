//
//  Integrator.cpp
//  Visualizer
//
//  Created by eschweickart on 6/10/14.
//
//

#include "Integrator.h"

Integrator::Integrator(Rod& r, std::vector<RodEnergy*>& energies) : r(r), energies(energies) { }

void Integrator::draw() {
  for (std::function<void(void)> f : drawFuncs) {
    f();
  }
  drawFuncs.clear();
}