//
//  Integrator.cpp
//  Visualizer
//
//  Created by eschweickart on 6/10/14.
//
//

#include "Integrator.h"

Integrator::Integrator(Yarn& y, std::vector<YarnEnergy*>& energies) : y(y), energies(energies) { }

void Integrator::draw() {
  for (std::function<void(void)> f : frames) {
    f();
  }
  frames.clear();
}