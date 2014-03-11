//
//  Simulator.cpp
//  Visualizer
//
//  Created by eschweickart on 2/10/14.
//
//

#include "Simulator.h"
#include "Clock.h"
#include "Constants.h"

// Implementations for Simulator
void Simulator::run() {
  std::cout << "Starting simulation... \n";
  
  // TODO: import initial model/constraints

  // TODO: create energies
  
  // TODO: init integrator
  
  /*
   * precompute omegaBar^i_j ((2) of Bergou)
   * set quasistatic material frame
   * while simulating do
   *   compute forces on centerline
   *   integrate centerline
   *   satisfy glue constraints
   *   enforce inextensibility
   *   collision detection/response
   *   update Bishop frame (parallel transport through time)
   *   update quasistatic material frame
   */
  
  bool done = false;

  // MAIN SIMULATION LOOP
  while (!done) {
    // Compute the forces on the centerline
    computeForces();
    
    // Integrate
    
    done = true;
  }
  
  std::cout << "Simulation completed successfully!\n";
  std::exit(0);
}



void Simulator::computeForces() {
  

  
}
