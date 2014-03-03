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
  // set restYarn
  Integrator i;
  
  // TODO: init workspace
  for (int i=1; i<=restYarn.numIntCPs(); i++) {
    Segment& ePrev = restYarn.segments[i-1];
    Segment& eNext = restYarn.segments[i];

    Vec3f curveBinorm = 2*ePrev.vec().cross(eNext.vec()) /
      (ePrev.length()*eNext.length() + ePrev.vec().dot(eNext.vec()));
    
    offVec<Vec2f> matCurvature(-(i-1));
    matCurvature.push_back(Vec2f(curveBinorm.dot(ePrev.m2()), -(curveBinorm.dot(ePrev.m1()))));
    matCurvature.push_back(Vec2f(curveBinorm.dot(eNext.m2()), -(curveBinorm.dot(eNext.m1()))));
    ws.restMatCurvature.push_back(matCurvature);
  }
  
  // init yarns
  curYarn = new Yarn(restYarn);
  nextYarn = new Yarn();
  
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
    
    
    // Switch yarns
    Yarn* tempYarn = curYarn;
    curYarn = nextYarn;
    nextYarn = tempYarn;
    
    done = true;
  }
  
  delete curYarn;
  delete nextYarn;
  
  std::cout << "Simulation completed successfully!\n";
  std::exit(0);
}



void Simulator::computeForces() {
  

  
}
