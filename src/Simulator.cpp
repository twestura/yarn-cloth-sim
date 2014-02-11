//
//  Simulator.cpp
//  Visualizer
//
//  Created by eschweickart on 2/10/14.
//
//

#include "Simulator.h"
#include <boost/thread/thread.hpp>

// Implementations for Segment
Segment::Segment(const CtrlPoint& a, const CtrlPoint& b) : first(a), second(b) {}
const Vec3f inline Segment::vec() const { return second.pos - first.pos; }


// Implementations for Yarn
const size_t inline Yarn::numCPs() const { return points.size(); }
const size_t inline Yarn::numIntCPs() const { return fmax(points.size()-2, 0); }
const size_t inline Yarn::numSegs() const { return segments.size(); }


// Implementations for Simulator
void Simulator::run() {
  std::cout << "Starting simulation... \n";
  
  // TODO: import initial model/constraints
  
  
  
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
    computeForces();
    
    done = true;
  }
  
  std::cout << "Simulation completed successfully!\n";
  std::exit(0);
}



void Simulator::computeForces() {
  /*
   * calculate bending energy
   * calculate twisting energy
   *
   */
  
  // Rename some variables to allow capture
  const Yarn& yarnRef = y;
  std::vector<Vec3f>& cbsRef = ws.curvatureBinormals;
  auto f = [&yarnRef, &cbsRef](size_t offset, size_t size) -> void {
    for (size_t i=offset; i<size; i++) {
      Vec3f ePrev = yarnRef.segments[i].vec();
      Vec3f eNext = yarnRef.segments[i+1].vec();
      cbsRef[i] = 2*ePrev.cross(eNext) /
                  (yarnRef.segments[i].restLength*yarnRef.segments[i+1].restLength +
                   ePrev.dot(eNext));
    }
  };
  
  if (y.numIntCPs() > 1000) { // Are threads really helpful here?
    for (int i=0; i<NUM_THREADS; i++) {
      ws.threads[i] = boost::thread(f,
                                    (i*y.numIntCPs())/NUM_THREADS,
                                    ((i+1)*y.numIntCPs())/NUM_THREADS);
    }
    for (int i=0; i<NUM_THREADS; i++) {
      ws.threads[i].join();
    }
  } else { // Just do the work on the main thread
    f(0, y.numIntCPs());
  }
  
}