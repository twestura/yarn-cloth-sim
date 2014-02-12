//
//  Simulator.cpp
//  Visualizer
//
//  Created by eschweickart on 2/10/14.
//
//

#include "Simulator.h"

// Implementations for Simulator
void Simulator::run() {
  std::cout << "Starting simulation... \n";
  
  // TODO: import initial model/constraints
  // set restYarn
  
  // TODO: init workspace
  for (int i=0; i<restYarn.numIntCPs(); i++) {
    offVec<Vec3f> a(-i);
    for (int j=0; j<3; j++)
      a.push_back(Vec3f(0, 0, 0));
    ws.gradCurveBinorm.push_back(a);
    ws.curvatureBinormals.push_back(Vec3f(0, 0, 0));
  }
  
  // TODO: init yarns
  curYarn = new Yarn();
  
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
  
  delete curYarn;
  
  std::cout << "Simulation completed successfully!\n";
  std::exit(0);
}



void Simulator::computeForces() {
  /*
   * calculate bending energy
   * calculate twisting energy
   *
   */
  
#ifndef PARALLEL
  // NON-PARALLEL:
  // Precomputation
  for (int i=1; i<=curYarn->numIntCPs(); i++) {
    Vec3f ePrev = curYarn->segments[i-1].vec();
    Vec3f eNext = curYarn->segments[i].vec();
    float denominator = restYarn.segments[i-1].length()*restYarn.segments[i].length() +
    ePrev.dot(eNext);
    ws.curvatureBinormals[i] = 2*ePrev.cross(eNext) / denominator;
    
    Eigen::Matrix3f bracketEprev;
    bracketEprev << 0, -ePrev[2], ePrev[1], ePrev[2], 0, -ePrev[0], -ePrev[1], ePrev[0], 0;
    Eigen::Matrix3f bracketEnext;
    bracketEnext << 0, -eNext[2], eNext[1], eNext[2], 0, -eNext[0], -eNext[1], eNext[0], 0;
    
    Eigen::Matrix3f outerProductWithNext = ws.curvatureBinormals[i]*(eNext.transpose());
    Eigen::Matrix3f outerProductWithPrev = ws.curvatureBinormals[i]*(ePrev.transpose());
    
    ws.gradCurveBinorm[i][i-1] = (2*bracketEnext+outerProductWithNext) / denominator;
    ws.gradCurveBinorm[i][i+1] = (2*bracketEprev+outerProductWithPrev) / denominator;
    ws.gradCurveBinorm[i][i] = -(ws.gradCurveBinorm[i][i-1] + ws.gradCurveBinorm[i][i+1]);
  }
  
  // Get forces for each point
  for (int i=0; i<curYarn->numCPs(); i++) {
    Vec3f force(0, 0, 0);
    
    // Bending forces
    for (int k=i-1; k<=i+1; k++) {
      if (k < 1 || k >= curYarn->numCPs()-1) continue;
      Vec3f forcePart(0, 0, 0);
      for (int j=k-1; j<=k; j++) {
        Vec3f m1 = curYarn->segments[j].m1();
        Vec3f m2 = curYarn->segments[j].m2();
        Eigen::Matrix<float,2,3> matMatrix;
        matMatrix << m1[0], m1[1], m1[2], m2[0], m2[1], m2[2];
        Vec2f gradMatCurvature = matMatrix * ws.gradCurveBinorm[i][k];
        
        // TODO: the bending matrix may depend on the edge.
        Eigen::Matrix2f bendMatrix;
        bendMatrix.Identity();
        
        Vec2f matCurvature;
        matCurvature << ws.curvatureBinormals[i].dot(m2), -(ws.curvatureBinormals[i].dot(m1));
        
        Vec2f restMatCurvature; // TODO: FIXME figure out how to set this
        
        forcePart += gradMatCurvature.transpose() * bendMatrix * (matCurvature - restMatCurvature);
      }
      force += forcePart / (curYarn->segments[k-1].length() + curYarn->segments[k].length());
    }
    
    // Twisting forces
    if (i > 0 && i < curYarn->numIntCPs()) {
      Vec3f prevRefTwist = ws.curvatureBinormals[i] / (2*restYarn.segments[i-1].length());
      Vec3f nextRefTwist = ws.curvatureBinormals[i] / (2*restYarn.segments[i].length());
      Vec3f derivRefTwist = -(prevRefTwist+nextRefTwist);
      
      // TODO: is this right?
      force += derivRefTwist;
    }
    
    // add in external forces
    
    nextYarn->points[i].force = force;
  }
  
  
  
#else
  // PARALLEL
  // Rename some variables to allow capture
  const Yarn& restYarnRef = restYarn;
  const Yarn& curYarnRef = *curYarn;
  std::vector<Vec3f>& cbsRef = ws.curvatureBinormals;
  auto f = [&restYarnRef, &curYarnRef, &cbsRef](size_t offset, size_t size) -> void {
    for (size_t i=offset; i<size; i++) {
      Vec3f ePrev = curYarnRef.segments[i].vec();
      Vec3f eNext = curYarnRef.segments[i+1].vec();
      cbsRef[i] = 2*ePrev.cross(eNext) /
                  (restYarnRef.segments[i].length()*restYarnRef.segments[i+1].length() +
                   ePrev.dot(eNext));
    }
  };
  
  MAYBE_PARALLEL(f, ws.threads, curYarnRef.numIntCPs(), 1000)
  /*
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
   */
#endif // ifndef PARALLEL
  
}