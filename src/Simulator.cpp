//
//  Simulator.cpp
//  Visualizer
//
//  Created by eschweickart on 2/10/14.
//
//

#include "Simulator.h"
#include "Constants.h"

// Implementations for Simulator
void Simulator::run() {
  std::cout << "Starting simulation... \n";
  
  // TODO: import initial model/constraints
  // set restYarn
  
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
    // Get unconstrained velocities
    integrate(Backward);
    
    
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
  
#ifndef PARALLEL
  // NON-PARALLEL:
  
  // Get forces for each point
  for (int i=0; i<curYarn->numCPs(); i++) {
    Vec3f force(0, 0, 0);
    
    if (i > 0 && i < curYarn->numCPs()-1) {
      Vec3f ePrev = curYarn->segments[i-1].vec();
      Vec3f eNext = curYarn->segments[i].vec();
      float denominator = restYarn.segments[i-1].length()*restYarn.segments[i].length() +
        ePrev.dot(eNext);
      Vec3f curveBinorm = 2*ePrev.cross(eNext) / denominator;
      
      Eigen::Matrix3f bracketEprev;
      bracketEprev << 0, -ePrev[2], ePrev[1], ePrev[2], 0, -ePrev[0], -ePrev[1], ePrev[0], 0;
      Eigen::Matrix3f bracketEnext;
      bracketEnext << 0, -eNext[2], eNext[1], eNext[2], 0, -eNext[0], -eNext[1], eNext[0], 0;
      
      Eigen::Matrix3f outerProductWithNext = curveBinorm*(eNext.transpose());
      Eigen::Matrix3f outerProductWithPrev = curveBinorm*(ePrev.transpose());
      
      offVec<Eigen::Matrix3f> gradCurveBinorm(-(i-1));
      for (int k=0; k<3; k++)
        gradCurveBinorm.push_back(Eigen::Matrix3f());
      
      gradCurveBinorm[i-1] = (2*bracketEnext+outerProductWithNext) / denominator;
      gradCurveBinorm[i+1] = (2*bracketEprev+outerProductWithPrev) / denominator;
      gradCurveBinorm[i] = -(gradCurveBinorm[i-1] + gradCurveBinorm[i+1]);
    
      // Bending forces
      for (int k=i-1; k<=i+1; k++) {
        if (k < 1 || k >= curYarn->numCPs()-1) continue;
        Vec3f forcePart(0, 0, 0);
        for (int j=k-1; j<=k; j++) {
          Vec3f m1 = curYarn->segments[j].m1();
          Vec3f m2 = curYarn->segments[j].m2();
          Eigen::Matrix<float,2,3> matMatrix;
          matMatrix << m1[0], m1[1], m1[2], m2[0], m2[1], m2[2];
          Eigen::Matrix<float,2,3> gradMatCurvature = matMatrix * gradCurveBinorm[k];
        
          // TODO: the bending matrix may depend on the edge.
          Eigen::Matrix2f bendMatrix = Eigen::Matrix2f::Identity();
        
          Vec2f matCurvature(curveBinorm.dot(m2), -(curveBinorm.dot(m1)));
          Vec2f& restMatCurvature = ws.restMatCurvature[k][j];
          
          // TODO: model yarn plasticity (deform restMatCurvature)
        
          forcePart += gradMatCurvature.transpose() * bendMatrix * (matCurvature-restMatCurvature);
        }
        force += forcePart / (curYarn->segments[k-1].length() + curYarn->segments[k].length());
      }
    
      // Twisting forces
      Vec3f prevRefTwist = curveBinorm / (2*ePrev.norm());
      Vec3f nextRefTwist = curveBinorm / (2*eNext.norm());
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

// TODO: implement this
void Simulator::integrate(IntegrationType it)
{
  switch (it) {
    case Forward:
      break;
      
    case Backward:
      break;
      
    case Symplectic:
      break;
      
    default:
      break;
  }
}
