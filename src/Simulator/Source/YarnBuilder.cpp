//
//  YarnBuilder.cpp
//  Visualizer
//
//  Created by eschweickart on 5/19/14.
//
//

#include "YarnBuilder.h"
#include "Eigen/Dense"
#include <fstream>
//#include "Simulator_Prefix.pch"

void YarnBuilder::buildBraid() {
  const real height = 20.0;
  const real width = 8.0;
  const real depth = 3.0;
  const int depthLevels = 3;
  const int vertLevels = 15;
  const int horizLevels = 7;
  const int totalPoints = 43;
  const Vec3c origin(0.0, 10.0, 0.0);
  
  std::vector<Vec3c> points;
  points.reserve(totalPoints);
  
  int vert = 14;
  int horiz = -3;
  int dep = 0;
  bool down = false;
  bool right = false;
  bool front = false;
  while (points.size() < totalPoints) {
    Vec3c out;
    
    out.x = dep * (depth / depthLevels);
    dep += front ? 1 : -1;
    if (dep == -depthLevels/2 || dep == depthLevels/2) front = !front;
    
    out.y = vert * (height / vertLevels);
    if (vert == 0 || vert == vertLevels-1) {
      down = !down;
      out.x = 0;
      dep = 0;
      front = !front;
    }
    vert += down ? -1 : 1;
    
    out.z = horiz * (width / horizLevels);
    if (horiz == -horizLevels/2 || horiz == horizLevels/2) {
      right = !right;
      out.x = 0;
      dep = 0;
      front = !front;
    }
    horiz += right ? 1 : -1;
    
    out += origin;
    points.push_back(out);
  }
  
  Vec3c u = (points[1] - points[0]).cross(Vec3c(1.0, 0.0, 0.0)).normalized();
  
  std::string filename = ci::app::getAppPath().string() + "braid.yarn";
  std::ofstream out(filename, std::ios::trunc);
  assert(out.is_open() && "Error: could not open file for output");
  
  out << totalPoints << "\n";
  
  for (int i=0; i<totalPoints; i++) {
    for (int j=0; j<3; j++) {
      out << points[i][j] << "\n";
    }
  }
  for (int i=0; i<3; i++) {
    out << u[i] << "\n";
  }
  
  out.close();
  std::cout << "braid file written to " << filename << "\n";
}


