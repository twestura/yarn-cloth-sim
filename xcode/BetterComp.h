//
//  BetterComp.h
//  Visualizer
//
//  Created by eschweickart on 1/21/14.
//
//

#ifndef __Visualizer__BetterComp__
#define __Visualizer__BetterComp__

#include <iostream>
#include "Eigen/Dense"
#include "Common.h"
#include "cinder/KdTree.h"


typedef unsigned char uchar;
typedef Eigen::Matrix<float, 3, 4> Matrixf34;
typedef Eigen::Matrix<double, 3, 3> Matrixd33;

struct Vertex {
  std::vector<int> proxyJoints;
  std::vector<float> weights;
};

const int MaxProxyJoints = 100;
const float ProxyJointInfluence = 1.5;

const double Resolution = .002;
const int MaxSteps = 127;
const double Threshold = Resolution * MaxSteps;
const std::string CompPath = "../../../../afghan/comp/";

class Compressor {
  void initialize();
  std::vector<Matrixf34> ls(const int);
  
  AppData& ad;
  KdTree<Vec3f, 3, NeighborLookupProc>& kdt;
  bool initialzed = false;
  std::vector<int> proxyJoints;
  std::vector<Vertex> vertices;
  
public:
  Compressor(AppData&, KdTree<Vec3f, 3, NeighborLookupProc>&);
  
  void compress(const int);
};


#endif /* defined(__Visualizer__BetterComp__) */
