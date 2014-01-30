//
//  Decompressor.h
//  Visualizer
//
//  Created by eschweickart on 1/27/14.
//
//

#ifndef __Visualizer__Decompressor__
#define __Visualizer__Decompressor__

#include <iostream>
#include "cinder/gl/GlslProg.h"

using namespace ci;

class Decompressor {
  bool isInit = false;
  std::pair<std::vector<Matrix33f*>, std::vector<Vec3f*>> frames;
  gl::GlslProg skinningProg;
  gl::VboMeshRef mesh;
  
  ///DEBUG
  gl::VboMeshRef debugMesh;
  std::vector<Vec3f> debugPoints;
  std::vector<Vec4f> debugpjis;
  std::vector<Vec4f> debugpjws;
  ///
  
  
public:
  int currentFrame;
  void init();
  void draw();
  void clear();
  ~Decompressor();
  
};

#endif /* defined(__Visualizer__Decompressor__) */
