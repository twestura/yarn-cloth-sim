//
//  Decompressor.h
//  Visualizer
//
//  Created by eschweickart on 1/27/14.
//
//

#ifndef __Visualizer__Decompressor__
#define __Visualizer__Decompressor__

// #define GPUDEBUG

#include <iostream>
#include "cinder/gl/GlslProg.h"

using namespace ci;

class Decompressor {
  bool isInit = false;
  int currentFrame;
  std::pair<std::vector<Matrix33f*>, std::vector<Vec3f*>> frames;
  gl::GlslProg skinningProg;
  gl::VboMeshRef mesh;
  
#ifdef GPUDEBUG
  gl::VboMeshRef debugMesh;
  std::vector<Vec3f> debugPoints;
  std::vector<Vec4f> debugpjis;
  std::vector<Vec4f> debugpjws;
#endif
  
  void viewFrame();
  
  
public:
  void init();
  void draw();
  void clear();
  
  void changeFrame(Direction);
  ~Decompressor();
  
};

#endif /* defined(__Visualizer__Decompressor__) */
