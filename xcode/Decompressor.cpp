//
//  Decompressor.cpp
//  Visualizer
//
//  Created by eschweickart on 1/27/14.
//
//

#include <fstream>
#include "cinder/gl/GlslProg.h"
#include "cinder/gl/Vbo.h"


#include "Resources.h"
#include "Common.h"
#include "Decompressor.h"
#include "BetterComp.h"

using namespace ci;

void Decompressor::init() {
  if (isInit) return;
  assert(MaxProxyJoints == 100); // Reprogram the vertex shader if you need to change this
  cout << "Starting compressed mode...\n";
  
  // Load compressed file
  stringstream inFileName;
  inFileName << CompPath << "afghan.comp"; // TODO: generalize
  ifstream inFile(inFileName.str(), ios::binary);
  if (!inFile) {
    cerr << "Error opening compressed file: " << inFileName.str() << "\n";
  }

  int numPoints;
  readBinary(numPoints, inFile);
  int numPJs;
  readBinary(numPJs, inFile);
  assert(numPJs == MaxProxyJoints);
  
  std::vector<Vec3f> points;
  std::vector<Vec4f> pji;
  std::vector<Vec4f> pjw;
  points.reserve(numPoints);
  pji.reserve(numPoints);
  pjw.reserve(numPoints);
  
  for (int i=0; i<numPoints; i++) {
    Vec3f newPoint;
    for (int j=0; j<3; j++) {
      readBinary(newPoint[j], inFile);
    }
    points.push_back(newPoint);
    
    uchar numPJs;
    readBinary(numPJs, inFile);
    Vec4f newpji;
    Vec4f newpjw;
    for (int j=0; j<numPJs; j++) {
      readBinary(newpji[j], inFile);
      newpji[j]++;
      readBinary(newpjw[j], inFile);
    }
    pji.push_back(newpji);
    pjw.push_back(newpjw);
  }
  
  int totalFrames;
  readBinary(totalFrames, inFile);
  std::vector<Matrix33f*>& mats = frames.first;
  std::vector<Vec3f*>& vecs = frames.second;
  assert(mats.empty() && vecs.empty());
  mats.reserve(totalFrames);
  vecs.reserve(totalFrames);
  for (int i=0; i<totalFrames; i++) {
    mats.push_back(new Matrix33f[MaxProxyJoints+1]);
    vecs.push_back(new Vec3f[MaxProxyJoints+1]);
    
    mats[i][0] = Matrix33f::identity();
    vecs[i][0] = Vec3f::zero();
    for (int j=1; j<MaxProxyJoints+1; j++) {
      for (int k=0; k<9; k++) {
        readBinary(mats[i][j][k], inFile);
      }
      for (int k=0; k<3; k++) {
        readBinary(vecs[i][j][k], inFile);
      }
    }
  }
  
  try {
    skinningProg = gl::GlslProg(cinder::app::loadResource(RES_VERT_GLSL), cinder::app::loadResource(RES_FRAG_GLSL));
  } catch (gl::GlslProgCompileExc e) {
    cerr << "Error compiling program: " << e.what();
    exit(1);
  }
  gl::VboMesh::Layout layout;
  layout.setStaticPositions();
  layout.setStaticIndices();
  mesh = gl::VboMesh::create(points.size(), 2*points.size()-1, layout, GL_LINES);
  
  mesh->bufferPositions(points);
  
  vector<uint32_t> indexBuffer;
  indexBuffer.reserve(2*(points.size()-1));
  for(int i=0; i<2*(points.size()-1); i++) {
    indexBuffer.push_back(i);
    indexBuffer.push_back(i+1);
  }
  mesh->bufferIndices(indexBuffer);

  gl::VboMesh::VertexIter vIter = mesh->mapVertexBuffer();
  for (int i=0; i<points.size(); i++) {
    vIter.setCustomVec4f(skinningProg.getAttribLocation("proxyJointIndices"), pji[i]);
    vIter.setCustomVec4f(skinningProg.getAttribLocation("proxyJointWeights"), pjw[i]);
    ++vIter;
  }
  
  currentFrame = 0;
  inFile.close();
  isInit = true;
}

void Decompressor::draw() {
  if (!isInit) init();
  if (currentFrame >= 0 && currentFrame < frames.first.size()) {
    if (currentFrame == 0) {
      Matrix33f mats[MaxProxyJoints+1];
      Vec3f vecs[MaxProxyJoints+1];
      skinningProg.uniform("worldProxyJoints", mats, MaxProxyJoints+1);
      skinningProg.uniform("worldProxyJointsT", vecs, MaxProxyJoints+1);
    } else {
      skinningProg.uniform("worldProxyJoints", frames.first[currentFrame], MaxProxyJoints+1);
      skinningProg.uniform("worldProxyJointsT", frames.second[currentFrame], MaxProxyJoints+1);
    }
    skinningProg.bind();
    gl::draw(mesh);
    skinningProg.unbind();
  } else {
    cerr << "Warning: compressed frame out of bounds\n";
  }
}

void Decompressor::clear() {
  for (Matrix33f* mat : frames.first) {
    delete[] mat;
  }
  frames.first.clear();
  
  for (Vec3f* vec : frames.second) {
    delete[] vec;
  }
  frames.second.clear();
  isInit = false;
}

Decompressor::~Decompressor() {
  clear();
}