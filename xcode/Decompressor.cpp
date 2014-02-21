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
#include "Compressor.h"

using namespace ci;

void Decompressor::init() {
  if (isInit) return;
  assert(MaxProxyJoints == 150); // Reprogram the vertex shader if you need to change this
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
      int index;
      readBinary(index, inFile);
      newpji[j] = index;
      readBinary(newpjw[j], inFile);
    }
    pji.push_back(newpji);
    pjw.push_back(newpjw);
    
#ifdef GPUDEBUG
    debugPoints.push_back(newPoint);
    debugpjis.push_back(newpji);
    debugpjws.push_back(newpjw);
#endif
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
    
    for (int j=0; j<MaxProxyJoints; j++) {
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
    cerr << "Error compiling GLSL program: " << e.what();
    exit(1);
  }
  gl::VboMesh::Layout layout;
  layout.setStaticPositions();
  layout.setStaticIndices();
  layout.addDynamicCustomVec4f(); // Proxy joint indices
  const int ProxyJointIndices = 0;
  layout.addDynamicCustomVec4f(); // Proxy joint weights
  const int ProxyJointWeights = 1;
  mesh = gl::VboMesh::create(points.size(), 2*points.size()-1, layout, GL_LINES);
  
  mesh->bufferPositions(points);
  
  vector<uint32_t> indexBuffer;
  indexBuffer.reserve(2*(points.size()-1));
  for(int i=0; i<2*(points.size()-1); i++) {
    indexBuffer.push_back(i);
    indexBuffer.push_back(i+1);
  }
  mesh->bufferIndices(indexBuffer);
  
  mesh->setCustomDynamicLocation(ProxyJointIndices, skinningProg.getAttribLocation("proxyJointIndices"));
  mesh->setCustomDynamicLocation(ProxyJointWeights, skinningProg.getAttribLocation("proxyJointWeights"));
  gl::VboMesh::VertexIter vIter = mesh->mapVertexBuffer();
  for (int i=0; i<points.size(); i++) {
    vIter.setCustomVec4f(ProxyJointIndices, pji[i]);
    vIter.setCustomVec4f(ProxyJointWeights, pjw[i]);
    ++vIter;
  }


#ifdef GPUDEBUG
  gl::VboMesh::Layout debugLayout;
  debugLayout.setStaticPositions();
  debugLayout.setStaticIndices();
  debugMesh = gl::VboMesh::create(points.size(), 2*points.size()-1, debugLayout, GL_LINES);
  debugMesh->bufferIndices(indexBuffer);
  debugMesh->bufferPositions(points);
#endif
  
  
  inFile.close();
  currentFrame = 0;
  viewFrame();
  isInit = true;
}

void Decompressor::draw() {
  if (!isInit) init();
  skinningProg.bind();
  gl::draw(mesh);
  skinningProg.unbind();
  
#ifdef GPUDEBUG
  gl::color(0, 1, 0, 0.6);
  gl::draw(debugMesh);
#endif
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

void Decompressor::changeFrame(Direction d) {
  if (d == Direction::Left) {
    if (currentFrame == 0)
      return;
    currentFrame--;
  } else {
    if (currentFrame == frames.first.size()-1)
      return;
    currentFrame++;
  }
  
  viewFrame();
}

void Decompressor::viewFrame() {
  skinningProg.bind();
  skinningProg.uniform("worldProxyJoints", frames.first[currentFrame], MaxProxyJoints);
  skinningProg.uniform("worldProxyJointsT", frames.second[currentFrame], MaxProxyJoints);
  skinningProg.unbind();
  
#ifdef GPUDEBUG
  std::vector<Vec3f> debugPosBuffer;
  for (int i=0; i<debugPoints.size(); i++) {
    if (currentFrame == 0) {
      debugPosBuffer.push_back(debugPoints[i]);
    } else {
      Vec3f pos = debugPoints[i] * (1-debugpjws[i][0]-debugpjws[i][1]-debugpjws[i][2]-debugpjws[i][3]);
      for (int j=0; j<4; j++) {
        Matrix33f& pjm = frames.first[currentFrame][(int)(debugpjis[i][j])];
        Vec3f& pjv = frames.second[currentFrame][(int)(debugpjis[i][j])];
        Matrix44f mat;
        mat.setColumn(0, Vec4f(pjm.at(0, 0), pjm.at(1, 0), pjm.at(2, 0), 0));
        mat.setColumn(1, Vec4f(pjm.at(0, 1), pjm.at(1, 1), pjm.at(2, 1), 0));
        mat.setColumn(2, Vec4f(pjm.at(0, 2), pjm.at(1, 2), pjm.at(2, 2), 0));
        mat.setColumn(3, Vec4f(pjv.x, pjv.y, pjv.z, 0));
        
        Vec3f temp = mat.transformVec(debugPoints[i]);
        pos += temp * debugpjws[i][j];
      }
      
      debugPosBuffer.push_back(pos);
    }
  }
  debugMesh->bufferPositions(debugPosBuffer);
#endif
}

Decompressor::~Decompressor() {
  clear();
}