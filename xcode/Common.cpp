//
//  Common.cpp
//  Visualizer
//
//  Created by eschweickart on 1/13/14.
//
//


#include "Common.h"
#include <iostream>
#include <fstream>
#include "Eigen/Dense"

using namespace std;
using namespace ci;

// Load a .pos file into the circular buffer.
void loadFrame(AppData& ad, const int frame, const bool back)
{
  stringstream posFilename;
  posFilename << POS_PATH << frame << ".pos";
  ifstream posFile(posFilename.str(), ios::binary);
  
  if (!posFile) {
    cerr << "Error: position file not found: " << posFilename.str() << "\n";
    exit(1);
  }
  
  stringstream residFilename;
  residFilename << RESID_PATH << frame << "-" << RADIUS << ".resid";
  ifstream residFile(residFilename.str(), ios::binary);
  
  if (!residFile) {
    cerr << "Warning: resid file not found: " << residFilename.str() << "\n";
  }
  
  vector<Vec3f> newPoints;
  vector<float> newResid;
  float minResid = INFINITY;
  float maxResid = -INFINITY;
  
  posFile.seekg(0, posFile.end);
  int length = posFile.tellg();
  posFile.seekg(0, posFile.beg);
  newPoints.reserve(length/3/sizeof(float));
  if (residFile) newResid.reserve(length/sizeof(float));
  for (int j=0; j<length/3/sizeof(float); j++) {
    float point[3];
    for (int i=0; i<3; i++) {
      readBinary(point[i], posFile);
    }
    newPoints.push_back(Vec3f(point[0], point[1], point[2]));
    
    if (residFile) {
      float resid;
      readBinary(resid, residFile);
      newResid.push_back(resid);
      
      if (newResid[newResid.size()-1] > maxResid) {
        maxResid = newResid[newResid.size()-1];
      }
      if (newResid[newResid.size()-1] < minResid) {
        minResid = newResid[newResid.size()-1];
      }
    }
  }
  
  cout << maxResid <<"\n";
  Frame f = (newResid.empty() ? Frame(newPoints) : Frame(newPoints, newResid, minResid, maxResid));
  
  if (back) {
    ad.frames.push_back(f);
  } else {
    ad.frames.push_front(f);
  }
  
  posFile.close();
  
}

// Given a group of points, compute the residual from the least squares
// best fit transformation from the current frame to the next. If retMax
// is true, the maximum error is returned; otherwise, the total error squared
// is returned.
float getResidual(const AppData& ad, const vector<uint32_t>& indices, const int curFrameNum, const int nextFrameNum, const bool retMax)
{
  using namespace Eigen;
  int n = indices.size();
  
  if (n==1) return 0;
  if (n==0) {
    cerr << "Error: getResidual called on group of size 0";
    exit(1);
  }
  
  // TODO: check that both frames are in memory
  const Frame* curFrame = &ad.frames[curFrameNum];
  const Frame* nextFrame = &ad.frames[nextFrameNum];
  
  Vec3f curAvg = Vec3f();
  Vec3f nextAvg = Vec3f();
  for (uint32_t index : indices) {
    curAvg += (*curFrame)[index];
    nextAvg += (*nextFrame)[index];
  }
  
  curAvg /= n;
  nextAvg /= n;
  
  Matrix<double, 3, 3> A = Array33d::Zero();
  Matrix<double, 3, 3> B = Array33d::Zero();
  for (uint32_t index : indices) {
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
        Vec3f p = (*curFrame)[index] - curAvg;
        A(i,j) += p[i]*p[j];
        Vec3f pNext = (*nextFrame)[index] - nextAvg;
        B(i,j) += p[i]*pNext[j];
      }
    }
  }
  
  // Q = MP
  // M = QP'[PP']^-1
  // P'M' = Q'
  // (PP')M' = PQ'
  Matrix<double, 3, 3> M = (A.fullPivHouseholderQr().solve(B)).transpose();
  
  // Compute residuals
  float ret = 0;
  Vec3f temp;
  Vector3d input;
  Vector3d output;
  for (uint32_t index : indices) {
    temp = (*curFrame)[index] - curAvg;
    input << temp.x, temp.y, temp.z;
    output = M * input;
    for (int j=0; j<3; j++) {
      float diff = abs(output(j) - ((*nextFrame)[index][j] - nextAvg[j]));
      if (retMax) {
        ret = max(ret, diff);
      } else {
        ret += diff * diff;
      }
    }
  }
  
  if (!retMax) {
    ret = sqrt(ret);
  }
  
  return ret;
  
}


// FIXME: this is broken.
float newGetResidual(const AppData& ad, const vector<uint32_t>& indices, const int curFrame, const int nextFrame, const bool retMax)
{
  using namespace Eigen;
  int n = indices.size();
  
  // Compute centroids
  Vec3f c = Vec3f::zero();
  Vec3f cTilde = Vec3f::zero();
  for (int index : indices) {
    c += ad.frames[nextFrame][index];
    cTilde += ad.frames[curFrame][index];
  }
  c /= n;
  cTilde /= n;
  
  // Create least squares system--see James and Twigg, Appendix A.
  Matrix<double, 12, 12> M;
  M << cTilde[0]*cTilde[0], cTilde[1]*cTilde[0], cTilde[2]*cTilde[0], 0, 0, 0, 0, 0, 0, cTilde[0], 0, 0,
  cTilde[0]*cTilde[1], cTilde[1]*cTilde[1], cTilde[2]*cTilde[1], 0, 0, 0, 0, 0, 0, 0, cTilde[0], 0,
  cTilde[0]*cTilde[2], cTilde[1]*cTilde[2], cTilde[2]*cTilde[2], 0, 0, 0, 0, 0, 0, 0, 0, cTilde[0],
  0, 0, 0, cTilde[0]*cTilde[0], cTilde[1]*cTilde[0], cTilde[2]*cTilde[0], 0, 0, 0, cTilde[1], 0, 0,
  0, 0, 0, cTilde[0]*cTilde[1], cTilde[1]*cTilde[1], cTilde[2]*cTilde[1], 0, 0, 0, 0, cTilde[1], 0,
  0, 0, 0, cTilde[0]*cTilde[2], cTilde[1]*cTilde[2], cTilde[2]*cTilde[2], 0, 0, 0, 0, 0, cTilde[1],
  0, 0, 0, 0, 0, 0, cTilde[0]*cTilde[0], cTilde[1]*cTilde[0], cTilde[2]*cTilde[0], cTilde[2], 0, 0,
  0, 0, 0, 0, 0, 0, cTilde[0]*cTilde[1], cTilde[1]*cTilde[1], cTilde[2]*cTilde[1], 0, cTilde[2], 0,
  0, 0, 0, 0, 0, 0, cTilde[0]*cTilde[2], cTilde[1]*cTilde[2], cTilde[2]*cTilde[2],  0, 0, cTilde[2],
  cTilde[0], cTilde[1], cTilde[2], 0, 0, 0, 0, 0, 0, 1, 0, 0,
  0, 0, 0, cTilde[0], cTilde[1], cTilde[2], 0, 0, 0, 0, 1, 0,
  0, 0, 0, 0, 0, 0, cTilde[0], cTilde[1], cTilde[2], 0, 0, 1;
  
  Matrix<double, 12, 1> B;
  B << c[0]*cTilde[0], c[1]*cTilde[0], c[2]*cTilde[0], c[0]*cTilde[1], c[1]*cTilde[1], c[2]*cTilde[1],
  c[0]*cTilde[2], c[1]*cTilde[2], c[2]*cTilde[2], c[0], c[1], c[2];
  
  // Solve least squares, convert to affine transform matrix
  Matrix<double, 12, 1> X = M.fullPivHouseholderQr().solve(B);
  
  Matrix<float, 3, 4> Transform;
  for (int i=0; i<9; i++) {
    Transform(i/3, i%3) = X(i);
  }
  Transform(0,3) = X(9);
  Transform(1,3) = X(10);
  Transform(2,3) = X(11);
  
  cout << Transform << "\n";
  
  // Compute residuals
  float ret = 0;
  Vec3f temp;
  Vector4f input;
  Vector3f output;
  for (int index : indices) {
    temp = ad.frames[curFrame][index];
    input << temp.x, temp.y, temp.z, 1;
    output = Transform * input;
    for (int j=0; j<3; j++) {
      float diff = abs(output(j) - ad.frames[nextFrame][index][j]);
      if (retMax) {
        ret = max(ret, diff);
      } else {
        ret += diff * diff;
      }
    }
  }
  
  if (!retMax) {
    ret = sqrt(ret);
  }
  
  return ret;
}

// A fairly unsafe method that writes contiguous memory to a given binary file.
// Abstracted here to confine the dragons and black magic.
void writeBinary(const void* data, const uint size, ofstream& outFile)
{
  char* out = (char*) data;
  for (int i=0; i<size; i++) {
    outFile << out[i];
  }
}

// A fairly unsafe method that reads from a given binary file.
// Abstracted here to confine the dragons and black magic.
template <class T>
void readBinary(T& data, ifstream& inFile)
{
  assert(!inFile.eof());
  char in[sizeof(T)];
  inFile.read(in, sizeof(T));
  data = *(T*)(&in);
}

// Actually generate symbols for specific instances of readBinary.
// FIXME: There has got to be a better way to do this. This kind of
// defeats the point of having a template function. Stupid linker.
void fakeyfakestupidmethoddontcallthis() {
  ifstream stupid("fakeyfakefile.txt");
  int wowthisisdumb;
  readBinary(wowthisisdumb, stupid);
  unsigned char sodumb;
  readBinary(sodumb, stupid);
}

// Load more frames into the circular buffer if necessary
void loadFramesIfNecessary(AppData& ad, const Direction d, const int frame) {
  if (d == Direction::Right) {
    if (ad.frames.left_buffer_size(frame) < 1 && END_FRAME != ad.currentFrame) {
      for (int i=1; i<NUM_FRAMES-2 && frame + i <= END_FRAME; i++) {
        loadFrame(ad, frame + i, true);
      }
    }
  } else {
    if (ad.frames.right_buffer_size(frame) < 1 && START_FRAME != frame && 1 != frame) {
      for (int i=1; i<NUM_FRAMES-2 && frame - i >= START_FRAME && frame - i >= 1; i++) {
        loadFrame(ad, frame - i, false);
      }
    }
  }
}

