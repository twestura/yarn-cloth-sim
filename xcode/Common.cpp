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
    char in[sizeof(float)];
    for (int i=0; i<3; i++) {
      posFile.read(in, sizeof(float));
      point[i] = *(float*)&in;
    }
    newPoints.push_back(Vec3f(point[0], point[1], point[2]));
    
    if (residFile) {
      residFile.read(in, sizeof(float));
      newResid.push_back(*(float*)&in);
      
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
float getResidual(const AppData& ad, const vector<uint32_t> indices, const int curFrameNum, const int nextFrameNum, const bool retMax)
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

// A fairly unsafe method that writes contiguous memory to a given file.
// Abstracted here to confine the dragons and black magic.
void writeBinary(const void* data, const uint size, ofstream* outFile)
{
  char* out = (char*) data;
  for (int i=0; i<size; i++) {
    (*outFile) << out[i];
  }
}