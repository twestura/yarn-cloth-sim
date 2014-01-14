//
//  Compressor.cpp
//  Visualizer
//
//  Created by eschweickart on 12/16/13.
//
//

#include <iostream>
#include <fstream>
#include "Compressor.h"
#include "Common.h"
#include <boost/smart_ptr.hpp>
#include "Eigen/Dense"


using namespace std;
using namespace ci;

#define COMP_PATH "../../../../afghan/comp/"

#define RESOLUTION (1E-7)
#define MAX_STEPS 32767
#define THRESHOLD (RESOLUTION*MAX_STEPS)

vector<Bvh> bvhVec;

// Set the bounding box of a bvh.
bool getBoundingBox(const AppData ad, Bvh* bvh)
{
  if (bvh->indices.empty()) return false;
  bvh->aabb[0] = Vec3f(ad.points[ad.currentFrame][bvh->indices[0]]); //min
  bvh->aabb[1] = Vec3f(ad.points[ad.currentFrame][bvh->indices[0]]); //max
  for (int i=1; i<bvh->indices.size(); i++) {
    for (int j=0; j<3; j++) {
      if (ad.points[ad.currentFrame][bvh->indices[i]][j] < bvh->aabb[0][j]) {
        bvh->aabb[0][j] = ad.points[ad.currentFrame][bvh->indices[i]][j];
      }
      if (ad.points[ad.currentFrame][bvh->indices[i]][j] > bvh->aabb[1][j]) {
        bvh->aabb[1][j] = ad.points[ad.currentFrame][bvh->indices[i]][j];
      }
    }
  }
  
  return true;
}

// Divide a bvh recursively.
void divideBvh(const AppData ad, const int index)
{
  float resid = getResidual(ad, bvhVec[index].indices, ad.currentFrame, ad.currentFrame+1, true);
  
  if (resid < THRESHOLD && bvhVec[index].indices.size() <= UCHAR_MAX) {
    bvhVec[index].ok = true;
    return;
  }
  
  assert(bvhVec[index].indices.size() > 1);
  
  // Get widest dimension
  int dim = (bvhVec[index].aabb[1][0] - bvhVec[index].aabb[0][0] > bvhVec[index].aabb[1][1] - bvhVec[index].aabb[0][1] ? 0 : 1);
  dim = (bvhVec[index].aabb[1][2] - bvhVec[index].aabb[0][2] > bvhVec[index].aabb[1][dim] - bvhVec[index].aabb[0][dim] ? 2 : dim);
  
  // TODO: sort points?
  // for now just split the volume in half. I think this makes more sense anyway.
  float split = (bvhVec[index].aabb[0][dim] + bvhVec[index].aabb[1][dim])/2;
  
  bvhVec.push_back(Bvh());
  int sibling = bvhVec.size()-1;
  assert(bvhVec[sibling].indices.empty());
  vector<uint32_t> newIndices;
  for (uint32_t index : bvhVec[index].indices) {
    if (ad.points[ad.currentFrame][index][dim] < split) {
      newIndices.push_back(index);
    } else {
      bvhVec[sibling].indices.push_back(index);
    }
  }
  bvhVec[index].indices = newIndices;
  
  getBoundingBox(ad, &bvhVec[index]);
  getBoundingBox(ad, &bvhVec[sibling]);
  
}

// Write a .key file.
int writeKeyframeFile(const AppData ad)
{
  if (bvhVec.empty()) {
    cerr << "Tried to write keyframe with empty bvh\n";
    exit(1);
  }
  
  stringstream outFileName;
  outFileName << COMP_PATH << ad.currentFrame << ".key";
  ofstream outFile(outFileName.str(), ios::trunc | ios::binary );
  
  if (!outFile) {
    cerr << "Unable to create outfile: " << outFileName.str() << "\n";
  }
  
  int bytesWritten = 0;
  
  for (Bvh bvh : bvhVec) {
    assert(bvh.indices.size() <= UCHAR_MAX);
    outFile << (char)bvh.indices.size();
    bytesWritten++;
    for (uint32_t index : bvh.indices) {
      for (int i=0; i<3; i++) {
        char* out = (char*)&(ad.points[ad.currentFrame][index][i]); // Here be dragons.
        for(int j=0; j<sizeof(float); j++) {
          outFile << out[j]; // Even more dragons...
          bytesWritten++;
        }
      }
    }
  }
  
  outFile.close();
  return bytesWritten;
}

// Write a .comp compressed file.
int writeCompressedFrame(const AppData ad)
{
  using namespace Eigen;
  
  stringstream outFileName;
  outFileName << COMP_PATH << ad.currentFrame << ".comp";
  ofstream outFile(outFileName.str(), ios::trunc | ios::binary );
  
  if (!outFile) {
    cerr << "Unable to create outfile: " << outFileName.str() << "\n";
  }
  
  int numbad = 0;
  int bytesWritten = 0;
  for (Bvh bvh : bvhVec) {
    // get transformation
    Vec3f avgPrev;
    Vec3f avg;
    for (uint32_t index : bvh.indices) {
      avgPrev += ad.points[ad.currentFrame-1][index];
      avg += ad.points[ad.currentFrame][index];
    }
    avg /= bvh.indices.size();
    avgPrev /= bvh.indices.size();
    Matrix<double, 3, 3> M;
    if (bvh.indices.size() == 1) {
      M.Identity();
    } else {
      
      Matrix<double, 3, 3> A = Array33d::Zero();
      Matrix<double, 3, 3> B = Array33d::Zero();
      
      for (uint32_t index : bvh.indices) {
        for (int i=0; i<3; i++) {
          for (int j=0; j<3; j++) {
            Vec3f pPrev = ad.points[ad.currentFrame-1][index] - avgPrev;
            A(i,j) += pPrev[i]*pPrev[j];
            Vec3f p = ad.points[ad.currentFrame][index] - avg;
            B(i,j) += pPrev[i]*p[j];
          }
        }
      }
      
      M = (A.fullPivHouseholderQr().solve(B)).transpose();
    }
    Matrix<float, 3, 3> Mf;
    
    // write transformation
    for (int i=0; i<9; i++) {
      Mf(i) = M(i);
      char* out = (char*)&(Mf(i)); // Here be dragons.
      for(int j=0; j<sizeof(float); j++) {
        outFile << out[j]; // Even more dragons...
        bytesWritten++;
      }
    }
    for (int i=0; i<3; i++) {
      char* out = (char*)&(avg[i]); // Here be dragons.
      for(int j=0; j<sizeof(float); j++) {
        outFile << out[j]; // Even more dragons...
        bytesWritten++;
      }
    }
    
    // transform points, get number of bad transformations
    vector<int> badIndices;
    vector<int16_t> corrections;
    vector<float> badPos;
    Vector3f p;
    Vector3f q;
    for (int i=0; i<bvh.indices.size(); i++) {
      for (int j=0; j<3; j++) {
        p(j) = ad.points[ad.currentFrame-1][bvh.indices[i]][j] - avgPrev[j];
      }
      q = Mf * p;
      for (int j=0; j<3; j++) {
        float qNew = q[j] + avg[j];
        float actual = ad.points[ad.currentFrame][bvh.indices[i]][j];
        if (abs(qNew - actual) > THRESHOLD) {
          numbad++;
          badIndices.push_back(3*i+j);
          badPos.push_back(actual);
        } else {
          int steps = round(((((double)actual)-qNew)/RESOLUTION));
          assert(abs(steps) <= MAX_STEPS);
          corrections.push_back((int16_t)steps);
        }
      }
    }
    
    if(badIndices.size() > UCHAR_MAX) {
      outFile.close();
      return INT32_MAX;
    }
    outFile << (char)badIndices.size();
    bytesWritten++;
    for (int index : badIndices) {
      if(index >= UCHAR_MAX) {
        outFile.close();
        return INT32_MAX;
      }
      outFile << (char)index;
      bytesWritten++;
    }
    
    auto correctionsIter = corrections.begin();
    auto badPosIter = badPos.begin();
    for (int i=0; i<corrections.size() + badPos.size(); i++) {
      if(std::binary_search(badIndices.begin(), badIndices.end(), i)) {
        char* out = (char*)&(*badPosIter); // Here be dragons.
        for(int j=0; j<sizeof(float); j++) {
          outFile << out[j]; // Even more dragons...
          bytesWritten++;
        }
        ++badPosIter;
      } else {
        char* out = (char*)&(*correctionsIter); // Here be dragons.
        for(int j=0; j<sizeof(int16_t); j++) {
          outFile << out[j]; // Even more dragons...
          bytesWritten++;
        }
        ++correctionsIter;
      }
    }
  }
  cout << "bad: " << numbad <<"\n";
  outFile.close();
  return bytesWritten;
}

// Compress a range of frames. Create a .comp file if possible; otherwise create a .key file.
void compress(AppData ad)
{
  while (ad.currentFrame <= END_FRAME) {
    if (ad.points.left_buffer_size(ad.currentFrame) < 1) {
      for (int i=1; i<NUM_FRAMES-2 && ad.currentFrame + i <= END_FRAME; i++) {
        loadFrame(&ad, ad.currentFrame + i, true); // FIXME
      }
    }
    int size;
    if (ad.currentFrame != START_FRAME) {
      cout << "Compressing frame " << ad.currentFrame << "...\n";
      size = writeCompressedFrame(ad);
    }
    if (ad.currentFrame == START_FRAME || size > ad.points[0].size()*3*sizeof(float)) {
      cout << "Compression failed. Making keyframe...";
      stringstream s;
      s << COMP_PATH << ad.currentFrame << ".comp";
      remove(s.str().c_str());
      
      bvhVec.clear();
      bvhVec.push_back(Bvh());
      for (int i=0; i<ad.points[0].size(); i++) {
        bvhVec[0].indices.push_back(i);
      }
      getBoundingBox(ad, &(bvhVec[0]));
      bool done = false;
      while (!done) {
        done = true;
        for (int i=0; i<bvhVec.size(); i++) {
          if (!bvhVec[i].ok) {
            done = false;
            divideBvh(ad, i);
          }
        }
      }
      cout << "num nodes: " << bvhVec.size() << "\n";
      writeKeyframeFile(ad);
    }
    ad.currentFrame++;
  }
  
  
}

// TODO: don't load this in; just convert it to a .pos file
// Decompress a given frame (.comp or .key) and load it into
// the circular buffer.
void decompressFrame(AppData ad, const int frame)
{
  cout << "Decompress frame " << frame << "\n";
  stringstream posFilename;
  posFilename << COMP_PATH << frame << ".key";
  ifstream posFile(posFilename.str(), ios::binary);
  
  if (posFile) { // Decompress keyframe
    bvhVec.clear();
    vector<Vec3f> newPoints;
    
    posFile.seekg(0, posFile.end);
    int length = posFile.tellg();
    posFile.seekg(0, posFile.beg);
    
    int counter = 0;
    int bytesRead = 0;
    while (bytesRead < length) {
      char in[3*sizeof(float)];
      posFile.read(in, sizeof(char));
      bytesRead++;
      int numPoints = in[0];
      if (numPoints < 0) { // Black magic to get around the absence of uchar
        numPoints &= 127;
        numPoints += 128;
      }
      bvhVec.push_back(Bvh());
      for (int i=0; i<numPoints; i++) {
        posFile.read(in, 3*sizeof(float));
        bytesRead += 3*sizeof(float);
        newPoints.push_back(Vec3f(*(float*)(&in[0]), *(float*)(&in[sizeof(float)]), *(float*)(&in[2*sizeof(float)])));
        bvhVec[bvhVec.size()-1].indices.push_back(counter);
        counter++;
      }
    }
    
    Frame f = Frame(newPoints);
    ad.points.push_back(f);
    
    posFile.close();
    return;
  }
  
  posFilename = stringstream();
  posFilename << COMP_PATH << frame << ".comp";
  posFile = ifstream(posFilename.str(), ios::binary);
  
  if (posFile) { // Decompress compressed file
    vector<Vec3f> newPoints;
    
    for (Bvh bvh : bvhVec) {
      Matrix33f T;
      char in[12*sizeof(float)];
      posFile.read(in, 12*sizeof(float));
      for (int i=0; i<9; i++) {
        T[i] = *(float*)(&in[i*sizeof(float)]);
      }
      Vec3f avg = Vec3f(*(float*)(&in[9*sizeof(float)]), *(float*)(&in[10*sizeof(float)]), *(float*)(&in[11*sizeof(float)]));
      
      Vec3f avgPrev;
      for (int index : bvh.indices) {
        avgPrev += ad.points[frame-1][index];
      }
      avgPrev /= bvh.indices.size();
      
      posFile.read(in, sizeof(char));
      int numBad = in[0];
      if (numBad < 0) { // Black magic to get around the absence of uchar
        numBad &= 127;
        numBad += 128;
      }
      vector<uint> badIndices;
      for (int i=0; i<numBad; i++) {
        posFile.read(in, sizeof(char));
        int badIndex = in[0];
        if (badIndex < 0) { // Black magic to get around the absence of uchar
          badIndex &= 127;
          badIndex += 128;
        }
        badIndices.push_back(badIndex);
      }
      
      int counter = 0;
      for (int index : bvh.indices) {
        Vec3f newPoint = ad.points[frame-1][index] - avgPrev;
        newPoint = T * newPoint;
        newPoint += avg;
        
        for (int i=0; i<3; i++) {
          if (binary_search(badIndices.begin(), badIndices.end(), counter)) {
            posFile.read(in, sizeof(float));
            newPoint[i] = *(float*)in;
          } else {
            posFile.read(in, sizeof(int16_t));
            int16_t correction = *(int16_t*)in;
            newPoint[i] += correction*RESOLUTION;
          }
          counter++;
        }
        newPoints.push_back(newPoint);
      }
    }
    
    Frame f = Frame(newPoints);
    ad.points.push_back(f);
    
    posFile.close();
  } else {
    cerr << "Error: no compressed frame found: " << posFilename.str() << "\n";
  }
}