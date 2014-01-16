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
#define DECOMP_PATH "../../../../afghan/comp/decomp/"

#define RESOLUTION (1E-7)
#define MAX_STEPS 32767
#define THRESHOLD (RESOLUTION*MAX_STEPS)

vector<AAGroup> groupVec;

// Set the bounding box of a AAGroup.
bool getBoundingBox(const AppData& ad, AAGroup* group)
{
  if (group->indices.empty()) return false;
  group->aabb[0] = Vec3f(ad.frames[ad.currentFrame][group->indices[0]]); //min
  group->aabb[1] = Vec3f(ad.frames[ad.currentFrame][group->indices[0]]); //max
  for (int i=1; i<group->indices.size(); i++) {
    for (int j=0; j<3; j++) {
      if (ad.frames[ad.currentFrame][group->indices[i]][j] < group->aabb[0][j]) {
        group->aabb[0][j] = ad.frames[ad.currentFrame][group->indices[i]][j];
      }
      if (ad.frames[ad.currentFrame][group->indices[i]][j] > group->aabb[1][j]) {
        group->aabb[1][j] = ad.frames[ad.currentFrame][group->indices[i]][j];
      }
    }
  }
  
  return true;
}

// Divide a group recursively.
void divideGroup(const AppData& ad, const int index)
{
  float resid = getResidual(ad, groupVec[index].indices, ad.currentFrame, ad.currentFrame+1, true);
  
  if (resid < THRESHOLD && groupVec[index].indices.size() <= UCHAR_MAX) {
    groupVec[index].ok = true;
    return;
  }
  
  assert(groupVec[index].indices.size() > 1);
  
  // Get widest dimension
  int dim = (groupVec[index].aabb[1][0] - groupVec[index].aabb[0][0] > groupVec[index].aabb[1][1] - groupVec[index].aabb[0][1] ? 0 : 1);
  dim = (groupVec[index].aabb[1][2] - groupVec[index].aabb[0][2] > groupVec[index].aabb[1][dim] - groupVec[index].aabb[0][dim] ? 2 : dim);
  
  // Split the volume in half.
  float split = (groupVec[index].aabb[0][dim] + groupVec[index].aabb[1][dim])/2;
  
  groupVec.push_back(AAGroup());
  int sibling = groupVec.size()-1;
  assert(groupVec[sibling].indices.empty());
  vector<uint32_t> newIndices;
  for (uint32_t index : groupVec[index].indices) {
    if (ad.frames[ad.currentFrame][index][dim] < split) {
      newIndices.push_back(index);
    } else {
      groupVec[sibling].indices.push_back(index);
    }
  }
  groupVec[index].indices = newIndices;
  
  getBoundingBox(ad, &groupVec[index]);
  getBoundingBox(ad, &groupVec[sibling]);
  
}

// Write a .key file.
int writeKeyframeFile(const AppData& ad)
{
  if (groupVec.empty()) {
    cerr << "Tried to write keyframe with empty group\n";
    exit(1);
  }
  
  stringstream outFileName;
  outFileName << COMP_PATH << ad.currentFrame << ".key";
  ofstream outFile(outFileName.str(), ios::trunc | ios::binary);
  
  if (!outFile) {
    cerr << "Unable to create outfile: " << outFileName.str() << "\n";
  }
  
  int bytesWritten = 0;
  
  for (AAGroup group : groupVec) {
    assert(group.indices.size() <= UCHAR_MAX);
    outFile << (unsigned char)group.indices.size();
    bytesWritten++;
    for (uint32_t index : group.indices) {
      for (int i=0; i<3; i++) {
        writeBinary(&(ad.frames[ad.currentFrame][index][i]), sizeof(float), &outFile);
        bytesWritten += sizeof(float);
      }
    }
  }
  
  outFile.close();
  return bytesWritten;
}

// Write a .comp compressed file.
int writeCompressedFrame(const AppData& ad)
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
  for (AAGroup group : groupVec) {
    // get transformation
    Vec3f avgPrev;
    Vec3f avg;
    for (uint32_t index : group.indices) {
      avgPrev += ad.frames[ad.currentFrame-1][index];
      avg += ad.frames[ad.currentFrame][index];
    }
    avg /= group.indices.size();
    avgPrev /= group.indices.size();
    Matrix<double, 3, 3> M;
    if (group.indices.size() == 1) {
      M.Identity();
    } else {
      
      Matrix<double, 3, 3> A = Array33d::Zero();
      Matrix<double, 3, 3> B = Array33d::Zero();
      
      for (uint32_t index : group.indices) {
        for (int i=0; i<3; i++) {
          for (int j=0; j<3; j++) {
            Vec3f pPrev = ad.frames[ad.currentFrame-1][index] - avgPrev;
            A(i,j) += pPrev[i]*pPrev[j];
            Vec3f p = ad.frames[ad.currentFrame][index] - avg;
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
      
      writeBinary(&(Mf(i)), sizeof(float), &outFile);
      bytesWritten += sizeof(float);
    }
    for (int i=0; i<3; i++) {
      writeBinary(&(avg[i]), sizeof(float), &outFile);
      bytesWritten += sizeof(float);
    }
    
    // transform points, get number of bad transformations
    vector<int> badIndices;
    vector<int16_t> corrections;
    vector<float> badPos;
    Vector3f p;
    Vector3f q;
    for (int i=0; i<group.indices.size(); i++) {
      for (int j=0; j<3; j++) {
        p(j) = ad.frames[ad.currentFrame-1][group.indices[i]][j] - avgPrev[j];
      }
      q = Mf * p;
      for (int j=0; j<3; j++) {
        float qNew = q[j] + avg[j];
        float actual = ad.frames[ad.currentFrame][group.indices[i]][j];
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
    outFile << (unsigned char)badIndices.size();
    bytesWritten++;
    for (int index : badIndices) {
      if(index >= UCHAR_MAX) {
        outFile.close();
        return INT32_MAX;
      }
      outFile << (unsigned char)index;
      bytesWritten++;
    }
    
    auto correctionsIter = corrections.begin();
    auto badPosIter = badPos.begin();
    for (int i=0; i<corrections.size() + badPos.size(); i++) {
      if(std::binary_search(badIndices.begin(), badIndices.end(), i)) {
        writeBinary(&(*badPosIter), sizeof(float), &outFile);
        bytesWritten += sizeof(float);
        ++badPosIter;
      } else {
        writeBinary(&(*correctionsIter), sizeof(int16_t), &outFile);
        bytesWritten += sizeof(int16_t);
        ++correctionsIter;
      }
    }
  }
  cout << "bad: " << numbad <<"\n";
  outFile.close();
  return bytesWritten;
}

// Compress a range of frames. Create a .comp file if possible; otherwise create a .key file.
void compress(AppData& ad)
{
  while (ad.currentFrame <= END_FRAME) {
    if (ad.frames.left_buffer_size(ad.currentFrame) < 1) {
      for (int i=1; i<NUM_FRAMES-2 && ad.currentFrame + i <= END_FRAME; i++) {
        loadFrame(ad, ad.currentFrame + i, true);
      }
    }
    int size;
    if (ad.currentFrame != START_FRAME) {
      cout << "Compressing frame " << ad.currentFrame << "...\n";
      size = writeCompressedFrame(ad);
    }
    if (ad.currentFrame == START_FRAME || size > ad.frames[0].size()*3*sizeof(float)) {
      cout << "Compression failed. Making keyframe...";
      stringstream s;
      s << COMP_PATH << ad.currentFrame << ".comp";
      remove(s.str().c_str());
      
      groupVec.clear();
      groupVec.push_back(AAGroup());
      for (int i=0; i<ad.frames[0].size(); i++) {
        groupVec[0].indices.push_back(i);
      }
      getBoundingBox(ad, &(groupVec[0]));
      bool done = false;
      while (!done) {
        done = true;
        for (int i=0; i<groupVec.size(); i++) {
          if (!groupVec[i].ok) {
            done = false;
            divideGroup(ad, i);
          }
        }
      }
      cout << "num nodes: " << groupVec.size() << "\n";
      writeKeyframeFile(ad);
    }
    ad.currentFrame++;
  }
  
  
}

// TODO: don't load this in; just convert it to a .pos file
// Decompress a given frame (.comp or .key) and load it into
// the circular buffer.
void decompressFrame(AppData& ad, const int frame)
{
  cout << "Decompress frame " << frame << "\n";
  stringstream posFilename;
  posFilename << COMP_PATH << frame << ".key";
  ifstream posFile(posFilename.str(), ios::binary);
  
  if (posFile) { // Decompress keyframe
    groupVec.clear();
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
      unsigned char numPoints = in[0];
      groupVec.push_back(AAGroup());
      for (int i=0; i<numPoints; i++) {
        posFile.read(in, 3*sizeof(float));
        bytesRead += 3*sizeof(float);
        newPoints.push_back(Vec3f(*(float*)(&in[0]), *(float*)(&in[sizeof(float)]), *(float*)(&in[2*sizeof(float)])));
        groupVec[groupVec.size()-1].indices.push_back(counter);
        counter++;
      }
    }
    
    Frame f = Frame(newPoints);
    ad.frames.push_back(f);
    
    posFile.close();
    return;
  }
  
  posFilename = stringstream();
  posFilename << COMP_PATH << frame << ".comp";
  posFile = ifstream(posFilename.str(), ios::binary);
  
  if (posFile) { // Decompress compressed file
    vector<Vec3f> newPoints;
    
    for (AAGroup group : groupVec) {
      Matrix33f T;
      char in[12*sizeof(float)];
      posFile.read(in, 12*sizeof(float));
      for (int i=0; i<9; i++) {
        T[i] = *(float*)(&in[i*sizeof(float)]);
      }
      Vec3f avg = Vec3f(*(float*)(&in[9*sizeof(float)]), *(float*)(&in[10*sizeof(float)]), *(float*)(&in[11*sizeof(float)]));
      
      Vec3f avgPrev;
      for (int index : group.indices) {
        avgPrev += ad.frames[frame-1][index];
      }
      avgPrev /= group.indices.size();
      
      posFile.read(in, sizeof(char));
      unsigned char numBad = in[0];
      vector<uint> badIndices;
      for (int i=0; i<numBad; i++) {
        posFile.read(in, sizeof(char));
        unsigned char badIndex = in[0];
        badIndices.push_back(badIndex);
      }
      
      int counter = 0;
      for (int index : group.indices) {
        Vec3f newPoint = ad.frames[frame-1][index] - avgPrev;
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
    ad.frames.push_back(f);
    
    posFile.close();
  } else {
    cerr << "Error: no compressed frame found: " << posFilename.str() << "\n";
  }
}