//
//  BetterComp.cpp
//  Visualizer
//
//  Created by eschweickart on 1/21/14.
//


#include <iostream>
#include <fstream>
#include <algorithm>
#include "BetterComp.h"
#include "Common.h"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "Eigen/SparseQR"
#include "cinder/KdTree.h"

using namespace ci;

Compressor::Compressor(AppData& inad, KdTree<Vec3f, 3, NeighborLookupProc>& inkdt) : ad(inad), kdt(inkdt) { }

void Compressor::compress(const int frame) {
  if (!initialzed) {
    initialize();
    initialzed = true;
  }

  // compress
  stringstream outFileName;
  outFileName << CompPath << "afghan-test.comp"; // TODO: generalize ///DEBUG
  ofstream outFile(outFileName.str(), ios::trunc | ios::binary);
  if (!outFile) {
    cerr << "Unable to create outfile: " << outFileName.str() << "\n";
  }
  
  // write first frame
  cout << "Compressing frame 0...\n";
  int numPoints = ad.frames[0].size();
  writeBinary(&numPoints, sizeof(int), outFile);
  int numPJs = proxyJoints.size();
  writeBinary(&numPJs, sizeof(int), outFile);
  for (int i=0; i<numPoints; i++) {
    for (int j=0; j<3; j++)
      writeBinary(&ad.frames[0][i][j], sizeof(float), outFile);
    uchar numPJs = vertices[i].proxyJoints.size();
    writeBinary(&numPJs, sizeof(uchar), outFile);
    for (int j=0; j<numPJs; j++) {
      writeBinary(&vertices[i].proxyJoints[j], sizeof(int), outFile); // TODO: compress this int
      writeBinary(&vertices[i].weights[j], sizeof(float), outFile);
    }
  }
  
  writeBinary(&frame, sizeof(int), outFile);
  // write compframe
  for (int f=1; f<=frame; f++) {
    loadFramesIfNecessary(ad, Direction::Right, f);
    cout << "Compressing frame " << f << "...\n";
    vector<Matrixf34> trans = ls(f);
    
    ///DEBUG
    cout << trans[4] << "\n";
    ///
    
    for (Matrixf34 m : trans) {
      for (int i=0; i<12; i++) {
        writeBinary(&m(i), sizeof(float), outFile);
      }
    }
  }

  cout << "Compression completed successfully!\n";
  outFile.close();
}

// TODO: clean up min/max calculations
void Compressor::initialize() {
  cout << "Initializing compressor...";
  
  // Distribute proxy-joints greedily
  proxyJoints.clear();
  proxyJoints.push_back(0);
  for (int i=1; i<MaxProxyJoints; i++) {
    float maxDist = 0;
    int bestCenter;
    for (int i=1; i<ad.frames[0].size(); i++) {
      float myMaxDist = 0;
      for (int pjIndex : proxyJoints) {
        float dist = (ad.frames[0][i] - ad.frames[0][pjIndex]).length();
        if (dist > myMaxDist)
          myMaxDist = dist;
      }
      if (myMaxDist > maxDist) {
        maxDist = myMaxDist;
        bestCenter = i;
      }
    }
    assert(maxDist > 0);
    proxyJoints.push_back(bestCenter);
  }
  
  // Assign vertex weights
  
  float maxDist = 0;
  for (int i=1; i<ad.frames[0].size(); i++) {
    float myMaxDist = 0;
    for (int pjIndex : proxyJoints) {
      float dist = (ad.frames[0][i] - ad.frames[0][pjIndex]).length();
      if (dist > myMaxDist)
        myMaxDist = dist;
    }
    if (myMaxDist > maxDist)
      maxDist = myMaxDist;
  }
  
  const float influence = maxDist * ProxyJointInfluence;
  
  vertices.clear();
  for (int i=0; i<ad.frames[0].size(); i++) {
    vertices.push_back(Vertex());
  }
  
  for (int i=0; i<proxyJoints.size(); i++) {
    int pjIndex = proxyJoints[i];
    NeighborLookupProc nlp;
    kdt.lookup(ad.frames[0][pjIndex], nlp, influence);
    for (int vIndex : nlp.neighbors) {
      float weight = (ad.frames[0][vIndex] - ad.frames[0][pjIndex]).length() / influence;
      if (vertices[vIndex].proxyJoints.size() < 4) {
        vertices[vIndex].proxyJoints.push_back(i);
        vertices[vIndex].weights.push_back(weight);
      } else {
        float min = INFINITY;
        int minIndex;
        for (int j=0; j<4; j++) {
          if (vertices[vIndex].weights[j] < min) {
            min = vertices[vIndex].weights[j];
            minIndex = j;
          }
        }
        if (min < weight) {
          vertices[vIndex].proxyJoints[minIndex] = i;
          vertices[vIndex].weights[minIndex] = weight;
        }
      }
    }
  }
  
  // Normalize weights
  for (Vertex v : vertices) {
    float sum = 0;
    for (float w : v.weights)
      sum += w;
    if (sum > 1) {
      for (int i=0; i<v.weights.size(); i++)
        v.weights[i] /= sum;
    }
  }
  
  cout << "Done!\n";
}

vector<Matrixf34> Compressor::ls(const int frame) {
  using namespace Eigen;
  vector<Triplet<double>> v;
  v.reserve(12*ad.frames[0].size()); // TODO: how much do we reserve?
  for (int i=0; i<vertices.size(); i++) {
    Vertex& vertex = vertices[i];
    for (int j=0; j<vertex.proxyJoints.size(); j++) {
      v.push_back(Triplet<double>(3*i, 12*vertex.proxyJoints[j], vertex.weights[j]*ad.frames[0][i][0]));
      v.push_back(Triplet<double>(3*i, 12*vertex.proxyJoints[j]+1, vertex.weights[j]*ad.frames[0][i][1]));
      v.push_back(Triplet<double>(3*i, 12*vertex.proxyJoints[j]+2, vertex.weights[j]*ad.frames[0][i][2]));
      v.push_back(Triplet<double>(3*i, 12*vertex.proxyJoints[j]+3, vertex.weights[j]));
      v.push_back(Triplet<double>(3*i+1, 12*vertex.proxyJoints[j]+4, vertex.weights[j]*ad.frames[0][i][0]));
      v.push_back(Triplet<double>(3*i+1, 12*vertex.proxyJoints[j]+5, vertex.weights[j]*ad.frames[0][i][1]));
      v.push_back(Triplet<double>(3*i+1, 12*vertex.proxyJoints[j]+6, vertex.weights[j]*ad.frames[0][i][2]));
      v.push_back(Triplet<double>(3*i+1, 12*vertex.proxyJoints[j]+7, vertex.weights[j]));
      v.push_back(Triplet<double>(3*i+2, 12*vertex.proxyJoints[j]+8, vertex.weights[j]*ad.frames[0][i][0]));
      v.push_back(Triplet<double>(3*i+2, 12*vertex.proxyJoints[j]+9, vertex.weights[j]*ad.frames[0][i][1]));
      v.push_back(Triplet<double>(3*i+2, 12*vertex.proxyJoints[j]+10, vertex.weights[j]*ad.frames[0][i][2]));
      v.push_back(Triplet<double>(3*i+2, 12*vertex.proxyJoints[j]+11, vertex.weights[j]));
    }
  }
  
  SparseMatrix<double> A(3*ad.frames[0].size(), 12*proxyJoints.size());
  A.setFromTriplets(v.begin(), v.end());
  SparseMatrix<double> AtA = A.transpose() * A;
  
  Matrix<double, Dynamic, 1> b(3*ad.frames[0].size(), 1);
  for (int i=0; i<ad.frames[0].size(); i++) {
    for (int j=0; j<3; j++) {
      b(3*i+j) = ad.frames[frame][i][j];
    }
  }
  VectorXd c = A.transpose() * b;
  
  SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solver;
  solver.compute(AtA);
  if(solver.info()!=Success) {
    cerr << "LS: A decomposition failed.\n";
    exit(1);
  }
  VectorXd x = solver.solve(c);
  if(solver.info()!=Success) {
    cerr << "LS: solving failed.\n";
    exit(1);
  }
  
  vector<Matrixf34> a;
  
  for (int i=0; i<proxyJoints.size(); i++) {
    Matrixf34 c;
    c << x(12*i), x(12*i+1), x(12*i+2), x(12*i+3),
    x(12*i+4), x(12*i+5), x(12*i+6), x(12*i+7),
    x(12*i+8), x(12*i+9), x(12*i+10), x(12*i+11);
    a.push_back(c);
  }
  
  return a;
}