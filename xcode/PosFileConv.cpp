//
//  PosFileConv.cpp
//  Visualizer
//
//  Created by eschweickart on 12/11/13.
//
//

#include "PosFileConv.h"
#include <iostream>
#include <fstream>
#include "cinder/Xml.h"

#define POS_PATH "../../../../afghan/pos/"
#define GARMENT_PATH "data/item_afghan_onplane"
#define DIR_PATH "../../../../afghan/"
#define XML_NAME "afghan_drop_step"

using namespace std;
using namespace ci;

int createPosFiles(int startFile, int endFile) {
  for (int i=startFile; i<=endFile; i++) {
    cout << "Creating pos file for frame " << i << "...\n";
    
    // Extract data from XML file for given frame and get the file of control point positions
    XmlTree tree;
    stringstream filename;
    filename << DIR_PATH << XML_NAME << i << "/" << XML_NAME << i << ".xml";
    
    try {
      tree = XmlTree(loadFile(filename.str()));
    } catch (exception e) {
      cerr << "Cannot load" << filename.str() << "\n";
      exit(1);
    }
    XmlTree garment = tree.getChild(GARMENT_PATH);
    string posfilename = garment.getChild("position").getAttribute("file");
    posfilename = DIR_PATH + posfilename;
    ifstream posfile;
    
    
    try {
      posfile.open(posfilename);
    } catch (exception e) {
      cout << e.what();
      exit(1);
    }
    
    if (!posfile.is_open()) {
      cerr << posfilename;
      exit(1);
    }
    
    vector<float> newPoints;
    while (!posfile.eof()) {
      string line;
      getline(posfile, line);
      if (!line.empty()) {
        newPoints.push_back(stof(line));
      }
    }
    
    posfile.close();
    
    stringstream ofilename;
    ofilename << POS_PATH << i << ".pos";
    
    ofstream posOutFile(ofilename.str(), ios::binary | ios::trunc);
    if (!posOutFile.is_open()) {
      cerr << "Error: failed to create file: " << ofilename.str() << "\n";
      exit(1);
    }
    for (int k=0; k<newPoints.size(); k++) {
      char* out = (char*)&(newPoints[k]); // Here be more dragons.
      for(int j=0; j<sizeof(float); j++) {
        posOutFile << out[j]; // Even more dragons...
      }
    }
    posOutFile.close();
    
  }

  return 0;
}
