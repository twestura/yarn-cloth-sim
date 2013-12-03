#include <iostream>
#include <fstream>
#include "cinder/app/AppNative.h"
#include "cinder/gl/gl.h"
#include "cinder/Xml.h"
#include "cinder/KdTree.h"
#include "cinder/Camera.h"
#include "cinder/gl/GlslProg.h"
#include "cinder/gl/Vbo.h"
#include "Resources.h"
#include "InfoPanel.h"
#include "Eigen/Dense"

using namespace ci;
using namespace ci::app;
using namespace std;

#define CAMERA_DELTA .3
#define GARMENT_PATH "data/item_afghan_onplane"
#define XML_PATH "../../../../afghan/afghan_drop_step"
#define XML_NAME "afghan_drop_step"
#define NUM_FRAMES 3

#define RESID_PATH "../../../../afghan/resid/"

struct ControlPoint {
  float x;
  float y;
  float z;
  Vec3f frameU;
  Vec3f frameV;
  float rotation;
};

struct NeighborLookupProc {
  vector<uint32_t> neighbors;
  void process( uint32_t id, float distSqrd, float &maxDistSqrd )
  {
    neighbors.push_back(id);
  }
};

class VisualizerApp : public AppNative {
  public:
	void setup();
  void keyDown( KeyEvent event );
  void keyUp( KeyEvent event );
	void mouseDown( MouseEvent event );
  void mouseWheel( MouseEvent event );
	void update();
	void draw();
  void drawInfoPanel();
  
  float getResidual(const NeighborLookupProc) const;
  void loadFrame(const int);
  void viewFrame(const int);
  
  InfoPanel	mInfoPanel;
  float mCounter = 0;
  
  vector<Vec3f> points[NUM_FRAMES];
  vector<ColorA> colors;
  int currentFrame = 0;
  
  gl::GlslProg threadShader;
  
  float radius = 1;
  KdTree<Vec3f, 3, NeighborLookupProc> kdtree;
  gl::VboMeshRef	mVboMesh;
  NeighborLookupProc selection = NeighborLookupProc();
  
  CameraPersp camera = CameraPersp();
  Vec3f eyePos = Vec3f( 100, 40, 0 );
  Vec3f targetPos = Vec3f( -65, 16, 36 );
  unsigned char camBitfield = 0;
  
};

void VisualizerApp::setup()
{
  // Load all frames
  for (int i=0; i<NUM_FRAMES; i++) {
    loadFrame(i);
  }
  
  // Push the points to the GPU
  gl::VboMesh::Layout layout;
  layout.setDynamicPositions();
  layout.setStaticIndices();
  layout.setDynamicColorsRGBA();
  mVboMesh = gl::VboMesh::create(points[0].size(), 2*(points[0].size()-1), layout, GL_LINES);
  
  vector<uint32_t> indexBuffer;
  indexBuffer.reserve(2*(points[0].size()-1));
  for(int i=0; i<2*(points[0].size()-1); i++) {
    indexBuffer.push_back(i);
    indexBuffer.push_back(i+1);
  }
  mVboMesh->bufferIndices(indexBuffer);
  
  viewFrame(0);
  
  kdtree.initialize(points[0]);
  
  camera.setPerspective(40, getWindowAspectRatio(), 1, 1000);
  camera.lookAt( eyePos, targetPos, Vec3f( 0, 1, 0 ) );
  
  mInfoPanel.createTexture();
  
  gl::Texture depthTex;
//  threadShader = gl::GlslProg(loadResource(RES_VERT_GLSL), loadResource(RES_FRAG_GLSL));
}

void VisualizerApp::mouseDown( MouseEvent event )
{
  
  if (event.isRight()) {
    
    stringstream filename;
    filename << RESID_PATH << currentFrame << "-" << radius << ".resid";
    ifstream residFile(filename.str());
    float residual[points[0].size()];
    float minRes = INFINITY;
    float maxRes = -INFINITY;
    
    if (residFile) {
      cout << "Cached file found. \n";
      for (int i=0; i<points[0].size(); i++) {
        string line;
        getline(residFile, line);
        residual[i] = stof(line);
        if (residual[i] > maxRes) {
          maxRes = residual[i];
        }
        if (residual[i] < minRes) {
          minRes = residual[i];
        }
      }
      residFile.close();
    } else {
      // compute all residuals, save to file
      ofstream residOutFile(filename.str());
      if (!residOutFile.is_open()) {
        cerr << "Warning: failed to create cache: " << filename.str() << "\n";
      }
      int percent = 0;
      cout << "No cached file. Computing residuals... \n";
      for (int i=0; i<points[0].size(); i++) {
        int p =(int)((float)i/points[0].size()*100);
        if (p != percent && p % 10 == 0) {
          cout << p << "%\n";
          percent = p;
        }
        NeighborLookupProc nlp = NeighborLookupProc();
        kdtree.lookup(points[0][i], nlp, radius);
        residual[i] = getResidual(nlp);
        residOutFile << residual[i] << "\n";
        if (residual[i] > maxRes) {
          maxRes = residual[i];
        }
        if (residual[i] < minRes) {
          minRes = residual[i];
        }
      }
      residOutFile.close();
      cout << "Done!\n";
    }
    
    // set colors
    colors.clear();
    for (int i=0; i<points[0].size(); i++) {
      float c = (residual[i]-minRes)/(maxRes-minRes);
      colors.push_back(ColorA(c, 0.4, 1-c, 0.4));
    }
    
    viewFrame(currentFrame);
    
  } else {
    Vec2f mouse = event.getPos();
    Vec2f winDem = getWindowSize();
    Ray mouseRay = camera.generateRay(mouse.x/winDem.x , (winDem.y-mouse.y)/winDem.y, getWindowAspectRatio());
    
    int index = -1;
    float minDist = INFINITY;
    for (int i=0; i<points[0].size(); i++) {
      Vec3f toPoint = points[currentFrame][i] - mouseRay.getOrigin();
      float thisDist = (toPoint - (toPoint.dot(mouseRay.getDirection()))*mouseRay.getDirection()).length();
      if (thisDist < minDist) {
        minDist = thisDist;
        index = i;
      }
    }
    
    if (!selection.neighbors.empty() && !colors.empty()) {
      for (int v : selection.neighbors) {
        colors[v] = ColorA(0.4, 0.4, 0.4, 0.4);
      }
    }
    
    selection.neighbors.clear();
    kdtree.lookup(points[currentFrame][index], selection, radius);
    
    if (colors.size() == 0) {
      sort(selection.neighbors.begin(), selection.neighbors.end());
      int i=0;
      for (int v : selection.neighbors) {
        while (i < v) {
          colors.push_back(ColorA(0.4, 0.4, 0.4, 0.4));
          ++i;
        }
        colors.push_back(ColorA(1, 0, 0, 0.6));
        ++i;
      }
      while (i < points[0].size()) {
        colors.push_back(ColorA(0.4, 0.4, 0.4, 0.4));
        ++i;
      }
    } else {
      for (int v : selection.neighbors) {
        colors[v] = ColorA(1, 0, 0, 0.6);
      }
    }
    
    viewFrame(currentFrame);
    targetPos = points[currentFrame][index];
    camera.lookAt(targetPos);
  }

}

void VisualizerApp::keyDown( KeyEvent event )
{
  switch (event.getCode()) {
    case event.KEY_w:
      if (!(camBitfield & 1)) {
        camBitfield += 1;
      }
      break;
    case event.KEY_s:
      if (!(camBitfield & 2)) {
        camBitfield += 2;
      }
      break;
    case event.KEY_a:
      if (!(camBitfield & 4)) {
        camBitfield += 4;
      }
      break;
    case event.KEY_d:
      if (!(camBitfield & 8)) {
        camBitfield += 8;
      }
      break;
    case event.KEY_e:
      if (!(camBitfield & 16)) {
        camBitfield += 16;
      }
      break;
    case event.KEY_q:
      if (!(camBitfield & 32)) {
        camBitfield += 32;
      }
      break;
      
    case event.KEY_LEFT:
      if(currentFrame > 0) {
        currentFrame--;
        viewFrame(currentFrame);
      }
      break;
    case event.KEY_RIGHT:
      if(currentFrame < NUM_FRAMES-1) {
        currentFrame++;
        viewFrame(currentFrame);
      }
      break;
      
    case event.KEY_ESCAPE:
      exit(0);
      break;
      
    default:
      break;
  }
  
}

void VisualizerApp::keyUp( KeyEvent event )
{
  switch (event.getCode()) {
    case event.KEY_w:
      camBitfield -= 1;
      break;
    case event.KEY_s:
      camBitfield -= 2;
      break;
    case event.KEY_a:
      camBitfield -= 4;
      break;
    case event.KEY_d:
      camBitfield -= 8;
      break;
    case event.KEY_e:
      camBitfield -= 16;
      break;
    case event.KEY_q:
      camBitfield -= 32;
      break;
      
    default:
      break;
  }
}

void VisualizerApp::mouseWheel( MouseEvent event )
{
  float inc = event.getWheelIncrement();
  if (event.isShiftDown()) {
    radius += inc / 10;
  } else {
    Vec3f los = camera.getViewDirection();
    eyePos += los * inc;
    camera.lookAt(eyePos, targetPos);
  }
}

void VisualizerApp::update()
{
  Vec3f delta = Vec3f::zero();
  if (camBitfield & 1) delta += Vec3f (1, 0, 0);
  if (camBitfield & 2) delta -= Vec3f (1, 0, 0);
  if (camBitfield & 4) delta += Vec3f (0, 0, 1);
  if (camBitfield & 8) delta -= Vec3f (0, 0, 1);
  if (camBitfield & 16) delta += Vec3f (0, 1, 0);
  if (camBitfield & 32) delta -= Vec3f (0, 1, 0);
  
  eyePos += delta * CAMERA_DELTA;
  camera.lookAt( eyePos, targetPos, Vec3f( 0, 1, 0 ) );
}

void VisualizerApp::draw()
{
  gl::enableAlphaBlending();
	gl::enableDepthRead( true );
	gl::enableDepthWrite( true );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  
  // Set projection/modelview matrices
  gl::setMatrices(camera);
	// clear out the window with black
	gl::clear( Color( 0, 0, 0 ) );
  // draw the scene
  gl::draw(mVboMesh);
  
  // drawInfoPanel();
  
  mCounter++;
}

void VisualizerApp::drawInfoPanel()
{
	glDisable( GL_LIGHTING );
	glEnable( GL_TEXTURE_2D );
	glColor4f( 1, 1, 1, 1 );
	
	gl::pushMatrices();
	gl::setMatricesWindow( getWindowSize() );
	mInfoPanel.update( Vec2f( getWindowWidth(), getWindowHeight() ), 0 );
	gl::popMatrices();
  
	glDisable( GL_TEXTURE_2D );
}


float VisualizerApp::getResidual(const NeighborLookupProc nlp) const {
  using namespace Eigen;
  int n = nlp.neighbors.size();
  Matrix<double, 3, Dynamic> P(3, n);
  Matrix<double, 3, Dynamic> Q(3, n);
  
  Vec3f avg = Vec3f();
  Vec3f avgNext = Vec3f();
  for (int i=0; i<n; i++) {
    avg += points[currentFrame][nlp.neighbors[i]];
    avgNext += points[currentFrame+1][nlp.neighbors[i]];
  }
  
  avg /= n;
  avgNext /= n;
  
  for(int i=0; i<n; i++) {
    P(0,i) = points[currentFrame][nlp.neighbors[i]].x - avg.x;
    P(1,i) = points[currentFrame][nlp.neighbors[i]].y - avg.y;
    P(2,i) = points[currentFrame][nlp.neighbors[i]].z - avg.z;
    Q(0,i) = points[currentFrame+1][nlp.neighbors[i]].x - avgNext.x;
    Q(1,i) = points[currentFrame+1][nlp.neighbors[i]].y - avgNext.y;
    Q(2,i) = points[currentFrame+1][nlp.neighbors[i]].z - avgNext.z;
  }
  
  
  // Q = MP
  // M = QP'[PP']^-1
  Matrix<double, 3, 3> A = P * P.transpose();
  Matrix<double, 3, 3> B = P * Q.transpose();
  Matrix<double, 3, 3> M = (A.fullPivHouseholderQr().solve(B)).transpose();
  return (float) (Q - M*P).norm();
}

void VisualizerApp::loadFrame(const int frame)
{
  assert(frame < NUM_FRAMES && frame >=0);
  // Extract data from XML file for given frame and get the file of control point positions
  XmlTree tree;
  stringstream filename;
  filename << XML_PATH << frame << "/" << XML_NAME << frame << ".xml";
  
  try {
    tree = XmlTree(loadFile(filename.str()));
  } catch (exception e) {
    cerr << "Cannot load" << filename.str() << "\n";
    exit(0);
  }
  XmlTree garment = tree.getChild(GARMENT_PATH);
  string posfilename = garment.getChild("position").getAttribute("file");
  ifstream posfile;
  try {
    posfile.open(posfilename);
  } catch (exception e) {
    cout << e.what();
    exit(1);
  }
  
  while (!posfile.eof()) {
    string line;
    float pos[3];
    bool set = true;
    
    for (int j=0; j<3; j++) {
      assert(!posfile.eof());
      try {
        getline(posfile, line);
        if (line.empty()) {
          set = false;
          break;
        }
        pos[j] = stof(line);
      } catch (exception e) {
        cout << line << " " << j;
        exit(1);
      }
    }
    
    if (set) {
      points[frame].push_back(Vec3f(pos[0], pos[1], pos[2]));
    }
  }
  
  posfile.close();
}

void VisualizerApp::viewFrame(const int frame)
{
  gl::VboMesh::VertexIter vIter = mVboMesh->mapVertexBuffer();
  for (int i=0; i<points[0].size(); i++) {
    vIter.setPosition(points[frame][i]);
    vIter.setColorRGBA(colors.size() == 0 ? ColorA(0.4, 0.4, 0.4, 0.4) : colors[i]);
    ++vIter;
  }
}

CINDER_APP_NATIVE( VisualizerApp, RendererGl )
