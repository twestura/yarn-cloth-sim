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
  
  InfoPanel	mInfoPanel;
  float mCounter = 0;
  
  float radius = 1;
  vector<Vec3f> points;
  vector<Vec3f> pointsNext;
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
  // Extract data from XML file for given frame and get the file of control point positions
  XmlTree tree;
  try {
    tree = XmlTree(loadFile(getOpenFilePath()));
  } catch (exception e) {
    exit(0);
  }
  XmlTree garment = tree.getChild(GARMENT_PATH);
  string posfilename = garment.getChild("position").getAttribute("file");
  if (posfilename.empty()) exit(0);
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
      points.push_back(Vec3f(pos[0], pos[1], pos[2]));
    }
  }
  
  posfile.close();
  
  // Push the points to the GPU
  gl::VboMesh::Layout layout;
  layout.setStaticPositions();
  layout.setStaticIndices();
  layout.setDynamicColorsRGBA();
  mVboMesh = gl::VboMesh::create(points.size(), 2*(points.size()-1), layout, GL_LINES);
  
  vector<uint32_t> indexBuffer;
  indexBuffer.reserve(2*(points.size()-1));
  for(int i=0; i<2*(points.size()-1); i++) {
    indexBuffer.push_back(i);
    indexBuffer.push_back(i+1);
  }
  mVboMesh->bufferPositions(points);
  mVboMesh->bufferIndices(indexBuffer);
  
  gl::VboMesh::VertexIter vIter = mVboMesh->mapVertexBuffer();
  while (vIter.isDone()) {
    vIter.setColorRGBA(ColorA(0.4, 0.4, 0.4, 0.4));
    ++vIter;
  }
  
  kdtree.initialize(points);
  
  camera.setPerspective(40, getWindowAspectRatio(), 1, 1000);
  camera.lookAt( eyePos, targetPos, Vec3f( 0, 1, 0 ) );
  
  mInfoPanel.createTexture();
}

void VisualizerApp::mouseDown( MouseEvent event )
{
  
  if (event.isRight()) {
    pointsNext.clear();
    pointsNext.reserve(points.size());
  
    XmlTree tree;
    try {
      tree = XmlTree(loadFile(getOpenFilePath()));
    } catch (exception e) {
      cout << e.what();
      return;
    }
    XmlTree garment = tree.getChild(GARMENT_PATH);
    string posfilename = garment.getChild("position").getAttribute("file");
    if (posfilename.empty()) return;
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
        pointsNext.push_back(Vec3f(pos[0], pos[1], pos[2]));
      }
    }
    
    posfile.close();
    assert(points.size() == pointsNext.size());
    cout << getResidual(selection) << "\n";
    
  } else {
    Vec2f mouse = event.getPos();
    Vec2f winDem = getWindowSize();
    Ray mouseRay = camera.generateRay(mouse.x/winDem.x , (winDem.y-mouse.y)/winDem.y, getWindowAspectRatio());
    
    int index = -1;
    float minDist = INFINITY;
    for (int i=0; i<points.size(); i++) {
      Vec3f toPoint = points[i] - mouseRay.getOrigin();
      float thisDist = (toPoint - (toPoint.dot(mouseRay.getDirection()))*mouseRay.getDirection()).length();
      if (thisDist < minDist) {
        minDist = thisDist;
        index = i;
      }
    }
    selection.neighbors.clear();
    kdtree.lookup(points[index], selection, radius);
    sort(selection.neighbors.begin(), selection.neighbors.end());
    
    gl::VboMesh::VertexIter vIter = mVboMesh->mapVertexBuffer();
    int i=0;
    for (int v : selection.neighbors) {
      while (i < v) {
        vIter.setColorRGBA(ColorA(0.4, 0.4, 0.4, 0.4));
        ++vIter;
        ++i;
      }
      vIter.setColorRGBA(ColorA(1, 0, 0, 0.6));
      ++vIter;
      ++i;
    }
    while (vIter.isDone()) {
      vIter.setColorRGBA(ColorA(0.4, 0.4, 0.4, 0.4));
      ++vIter;
    }
    targetPos = points[index];
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
  Matrix<float, 4, Dynamic> P(4, n);
  Matrix<float, 3, Dynamic> Q(3, n);
  
  for(int i=0; i<n; i++) {
    P(0,i) = points[nlp.neighbors[i]].x;
    P(1,i) = points[nlp.neighbors[i]].y;
    P(2,i) = points[nlp.neighbors[i]].z;
    P(3,i) = 1;
    Q(0,i) = pointsNext[nlp.neighbors[i]].x;
    Q(1,i) = pointsNext[nlp.neighbors[i]].y;
    Q(2,i) = pointsNext[nlp.neighbors[i]].z;
  }
  // Q = MP
  // M = QP'[PP']^-1
  Matrix<float, 4, 4> A = P * P.transpose();
  Matrix<float, 4, 3> B = P * Q.transpose();
  Matrix<float, 3, 4> M = (A.fullPivHouseholderQr().solve(B)).transpose();
  return (Q - M*P).norm();
}

CINDER_APP_NATIVE( VisualizerApp, RendererGl )
