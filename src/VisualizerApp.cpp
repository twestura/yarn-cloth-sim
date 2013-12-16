#include <iostream>
#include <fstream>
#include <thread>
#include "cinder/app/AppNative.h"
#include "cinder/gl/gl.h"
#include "cinder/Xml.h"
#include "cinder/KdTree.h"
#include "cinder/Camera.h"
#include "cinder/gl/GlslProg.h"
#include "cinder/gl/Vbo.h"
#include "Resources.h"
#include "InfoPanel.h"
#include "PosFileConv.h"
#include "Eigen/Dense"
#include <boost/circular_buffer.hpp>

using namespace ci;
using namespace ci::app;
using namespace std;

#define CAMERA_DELTA .3
#define GARMENT_PATH "data/item_afghan_onplane"
#define DIR_PATH "../../../../afghan/"
#define XML_NAME "afghan_drop_step"
#define NUM_FRAMES 7
#define START_FRAME 0
#define END_FRAME 60

#define RESID_PATH "../../../../afghan/resid/"
#define POS_PATH "../../../../afghan/pos/"
#define MIN_RES 0

//#define CREATE_POS_FILES
//#define CREATE_RESID_FILES
//#define RADIUS_ONE_HALF

#ifdef RADIUS_ONE_HALF
#define RADIUS 0.5
#define MAX_RES 0.08 // TODO: This should not be hard-coded.
#else
#define RADIUS 1
#define MAX_RES 2.07
#endif // ifdef RADIUS_ONE_HALF

struct Frame {
  float minResid;
  float maxResid;
  vector<float> resid;
  vector<Vec3f> pos;
  
  Frame () {}
  
  Frame (vector<Vec3f> p) {
    pos = p;
    resid.clear();
    minResid = maxResid = 0;
  }
  
  Frame (vector<Vec3f> p, vector<float> r, float min, float max) {
    pos = p;
    resid = r;
    minResid = min;
    maxResid = max;
  }
  
  Vec3f& operator[]( const uint32_t index )
  {
    return pos[index];
  }
  
  const Vec3f& operator[]( const uint32_t index ) const
  {
    return pos[index];
  }
  
  uint32_t size()
  {
    return pos.size();
  }
};

struct NeighborLookupProc {
  vector<uint32_t> neighbors;
  void process( uint32_t id, float distSqrd, float &maxDistSqrd )
  {
    neighbors.push_back(id);
  }
};

template <class T>
struct ModCircularBuffer {
private:
  bool is_empty = true;
  int offset = 0;
  T zero;
  boost::circular_buffer<T> cb;
public:
  
  void init( const uint32_t size, int start_frame )
  {
    cb = boost::circular_buffer<T>(size);
    offset = (start_frame == 0 ? 0 : -start_frame+1);
  }
  
  T& operator[]( const uint32_t index )
  {
    if (index == 0) {
      return zero;
    } else {
      return cb[(index-1) + offset];
    }
  }
  
  const T& operator[]( const uint32_t index ) const
  {
    if (index == 0) {
      return zero;
    } else {
      return cb[(index-1) + offset];
    }
  }
  
  void push_back( T elt )
  {
    if (is_empty) {
      zero = elt;
      is_empty = false;
    } else {
      if (cb.full()) offset--;
      cb.push_back(elt);
    }
  }
  
  void push_front( T elt )
  {
    if (cb.full()) offset++;
    cb.push_front(elt);
  }
  
  int left_buffer_size( int frame )
  {
    return cb.capacity() - frame - offset;
  }
  
  int right_buffer_size( int frame )
  {
    return frame-1 + offset;
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
  void readHandlerBack( const boost::system::error_code& ec,
                       std::size_t bytes_transferred );
  void readHandlerFront( const boost::system::error_code& ec,
                        std::size_t bytes_transferred );
  
  float getResidual(const NeighborLookupProc) const;
  void loadFrame(const int, const bool);
  void viewFrame(const int);
  void writeResidFile();
  void createResidFiles(const int);
  
  InfoPanel	mInfoPanel; // TODO: remove this?
  float mCounter = 0;
  
//  mutable_buffers_1 mbs;
  
  ModCircularBuffer<Frame> points;
  int currentFrame = START_FRAME;
  bool globalRelativeColoring = false;
  
  gl::GlslProg threadShader;
  
  float radius = RADIUS;
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
#ifdef CREATE_POS_FILES
  createPosFiles(START_FRAME, END_FRAME);
  exit(0);
#endif
  
  // Load starting frames
  points.init(NUM_FRAMES, START_FRAME);
  loadFrame(0, true);
  
  int start =START_FRAME == 0 ? 1 : START_FRAME;
  for (int i=start; i<start+NUM_FRAMES; i++) {
    loadFrame(i, true);
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
  
  viewFrame(currentFrame);
  
  kdtree.initialize(points[0].pos);
  
  camera.setPerspective(40, getWindowAspectRatio(), 1, 1000);
  camera.lookAt( eyePos, targetPos, Vec3f( 0, 1, 0 ) );
  
  mInfoPanel.createTexture();
  
#ifdef CREATE_RESID_FILES
  createResidFiles(END_FRAME-1);
  exit(0);
#endif
  
  // TODO: remove this
//  gl::Texture depthTex;
//  threadShader = gl::GlslProg(loadResource(RES_VERT_GLSL), loadResource(RES_FRAG_GLSL));
}

void VisualizerApp::mouseDown( MouseEvent event )
{
  
  if (event.isRight()) {
    stringstream filename;
    filename << RESID_PATH << currentFrame << "-" << radius << ".resid";
    ifstream residFile(filename.str(), ios::binary);
    
    if (!residFile) {
      writeResidFile();
    }
    
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
      if(currentFrame > START_FRAME) {
        currentFrame--;
        viewFrame(currentFrame);
        if (points.right_buffer_size(currentFrame) < 1 && START_FRAME != currentFrame && 1 != currentFrame) {
          for (int i=1; i<NUM_FRAMES-2 && currentFrame - i >= START_FRAME && currentFrame - i >= 1; i++) {
            loadFrame(currentFrame - i, false);
          }
        }
      }
      break;
    case event.KEY_RIGHT:
      if(currentFrame < END_FRAME) {
        currentFrame++;
        viewFrame(currentFrame);
        if (points.left_buffer_size(currentFrame) < 1 && END_FRAME != currentFrame) {
          for (int i=1; i<NUM_FRAMES-2 && currentFrame + i <= END_FRAME; i++) {
            loadFrame(currentFrame + i, true);
          }
        }
      }
      break;
      
    case event.KEY_m:
      globalRelativeColoring = !globalRelativeColoring;
      viewFrame(currentFrame);
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


float VisualizerApp::getResidual(const NeighborLookupProc nlp) const
{
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

void VisualizerApp::loadFrame(const int frame, const bool back)
{
 // cout << "load frame " << frame << " to " << (back ? "back" : "front") << "\n";
  
  stringstream posFilename;
  posFilename << POS_PATH << frame << ".pos";
  ifstream posFile(posFilename.str(), ios::binary);
  
  if (!posFile) {
    cerr << "Error: position file not found: " << posFilename.str() << "\n";
    exit(1);
  }
  
  stringstream residFilename;
  residFilename << RESID_PATH << frame << "-" << radius << ".resid";
  ifstream residFile(residFilename.str(), ios::binary);
  
  if (!residFile) {
    cerr << "Warning: resid file not found: " << residFilename.str() << "\n";
  }
  
  vector<Vec3f> newPoints;
  vector<float> newResid;
  float minResid = INFINITY;
  float maxResid = -INFINITY;
  
  while(!posFile.eof()) {
    float point[4];
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
    points.push_back(f);
  } else {
    points.push_front(f);
  }
  
  posFile.close();
  
  /*
  if (back) {
    async_read(posfile, mbs, boost::bind(&VisualizerApp::readHandlerBack,
                                         this, boost::asio::placeholders::error,
                                         boost::asio::placeholders::bytes_transferred));
  } else {
    async_read(posfile, mbs, boost::bind(&VisualizerApp::readHandlerFront,
                                         this, boost::asio::placeholders::error,
                                         boost::asio::placeholders::bytes_transferred));
  }
  */
}

/*

void VisualizerApp::readHandlerBack(const boost::system::error_code& ec,
                 size_t bytes_transferred)
{
  if (ec) {
    cout << "Error reading pos file: " << ec << "\n";
    exit(1);
  }
  vector<Vec3f> newPoints;
  // TODO: reserve space
  int i=0;
  for (i=0; i<bytes_transferred; i++) {
    bool set = true;
    float pos[3];
    for (int j=0; j<3; j++) {
      stringstream s;
      while (mbs[i] != '\n') {
        s << mbs[i];
        i++;
      }
      string line = s.str();
      if (line.empty()) {
        set = false;
        break;
      }
      pos[j] = stof(line);
    }
    if (set) {
      newPoints.push_back(Vec3f(pos[0], pos[1], pos[2]));
    }
  }
  
  points.push_back(newPoints);
}
 
*/

void VisualizerApp::viewFrame(const int frame)
{
  float minResid = (globalRelativeColoring ? MIN_RES : points[frame].minResid);
  float maxResid = (globalRelativeColoring ? MAX_RES : points[frame].maxResid);
  gl::VboMesh::VertexIter vIter = mVboMesh->mapVertexBuffer();
  for (int i=0; i<points[0].size(); i++) {
    vIter.setPosition(points[frame][i]);
    float c = (points[frame].resid.empty() ? 0 : (points[frame].resid[i]-minResid)/(maxResid-minResid));
    vIter.setColorRGBA(ColorA(c, 0.4, 1-c, 0.6));
    ++vIter;
  }
}

void VisualizerApp::writeResidFile()
{
  stringstream filename;
  filename << RESID_PATH << currentFrame << "-" << radius << ".resid";
  
  cout << "Creating cache for frame " << currentFrame << "...\n";
  // compute all residuals, save to file
  ofstream residOutFile(filename.str(), ios::binary | ios::trunc);
  if (!residOutFile.is_open()) {
    cerr << "Warning: failed to create cache: " << filename.str() << "\n";
  }
  int percent = 0;
  points[currentFrame].minResid = INFINITY;
  points[currentFrame].maxResid = -INFINITY;
  points[currentFrame].resid.clear();
  for (int i=0; i<points[0].size(); i++) {
    int p =(int)((float)i/points[0].size()*100);
    if (p != percent && p % 10 == 0) {
      cout << p << "%\n";
      percent = p;
    }
    NeighborLookupProc nlp = NeighborLookupProc();
    kdtree.lookup(points[0][i], nlp, radius);
    points[currentFrame].resid.push_back(getResidual(nlp));
    
    if (points[currentFrame].resid[i] > points[currentFrame].maxResid) {
      points[currentFrame].maxResid = points[currentFrame].resid[i];
    }
    if (points[currentFrame].resid[i] < points[currentFrame].minResid) {
      points[currentFrame].minResid = points[currentFrame].resid[i];
    }
    
    char* out = (char*)&(points[currentFrame].resid[i]); // Here be more dragons.
    for(int j=0; j<sizeof(float); j++) {
      residOutFile << out[j]; // Even more dragons...
    }
  }
  residOutFile.close();
  cout << "Done!\n";
  
}

void VisualizerApp::createResidFiles(const int endFile)
{
  while (currentFrame <= endFile) {
    writeResidFile();
    currentFrame++;
    if (points.left_buffer_size(currentFrame) < 1) {
      for (int i=1; i<NUM_FRAMES-2 && currentFrame + i <= END_FRAME; i++) {
        loadFrame(currentFrame + i, true);
      }
    }
  }
}

CINDER_APP_NATIVE( VisualizerApp, RendererGl )
