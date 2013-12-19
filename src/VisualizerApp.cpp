#include <iostream>
#include <fstream>
#include "cinder/app/AppNative.h"
#include "cinder/gl/gl.h"
#include "cinder/Xml.h"
#include "cinder/KdTree.h"
#include "cinder/Camera.h"
#include "cinder/gl/Vbo.h"
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
#define COMP_PATH "../../../../afghan/comp/"

#define RESOLUTION (1E-7)
#define MAX_STEPS 32767
#define THRESHOLD (RESOLUTION*MAX_STEPS)

//#define CREATE_POS_FILES
//#define CREATE_RESID_FILES
//#define COMPRESS
//#define DECOMPRESS
//#define RADIUS_ONE_HALF

#define MIN_RES 0
#ifdef RADIUS_ONE_HALF
#define RADIUS 0.5
#define MAX_RES 0.08 // TODO: This should not be hard-coded.
#else
#define RADIUS 1
#define MAX_RES 2.07
#endif // ifdef RADIUS_ONE_HALF

// A struct representing one frame of the animation.
// Optionally stores residual information.
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
  
  const uint32_t size() const
  {
    return pos.size();
  }
};

// A struct for use in the KD-tree. Finds all control points within
// a given squared distance.
struct NeighborLookupProc {
  vector<uint32_t> neighbors;
  void process( uint32_t id, float distSqrd, float &maxDistSqrd )
  {
    neighbors.push_back(id);
  }
};

// A group of control points within an AABB.
// TODO: This is not a bounding volume heirchy. Rename this to something
// more accurate.
struct Bvh {
  ci::Vec3f aabb[2];
  std::vector<uint32_t> indices;
  bool ok = false;
};

// A modified circular buffer to load in frames dynamically. Always allows
// access to the first element of the buffer.
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
  
  float getResidual(const NeighborLookupProc, const bool) const;
  float newGetResidual(const vector<uint32_t>, const bool) const;
  void loadFrame(const int, const bool);
  void viewFrame(const int);
  void writeResidFile();
  void createResidFiles(const int);
  
  void compress();
  void divideBvh(const int index);
  bool getBoundingBox(Bvh* bvh);
  int writeKeyframeFile() const;
  int writeCompressedFrame();
  void decompressFrame(const int frame);
  
  ModCircularBuffer<Frame> points;
  int currentFrame = START_FRAME;
  bool globalRelativeColoring = false;
  
  float radius = RADIUS;
  KdTree<Vec3f, 3, NeighborLookupProc> kdtree;
  gl::VboMeshRef	mVboMesh;
  NeighborLookupProc selection = NeighborLookupProc();
  
  vector<Bvh> bvhVec;
  
  CameraPersp camera = CameraPersp();
  Vec3f eyePos = Vec3f( 100, 40, 0 );
  Vec3f targetPos = Vec3f( -65, 16, 36 );
  unsigned char camBitfield = 0;
};

void VisualizerApp::setup()
{
  
#ifdef CREATE_POS_FILES
  // create binary position files
  createPosFiles(START_FRAME, END_FRAME);
  exit(0);
#endif
  
#ifdef DECOMPRESS
  // decompress and load compressed frames
  points.init(NUM_FRAMES, START_FRAME);
  decompressFrame(0);
  decompressFrame(1);
//  decompressFrame(2);
//  decompressFrame(3);
//  decompressFrame(4);
//  decompressFrame(5);
#else
  // Load starting frames
  points.init(NUM_FRAMES, START_FRAME);
  loadFrame(0, true);
  
  int start =START_FRAME == 0 ? 1 : START_FRAME;
  for (int i=start; i<start+NUM_FRAMES; i++) {
    loadFrame(i, true);
  }
#endif // ifdef DECOMPRESS
  
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
  
#ifdef CREATE_RESID_FILES
  createResidFiles(END_FRAME-1);
  exit(0);
#endif
  
#ifdef COMPRESS
  compress();
#endif
}

// Triggered when the mouse is clicked.
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
    
    selection.neighbors.clear();
    targetPos = points[currentFrame][index];
    kdtree.lookup(targetPos, selection, radius);
    camera.lookAt(targetPos);
//    float n = getResidual(selection, false);
//    cout << "old resid: " << n << "\n";
  }

}

// Triggered when a key is pressed.
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

// Triggered when a key is released.
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

// Triggered on a mouse wheel event.
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

// Called when the application is running.
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

// Draw the scene.
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
  
  gl::color(0, 1, 0, 0.6);
  /*
  for (Bvh bvh : bvhVec) {
    gl::drawStrokedCube(AxisAlignedBox3f(bvh.aabb[0], bvh.aabb[1]));
  }
   */
}

// Given a group of points, compute the residual from the least squares
// best fit transformation from the current frame to the next. If retMax
// is true, the maximum error is returned; otherwise, the total error squared
// is returned.
float VisualizerApp::getResidual(const NeighborLookupProc nlp, bool retMax) const
{
  using namespace Eigen;
  int n = nlp.neighbors.size();
  
  if (n==1) return 0;
  
  Vec3f avg = Vec3f();
  Vec3f avgNext = Vec3f();
  for (int i=0; i<n; i++) {
    avg += points[currentFrame][nlp.neighbors[i]];
    avgNext += points[currentFrame+1][nlp.neighbors[i]];
  }
  
  avg /= n;
  avgNext /= n;
  
  Matrix<double, 3, 3> A = Array33d::Zero();
  Matrix<double, 3, 3> B = Array33d::Zero();
  for (uint32_t index : nlp.neighbors) {
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
        Vec3f p = points[currentFrame][index] - avg;
        A(i,j) += p[i]*p[j];
        Vec3f pNext = points[currentFrame+1][index] - avgNext;
        B(i,j) += p[i]*pNext[j];
      }
    }
  }
  
  // Q = MP
  // M = QP'[PP']^-1
  Matrix<double, 3, 3> M = (A.fullPivHouseholderQr().solve(B)).transpose();
  
  // Compute residuals
  float ret = 0;
  Vec3f temp;
  Vector3d input;
  Vector3d output;
  for (int index : nlp.neighbors) {
    temp = points[currentFrame][index] - avg;
    input << temp.x, temp.y, temp.z;
    output = M * input;
    for (int j=0; j<3; j++) {
      float diff = abs(output(j) - (points[currentFrame+1][index][j] - avgNext[j]));
      if (retMax) {
        ret = max(ret, diff);
      } else {
        ret += diff * diff;
      }
    }
  }
  
  return ret;
  
}

// Similar to the method above, in an attempt to get things to be more accurate.
// Not entirely successful.
float VisualizerApp::newGetResidual(const vector<uint32_t> indices, const bool retMax) const
{
  using namespace Eigen;
  int n = indices.size();
  
  // Compute centroids
  Vec3f c = Vec3f();
  Vec3f cTilde = Vec3f();
  for (int index : indices) {
    c += points[currentFrame][index];
    cTilde += points[currentFrame][index]; //here
  }
  c /= n;
  cTilde /= n;
  
  // Create least squares system--see James and Twigg, Appendix A.
  Matrix<float, 12, 12> M;
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
  
  Matrix<float, 12, 1> B;
  B << c[0]*cTilde[0], c[1]*cTilde[0], c[2]*cTilde[0], c[0]*cTilde[1], c[1]*cTilde[1], c[2]*cTilde[1],
  c[0]*cTilde[2], c[1]*cTilde[2], c[2]*cTilde[2], c[0], c[1], c[2];
  
  // Solve least squares, convert to affine transform matrix
  Matrix<float, 12, 1> X = M.fullPivLu().solve(B);
  
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
    temp = points[currentFrame][index] - c;
    input << temp.x, temp.y, temp.z, 1;
    output = Transform * input;
    for (int j=0; j<3; j++) {
      float diff = abs(output(j) + cTilde[j] - points[currentFrame][index][j]); //here
      if (retMax) {
        ret = max(ret, diff);
      } else {
        ret += diff * diff;
      }
    }
  }
  
  return ret;
}

// Decompress a given frame (.comp or .key) and load it into
// the circular buffer.
void VisualizerApp::decompressFrame(const int frame)
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
    points.push_back(f);
    
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
        avgPrev += points[frame-1][index];
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
        Vec3f newPoint = points[frame-1][index] - avgPrev;
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
    points.push_back(f);
    
    posFile.close();
  } else {
    cerr << "Error: no compressed frame found: " << posFilename.str() << "\n";
  }
}

// Load a .pos file into the circular buffer.
void VisualizerApp::loadFrame(const int frame, const bool back)
{
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
    points.push_back(f);
  } else {
    points.push_front(f);
  }
  
  posFile.close();
  
}

// Prime a frame to be drawn.
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

// Write a .resid file containing residual information for this frame.
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
    points[currentFrame].resid.push_back(getResidual(nlp, false));
    
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

// Create .resid files for a range of frames.
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

// Compress a range of frames. Create a .comp file if possible; otherwise create a .key file.
void VisualizerApp::compress()
{
  while (currentFrame <= END_FRAME) {
    if (points.left_buffer_size(currentFrame) < 1) {
      for (int i=1; i<NUM_FRAMES-2 && currentFrame + i <= END_FRAME; i++) {
        loadFrame(currentFrame + i, true);
      }
    }
    int size;
    if (currentFrame != START_FRAME) {
      cout << "Compressing frame " << currentFrame << "...\n";
      size = writeCompressedFrame();
    }
    if (currentFrame == START_FRAME || size > points[0].size()*3*sizeof(float)) {
      cout << "Compression failed. Making keyframe...";
      stringstream s;
      s << COMP_PATH << currentFrame << ".comp";
      remove(s.str().c_str());
      
      bvhVec.clear();
      bvhVec.push_back(Bvh());
      for (int i=0; i<points[0].size(); i++) {
        bvhVec[0].indices.push_back(i);
      }
      getBoundingBox(&(bvhVec[0]));
      bool done = false;
      while (!done) {
        done = true;
        for (int i=0; i<bvhVec.size(); i++) {
          if (!bvhVec[i].ok) {
            done = false;
            divideBvh(i);
          }
        }
      }
      cout << "num nodes: " << bvhVec.size() << "\n";
      writeKeyframeFile();
    }
    currentFrame++;
  }
  
  
}

// Set the bounding box of a bvh.
bool VisualizerApp::getBoundingBox(Bvh* bvh)
{
  if (bvh->indices.empty()) return false;
  bvh->aabb[0] = Vec3f(points[currentFrame][bvh->indices[0]]); //min
  bvh->aabb[1] = Vec3f(points[currentFrame][bvh->indices[0]]); //max
  for (int i=1; i<bvh->indices.size(); i++) {
    for (int j=0; j<3; j++) {
      if (points[currentFrame][bvh->indices[i]][j] < bvh->aabb[0][j]) {
        bvh->aabb[0][j] = points[currentFrame][bvh->indices[i]][j];
      }
      if (points[currentFrame][bvh->indices[i]][j] > bvh->aabb[1][j]) {
        bvh->aabb[1][j] = points[currentFrame][bvh->indices[i]][j];
      }
    }
  }
  
  return true;
}

// Divide a bvh recursively.
void VisualizerApp::divideBvh(const int index)
{
  NeighborLookupProc nlp;
  nlp.neighbors = bvhVec[index].indices;
  float resid = getResidual(nlp, true);
  
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
    if (points[currentFrame][index][dim] < split) {
      newIndices.push_back(index);
    } else {
      bvhVec[sibling].indices.push_back(index);
    }
  }
  bvhVec[index].indices = newIndices;
  
  getBoundingBox(&bvhVec[index]);
  getBoundingBox(&bvhVec[sibling]);
  
}

// Write a .key file.
int VisualizerApp::writeKeyframeFile() const
{
  if (bvhVec.empty()) {
    cerr << "Tried to write keyframe with empty bvh\n";
    exit(1);
  }
  
  stringstream outFileName;
  outFileName << COMP_PATH << currentFrame << ".key";
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
        char* out = (char*)&(points[currentFrame][index][i]); // Here be dragons.
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
int VisualizerApp::writeCompressedFrame()
{
  using namespace Eigen;
  
  stringstream outFileName;
  outFileName << COMP_PATH << currentFrame << ".comp";
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
      avgPrev += points[currentFrame-1][index];
      avg += points[currentFrame][index];
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
            Vec3f pPrev = points[currentFrame-1][index] - avgPrev;
            A(i,j) += pPrev[i]*pPrev[j];
            Vec3f p = points[currentFrame][index] - avg;
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
        p(j) = points[currentFrame-1][bvh.indices[i]][j] - avgPrev[j];
      }
      q = Mf * p;
      for (int j=0; j<3; j++) {
        float qNew = q[j] + avg[j];
        float actual = points[currentFrame][bvh.indices[i]][j];
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

CINDER_APP_NATIVE( VisualizerApp, RendererGl )
