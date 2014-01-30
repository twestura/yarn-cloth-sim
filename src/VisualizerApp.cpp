#include <iostream>
#include <fstream>

#include "cinder/app/AppNative.h"
#include "cinder/gl/gl.h"
#include "cinder/Xml.h"
#include "cinder/KdTree.h"
#include "cinder/Camera.h"
#include "cinder/gl/Vbo.h"
#include "cinder/Rand.h"

#include "Common.h"
#include "BetterComp.h"
#include "PosFileConv.h"
#include "Decompressor.h"

using namespace ci;
using namespace ci::app;
using namespace std;

#define CAMERA_DELTA .3
#define GARMENT_PATH "data/item_afghan_onplane"
#define DIR_PATH "../../../../afghan/"
#define XML_NAME "afghan_drop_step"

class VisualizerApp : public AppNative {
  public:
	void setup();
  void keyDown( KeyEvent event );
  void keyUp( KeyEvent event );
	void mouseDown( MouseEvent event );
  void mouseWheel( MouseEvent event );
	void update();
	void draw();
  
  void perturb(const float p);
  void viewFrame(const int);
  void writeResidFile();
  void createResidFiles(const int);
  
  AppData ad;
  Decompressor decomp;
  
  bool globalRelativeColoring = false;
  bool compressedMode = false;
  
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
  // create binary position files
  createPosFiles(START_FRAME, END_FRAME);
  exit(0);
#endif
  
  // Load starting frames
  ad.frames.init(NUM_FRAMES, START_FRAME);
  loadFrame(ad, 0, true);
  
  int start = START_FRAME == 0 ? 1 : START_FRAME;
  for (int i=start; i<start+NUM_FRAMES; i++) {
    loadFrame(ad, i, true);
  }
  
  // Push the points to the GPU
  gl::VboMesh::Layout layout;
  layout.setDynamicPositions();
  layout.setStaticIndices();
  layout.setDynamicColorsRGBA();
  mVboMesh = gl::VboMesh::create(ad.frames[0].size(), 2*(ad.frames[0].size()-1), layout, GL_LINES);
  
  vector<uint32_t> indexBuffer;
  indexBuffer.reserve(2*(ad.frames[0].size()-1));
  for(int i=0; i<2*(ad.frames[0].size()-1); i++) {
    indexBuffer.push_back(i);
    indexBuffer.push_back(i+1);
  }
  mVboMesh->bufferIndices(indexBuffer);
  
  viewFrame(ad.currentFrame);
  
  kdtree.initialize(ad.frames[0].pos);
  
  camera.setPerspective(40, getWindowAspectRatio(), 1, 1000);
  camera.lookAt( eyePos, targetPos, Vec3f( 0, 1, 0 ) );
  
#ifdef CREATE_RESID_FILES
  createResidFiles(END_FRAME-1);
  exit(0);
#endif
  
  cout << "# vertices: " << ad.frames[0].size() <<"\n";

#ifdef COMPRESS
  Compressor c(ad, kdtree);
  c.compress(10);
  exit(0);
#endif
  
}

// Triggered when the mouse is clicked.
void VisualizerApp::mouseDown( MouseEvent event )
{
  
  if (event.isRight()) {
    stringstream filename;
    filename << RESID_PATH << ad.currentFrame << "-" << radius << ".resid";
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
    for (int i=0; i<ad.frames[0].size(); i++) {
      Vec3f toPoint = ad.frames[ad.currentFrame][i] - mouseRay.getOrigin();
      float thisDist = (toPoint - (toPoint.dot(mouseRay.getDirection()))*mouseRay.getDirection()).length();
      if (thisDist < minDist) {
        minDist = thisDist;
        index = i;
      }
    }
    
    selection.neighbors.clear();
    targetPos = ad.frames[ad.currentFrame][index];
    kdtree.lookup(targetPos, selection, radius);
    camera.lookAt(targetPos);
    float n = getResidual(ad, selection.neighbors, ad.currentFrame, ad.currentFrame+1, false);
    float newn = newGetResidual(ad, selection.neighbors, ad.currentFrame, ad.currentFrame+1, false);
    cout << "old: " << n << "\nnew: " << newn << "\n";
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
      if (compressedMode) {
        decomp.currentFrame--;
      } else if(ad.currentFrame > START_FRAME) {
        ad.currentFrame--;
        viewFrame(ad.currentFrame);
        loadFramesIfNecessary(ad, Direction::Left, ad.currentFrame);
      }
      break;
    case event.KEY_RIGHT:
      if (compressedMode) {
        decomp.currentFrame++;
      } else if(ad.currentFrame < END_FRAME) {
        ad.currentFrame++;
        viewFrame(ad.currentFrame);
        loadFramesIfNecessary(ad, Direction::Right, ad.currentFrame);
      }
      break;
      
    case event.KEY_m:
      globalRelativeColoring = !globalRelativeColoring;
      viewFrame(ad.currentFrame);
      break;
      
    case event.KEY_p:
      perturb(.005);
      break;
      
    case event.KEY_c:
      if (compressedMode) {
        compressedMode = false;
        decomp.clear();
      } else {
        compressedMode = true;
        mVboMesh->unbindBuffers();
        decomp.init();
        mVboMesh->bindAllData();
      }
      break;
      
    case event.KEY_ESCAPE:
      exit(0);
      break;
      
    default:
      break;
  }
  
}

void VisualizerApp::perturb(const float p)
{
  float minResid = (globalRelativeColoring ? MIN_RES : ad.frames[ad.currentFrame].minResid);
  float maxResid = (globalRelativeColoring ? MAX_RES : ad.frames[ad.currentFrame].maxResid);
  gl::VboMesh::VertexIter vIter = mVboMesh->mapVertexBuffer();
  Rand r(getElapsedFrames());
  for (int i=0; i<ad.frames[0].size(); i++) {
    Vec3f temp = ad.frames[ad.currentFrame][i];
    for (int j=0; j<3; j++)
      temp[j] += r.nextFloat(-p, p);
    vIter.setPosition(temp);
    float c = (ad.frames[ad.currentFrame].resid.empty() ? 0 : (ad.frames[ad.currentFrame].resid[i]-minResid)/(maxResid-minResid));
    vIter.setColorRGBA(ColorA(c, 0.4, 1-c, 0.6));
    ++vIter;
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
  if (compressedMode) {
    decomp.draw();
  } else {
    mVboMesh->bindAllData();
    gl::draw(mVboMesh);
    mVboMesh->unbindBuffers();
  }
  
  gl::color(0, 1, 0, 0.6);

}

// Prime a frame to be drawn.
void VisualizerApp::viewFrame(const int frame)
{
  float minResid = (globalRelativeColoring ? MIN_RES : ad.frames[frame].minResid);
  float maxResid = (globalRelativeColoring ? MAX_RES : ad.frames[frame].maxResid);
  gl::VboMesh::VertexIter vIter = mVboMesh->mapVertexBuffer();
  for (int i=0; i<ad.frames[0].size(); i++) {
    vIter.setPosition(ad.frames[frame][i]);
    float c = (ad.frames[frame].resid.empty() ? 0 : (ad.frames[frame].resid[i]-minResid)/(maxResid-minResid));
    vIter.setColorRGBA(ColorA(c, 0.4, 1-c, 0.6));
    ++vIter;
  }
}

// Write a .resid file containing residual information for this frame.
void VisualizerApp::writeResidFile()
{
  stringstream filename;
  filename << RESID_PATH << ad.currentFrame << "-" << radius << ".resid";
  
  cout << "Creating resid cache for frame " << ad.currentFrame << ", radius " << radius << "...\n";
  // compute all residuals, save to file
  ofstream residOutFile(filename.str(), ios::binary | ios::trunc);
  if (!residOutFile.is_open()) {
    cerr << "Warning: failed to create resid cache: " << filename.str() << "\n";
    return;
  }
  int percent = 0;
  
  Frame* curFrame = &ad.frames[ad.currentFrame];
  
  curFrame->minResid = INFINITY;
  curFrame->maxResid = -INFINITY;
  curFrame->resid.clear();
  for (int i=0; i<ad.frames[0].size(); i++) {
    int p =(int)((float)i/ad.frames[0].size()*100);
    if (p != percent && p % 10 == 0) {
      cout << p << "%\n";
      percent = p;
    }
    
    NeighborLookupProc nlp = NeighborLookupProc();
    kdtree.lookup(ad.frames[0][i], nlp, radius);
    
    curFrame->resid.push_back(getResidual(ad, nlp.neighbors, ad.currentFrame, ad.currentFrame+1, false));
    
    if (curFrame->resid[i] > curFrame->maxResid) {
      curFrame->maxResid = curFrame->resid[i];
    }
    if (curFrame->resid[i] < curFrame->minResid) {
      curFrame->minResid = curFrame->resid[i];
    }
    
    writeBinary(&(curFrame->resid[i]), sizeof(float), residOutFile);
  }
  residOutFile.close();
  cout << "Done!\n";
  
}

// Create .resid files for a range of frames.
void VisualizerApp::createResidFiles(const int endFile)
{
  while (ad.currentFrame <= endFile) {
    writeResidFile();
    ad.currentFrame++;
    if (ad.frames.left_buffer_size(ad.currentFrame) < 1) {
      for (int i=1; i<NUM_FRAMES-2 && ad.currentFrame + i <= END_FRAME; i++) {
        loadFrame(ad, ad.currentFrame + i, true);
      }
    }
  }
}


CINDER_APP_NATIVE( VisualizerApp, RendererGl )
