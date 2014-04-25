//
//  ThreadTestApp.cpp
//  Visualizer
//
//  Created by eschweickart on 3/5/14.
//
//

#include "ThreadTestApp.h"
#include "cinder/app/AppNative.h"
#include "cinder/gl/gl.h"
#include "cinder/Camera.h"
#include <fstream>

#include "Util.h"
#include "Yarn.h"
#include "Clock.h"
#include "Integrator.h"

using namespace ci;
using namespace ci::app;
using namespace Eigen;

#define ANG_VEL (10*3.14159) // Constant angular velocity of the material frame.

class ThreadTestApp : public AppNative {
  public:
	void setup();
	void mouseDown( MouseEvent event );
  void mouseDrag( MouseEvent event );
  void mouseUp( MouseEvent event );
  void keyDown( KeyEvent event );
	void update();
	void draw();
  void setMaterialFrame();
  
  void loadYarnFile(std::string filename);
  void loadDefaultYarn(int numPoints);
  
  // Set to false to pause the simulation
  bool running = true;
  
  // Camera for the scene, along with its position and orientation
  CameraOrtho cam = CameraOrtho();
  ci::Vec4f cameraBounds = ci::Vec4f(-15, 15, -15, 15);
  ci::Vec3f eyePos = ci::Vec3f( 50, 0, 0 );
  ci::Vec3f targetPos = ci::Vec3f( 0, 0, 0 );
  
  Yarn* y = 0;
  Clock c;
  Integrator* integrator;
  std::vector<YarnEnergy*> energies;
  MouseSpring* mouseSpring;
  
  Spring* testSpring1;
  Spring* testSpring2;
  Eigen::Vector3f testSpring1Clamp = Eigen::Vector3f(10, 5, 5);
  Eigen::Vector3f testSpring2Clamp = Eigen::Vector3f(-10, 5, 5);
  
  float twist = 0;
  
  // Interactive stuff
  bool isMouseDown = false;
  ci::Vec3f mousePosition = ci::Vec3f(0, 0, 0);
  bool isRotate = false;
  
  Eigen::Vector3f p0 = Eigen::Vector3f(0, 0, -10);
  
};


void ThreadTestApp::setup()
{
  // Setup scene
  cam.setOrtho(cameraBounds[0], cameraBounds[1], cameraBounds[2], cameraBounds[3], 1, 100);
  cam.lookAt(eyePos, targetPos, ci::Vec3f( 0, 1, 0 ));
  
  loadDefaultYarn(42);
  
  // Create Yarn Energies
  YarnEnergy* gravity = new Gravity(*y, Explicit, Eigen::Vector3f(0, -9.8, 0));
  energies.push_back(gravity);
  
  mouseSpring = new MouseSpring(*y, Explicit, y->numCPs()-1, 10);
  energies.push_back(mouseSpring);
  
  YarnEnergy* bending = new Bending(*y, Implicit);
  energies.push_back(bending);
  
  YarnEnergy* twisting = new Twisting(*y, Explicit);
  energies.push_back(twisting);
  
  YarnEnergy* intContact = new IntContact(*y, Implicit);
  energies.push_back(intContact);
  
  Spring* clamp1 = new Spring(*y, Implicit, 0, 1000);
  clamp1->setClamp(y->rest().points[0].pos);
  Spring* clamp2 = new Spring(*y, Implicit, 1, 1000);
  clamp2->setClamp(y->rest().points[1].pos);
  energies.push_back(clamp1);
  energies.push_back(clamp2);
  
  testSpring1 = new Spring(*y, Explicit, 2*y->numCPs()/3, 50);
  testSpring1->setClamp(testSpring1Clamp);
  testSpring2 = new Spring(*y, Explicit, y->numCPs()-1, 50);
  testSpring2->setClamp(testSpring2Clamp);
  energies.push_back(testSpring1);
  energies.push_back(testSpring2);
  
  YarnEnergy* stretch = new Stretching(*y, Implicit);
  energies.push_back(stretch);
  
  // Create integrator
  integrator = new Integrator(energies);
  
}

void ThreadTestApp::mouseDown( MouseEvent event )
{
  if (!running) return;
  isMouseDown = true;
  Vec2i mouse = event.getPos();
  Vec2i windowSize = getWindowSize();
  mousePosition.z = cameraBounds[1] * (1.0 - (float) mouse.x / windowSize.x)  + cameraBounds[0] * ((float) mouse.x / windowSize.x);
  mousePosition.y = cameraBounds[3] * (1.0 - (float) mouse.y / windowSize.y)  + cameraBounds[2] * ((float) mouse.y / windowSize.y);
}

void ThreadTestApp::mouseDrag( MouseEvent event )
{
  if (!running) return;
  Vec2i mouse = event.getPos();
  Vec2i windowSize = getWindowSize();
  mousePosition.z = cameraBounds[1] * (1.0 - (float) mouse.x / windowSize.x)  + cameraBounds[0] * ((float) mouse.x / windowSize.x);
  mousePosition.y = cameraBounds[3] * (1.0 - (float) mouse.y / windowSize.y)  + cameraBounds[2] * ((float) mouse.y / windowSize.y);
}

void ThreadTestApp::mouseUp( MouseEvent event )
{
  if (!running) return;
  isMouseDown = false;
  
}

void ThreadTestApp::keyDown( KeyEvent event )
{
  char input = event.getChar();
  switch (input) {
    case 'r':
    isRotate = !isRotate;
    break;
    case 'R':
    twist = 0;
    break;
    case 'p':
    running = !running;
    break;
      
      case 'w':
      // p0.y() += 1;
      testSpring1Clamp.z() += 1;
      testSpring2Clamp.z() += 1;
      testSpring1->setClamp(testSpring1Clamp);
      testSpring2->setClamp(testSpring2Clamp);
      break;
      case 's':
      // p0.y() -= 1;
      testSpring1Clamp.z() -= 1;
      testSpring2Clamp.z() -= 1;
      testSpring1->setClamp(testSpring1Clamp);
      testSpring2->setClamp(testSpring2Clamp);
      break;
      
    default:;
  }
}

void ThreadTestApp::update()
{
  if (!running) return;
  
  Eigen::Vector3f mp;
  if (isMouseDown) {
    mp << mousePosition.x+1, mousePosition.y, mousePosition.z;
  } else {
    mp = y->cur().points[y->numCPs()-1].pos;
  }
  mouseSpring->setMouse(mp, isMouseDown);
  
  integrator->integrate(*y, c);
  
  /// Update Bishop frame
  for(int i=1; i<y->numSegs(); i++) {
    Segment& prevSeg = y->cur().segments[i];
    Segment& seg = y->next().segments[i];
    
    if (i < y->numSegs()-1) {
      Segment& refSeg = y->next().segments[i+1];
      seg.parallelTransport(prevSeg, refSeg);
    } else {
      seg.parallelTransport(prevSeg);
    }
  }
  
  /// Update material frame rotation
  if (isRotate) {
    twist += ANG_VEL*c.timestep();
  }
  
  float delta = twist/y->numCPs(); // Twist is constant and instant
  for(int i=1; i<y->numSegs(); i++) {
    y->next().segments[i].setRot(y->next().segments[i-1].getRot() + delta);
  }
  
  // Swap Yarns
  y->swapYarns();

  c.increment();
}

void ThreadTestApp::draw()
{
	// Clear out the window with grey
	gl::clear( Color( 0.5, 0.5, 0.6 ) );
  
  // Enable alpha blending and depth testing
  gl::enableAlphaBlending();
	gl::enableDepthRead( true );
	gl::enableDepthWrite( true );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  
  // Draw framerate counter
  gl::setMatricesWindow(getWindowSize());
  std::stringstream ss;
  ss << getAverageFps();
  gl::drawStringRight(ss.str(), ci::Vec2f(getWindowWidth() - toPixels(10), getWindowHeight()- toPixels(20)), Color(0, 0, 0), Font("Arial", toPixels(12)));
  
  // Set projection/modelview matrices
  gl::setMatrices(cam);
  
  // Draw the rod and the normal of the bishop frame
  for(int i=1; i<y->numSegs(); i++) {
    ci::Color c(0.4, -cos(y->cur().segments[i].getRot()), cos(y->cur().segments[i].getRot()));
    c *= (y->cur().segments[i].getFirst().pos.x() + 5) / 10;
    gl::color(c);
    gl::lineWidth(2);
    ci::Vec3f p0(y->cur().points[i].pos.x(), y->cur().points[i].pos.y(), y->cur().points[i].pos.z());
    ci::Vec3f p1(y->cur().points[i+1].pos.x(), y->cur().points[i+1].pos.y(), y->cur().points[i+1].pos.z());
    gl::drawLine(p0, p1);
    gl::color(1, 1, 0);
    gl::lineWidth(1);
    Eigen::Vector3f eu = y->cur().segments[i].getU();
    ci::Vec3f u(eu.x(), eu.y(), eu.z());
    gl::drawLine((p0+p1)/2, (p0+p1)/2+u);
  }
  
  if (isMouseDown) {
    gl::color(1, 0, 0);
    ci::Vec3f p(y->cur().points[y->numCPs()-1].pos.x(), y->cur().points[y->numCPs()-1].pos.y(), y->cur().points[y->numCPs()-1].pos.z());
    gl::drawLine(p, mousePosition);
  }

  /*
  gl::lineWidth(5);
  gl::color(1, 0, 0);
  CtrlPoint splinepts[4];
  splinepts[0].pos = p0;
  splinepts[1].pos = Eigen::Vector3f(0, 10, -8);
  splinepts[2].pos = Eigen::Vector3f(0, 8, 0);
  splinepts[3].pos = Eigen::Vector3f(0, -2, 4);
  
  for (int i=0; i<4; i++) {
    ci::Vec3f p(splinepts[i].pos.x(), splinepts[i].pos.y(), splinepts[i].pos.z());
    gl::drawLine(p, p + ci::Vec3f(0, 0.1, 0));
  }
  
  Spline s(splinepts[0], splinepts[1], splinepts[2], splinepts[3]);
  std::vector<ci::Vec3f> pl;
  std::vector<ci::Vec3f> plTest;
  for (int i=0; i<=12; i++) {
    float t = (((float) i) / 12);
    Eigen::Vector3f p = s.eval(t, false);
    ci::Vec3f v(p.x(), p.y(), p.z());
    pl.push_back(v);
    gl::drawSphere(v, .05);
    p = s.eval(t, true);
    ci::Vec3f v1(p.x(), p.y(), p.z());
    plTest.push_back(v1);
    gl::drawSphere(v1, .05);
  }
  gl::lineWidth(1);
  

  for (int i=0; i<pl.size()-1; i++) {
    gl::color(0, 1, 0);
    gl::drawLine(pl[i], pl[i+1]);
    gl::color(0, 1, 1);
    gl::drawLine(plTest[i], plTest[i+1]);
  }
  */
  
}

void ThreadTestApp::loadYarnFile(std::string filename) {
  if (y != 0) delete y;
  if (filename.empty()) filename = getOpenFilePath().string();
  
  std::ifstream yarnFile;

  try {
    yarnFile.open(filename);
  } catch (Exception e) {
    std::cerr << e.what() << "\n";
    exit(1);
  }
  
  if (!yarnFile.is_open()) {
    std::cerr << filename << " failed to open!\n";
    exit(1);
  }
  
  std::vector<Eigen::Vector3f> yarnPoints;
  while (!yarnFile.eof()) {
    Eigen::Vector3f p;
    for (int i=0; i<3; i++) {
      std::string line;
      std::getline(yarnFile, line);
      if(!line.empty()) {
        p(i) = std::stof(line);
      }
    }
    yarnPoints.push_back(p);
  }
  
  y = new Yarn(yarnPoints, Eigen::Vector3f(0, 0, 1));
}

void ThreadTestApp::loadDefaultYarn(int numPoints) {
  if (y != 0) delete y;
  
  std::vector<Eigen::Vector3f> yarnPoints;
  for(int i=0; i < numPoints; i++) {
    Eigen::Vector3f p(0, (numPoints-i)*20.0f/numPoints-10, 0);
    yarnPoints.push_back(p);
  }
  
  y = new Yarn(yarnPoints, Eigen::Vector3f(0, 0, 1));
  
}




CINDER_APP_NATIVE( ThreadTestApp, RendererGl )
