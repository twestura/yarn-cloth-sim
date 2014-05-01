//
//  SimulatorApp.cpp
//  Visualizer
//
//  Created by eschweickart on 3/5/14.
//
//

#include "SimulatorApp.h"
#include "cinder/app/AppNative.h"
#include "cinder/gl/gl.h"
#include "cinder/gl/Light.h"
#include "cinder/gl/DisplayList.h"
#include "cinder/gl/GlslProg.h"
#include "cinder/Camera.h"
#include "cinder/DataSource.h"

#include "Util.h"
#include "Yarn.h"
#include "Clock.h"
#include "Integrator.h"

using namespace ci;
using namespace ci::app;
using namespace Eigen;

#define ANG_VEL (10*3.14159) // Constant angular velocity of the material frame.

class SimulatorApp : public AppNative {
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
  
  // Rendering stuff
  gl::GlslProg meshProg;
  gl::DisplayList* spheredl;
  gl::DisplayList* cylinderdl;
  gl::Material m = gl::Material(Color(.2, .2, .2), Color(.8, .9, .9));
  gl::Light* l;
  
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
};


void SimulatorApp::setup()
{
  // Setup scene
  cam.setOrtho(cameraBounds[0], cameraBounds[1], cameraBounds[2], cameraBounds[3], 1, 100);
  cam.lookAt(eyePos, targetPos, ci::Vec3f( 0, 1, 0 ));
  
  // Load the yarn
  loadDefaultYarn(42);
  
  // Setup rendering stuff
  spheredl = new gl::DisplayList(GL_COMPILE);
  spheredl->newList();
  gl::drawSphere(ci::Vec3f::zero(), constants::radius);
  spheredl->endList();
  
  cylinderdl = new gl::DisplayList(GL_COMPILE);
  cylinderdl->newList();
  gl::drawCylinder(constants::radius, constants::radius, 1);
  cylinderdl->endList();
  
  l = new gl::Light(gl::Light::POINT, 0);
  
  try {
    // FIXME: Make this portable
    DataSourcePathRef vert = DataSourcePath::create(fs::path("/Users/eschweickart/Research/Visualizer/xcode/Simulator/vert.glsl"));
    DataSourcePathRef frag = DataSourcePath::create(fs::path("/Users/eschweickart/Research/Visualizer/xcode/Simulator/frag.glsl"));
    meshProg = gl::GlslProg(vert, frag);
  } catch (gl::GlslProgCompileExc e) {
    std::cerr << "Error compiling GLSL program: " << e.what();
    exit(1);
  } catch (ci::StreamExc e) {
    std::cerr << "Error loading shaders: " << e.what();
    exit(1);
  }
  
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

void SimulatorApp::mouseDown( MouseEvent event )
{
  if (!running) return;
  isMouseDown = true;
  Vec2i mouse = event.getPos();
  Vec2i windowSize = getWindowSize();
  mousePosition.z = cameraBounds[1] * (1.0 - (float) mouse.x / windowSize.x)  + cameraBounds[0] * ((float) mouse.x / windowSize.x);
  mousePosition.y = cameraBounds[3] * (1.0 - (float) mouse.y / windowSize.y)  + cameraBounds[2] * ((float) mouse.y / windowSize.y);
}

void SimulatorApp::mouseDrag( MouseEvent event )
{
  if (!running) return;
  Vec2i mouse = event.getPos();
  Vec2i windowSize = getWindowSize();
  mousePosition.z = cameraBounds[1] * (1.0 - (float) mouse.x / windowSize.x)  + cameraBounds[0] * ((float) mouse.x / windowSize.x);
  mousePosition.y = cameraBounds[3] * (1.0 - (float) mouse.y / windowSize.y)  + cameraBounds[2] * ((float) mouse.y / windowSize.y);
}

void SimulatorApp::mouseUp( MouseEvent event )
{
  if (!running) return;
  isMouseDown = false;
  
}

void SimulatorApp::keyDown( KeyEvent event )
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
      testSpring1Clamp.z() += 1;
      testSpring2Clamp.z() += 1;
      testSpring1->setClamp(testSpring1Clamp);
      testSpring2->setClamp(testSpring2Clamp);
      break;
      case 's':
      testSpring1Clamp.z() -= 1;
      testSpring2Clamp.z() -= 1;
      testSpring1->setClamp(testSpring1Clamp);
      testSpring2->setClamp(testSpring2Clamp);
      break;
      
    default:;
  }
}

void SimulatorApp::update()
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

void SimulatorApp::draw()
{
	// Clear out the window with grey
	gl::clear(Color( 0.5, 0.5, 0.6 ));
  
  // Enable alpha blending and depth testing
  gl::enableAlphaBlending();
	gl::enableDepthRead( true );
	gl::enableDepthWrite( true );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  
  // Draw framerate counter
  gl::setMatricesWindow(getWindowSize());
  std::stringstream ss;
  ss << getAverageFps();
  gl::drawStringRight(ss.str(),
                      ci::Vec2f(getWindowWidth()-toPixels(10), getWindowHeight()-toPixels(20)),
                      Color(0, 0, 0),
                      Font("Arial", toPixels(12)));
  
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
  

  m.apply();
  
  l->setDiffuse(Color::white());
  l->setAmbient(Color::white());
  l->setPosition(ci::Vec3f(0, 50, 0));
  l->enable();
  
  meshProg.bind();
  
  for (int i=0; i<y->numCPs(); i++) {
    gl::pushModelView();
    gl::translate(y->cur().points[i].pos.x(), y->cur().points[i].pos.y(), y->cur().points[i].pos.z());
    spheredl->draw();
    gl::popModelView();
  }
  
  for (int i=0; i<y->numSegs(); i++) {
    gl::pushModelView();
    Segment& s = y->cur().segments[i];
    gl::translate(s.getFirst().pos.x(), s.getFirst().pos.y(), s.getFirst().pos.z());
    Eigen::Vector3f seg = s.vec().normalized();
    Quatf q(ci::Vec3f(0,1,0), ci::Vec3f(seg.x(), seg.y(), seg.z()));
    gl::rotate(q);
    gl::scale(1, s.length(), 1);
    cylinderdl->draw();
    gl::popModelView();
  }

  meshProg.unbind();
  
  for (YarnEnergy* e : energies) {
    e->draw();
  }
  
}

void SimulatorApp::loadYarnFile(std::string filename) {
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

void SimulatorApp::loadDefaultYarn(int numPoints) {
  if (y != 0) delete y;
  
  std::vector<Eigen::Vector3f> yarnPoints;
  for(int i=0; i < numPoints; i++) {
    Eigen::Vector3f p(0, (numPoints-i)*20.0f/numPoints-10, 0);
    yarnPoints.push_back(p);
  }
  
  y = new Yarn(yarnPoints, Eigen::Vector3f(0, 0, 1));
  
}




CINDER_APP_NATIVE( SimulatorApp, RendererGl )
