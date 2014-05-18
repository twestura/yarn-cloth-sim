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
#include "cinder/gl/Texture.h"
#include "cinder/TriMesh.h"
#include "cinder/Sphere.h"
#include "cinder/Camera.h"

#include "Resources.h"
#include "Util.h"
#include "Yarn.h"
#include "Clock.h"
#include "Integrator.h"

using namespace ci;
using namespace ci::app;
using namespace Eigen;

class SimulatorApp : public AppNative {
  public:
	void setup();
	void mouseDown( MouseEvent event );
  void mouseDrag( MouseEvent event );
  void mouseUp( MouseEvent event );
  void mouseWheel( MouseEvent event );
  void keyDown( KeyEvent event );
  void resize();
	void update();
	void draw();
  
  void loadYarnFile(std::string filename);
  void loadDefaultYarn(int numPoints);
  
  // Set to false to pause the simulation
  bool running = true;
  
  // Camera for the scene, along with its position and orientation
  CameraPersp cam;
  ci::Vec3f eyePos = ci::Vec3f( 50, 0, 0 );
  ci::Vec3f targetPos = ci::Vec3f( 0, 0, 0 );
  
  // Rendering stuff
  gl::GlslProg yarnProg;
  gl::GlslProg diffuseProg;
  gl::Texture yarnTex;
  gl::Texture floorTex;
  gl::DisplayList* spheredl;
  gl::DisplayList* cylinderdl;
  gl::Material m = gl::Material(Color(.3, .3, .3), Color(.8, .9, .9));
  gl::Light* l;
  TriMesh floor;
  
  Yarn* y = 0;
  Clock c;
  Integrator* integrator;
  std::vector<YarnEnergy*> energies;
  MouseSpring* mouseSpring;
  
  Spring* testSpring1;
  Spring* testSpring2;
  Eigen::Vector3f testSpring1Clamp = Eigen::Vector3f(10, 15, 5);
  Eigen::Vector3f testSpring2Clamp = Eigen::Vector3f(-10, 15, 5);
  
  float twist = 0;
  
  // Interactive stuff
  bool isMouseDown = false;
  ci::Vec3f mousePosition;
  bool isRotate = false;
};


void SimulatorApp::setup()
{
  // Setup scene
  cam.setPerspective(40, getWindowAspectRatio(), .1, 1000);
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
    yarnTex = loadImage(loadResource(RES_SIM_YARN_TEX));
  } catch (ImageIoException e) {
    std::cerr << "Error loading textures: " << e.what();
    exit(1);
  }
  
  // Load and compile shaders
  try {
    diffuseProg = gl::GlslProg(loadResource(RES_SIM_VERT_GLSL), loadResource(RES_SIM_FRAG_GLSL));
    yarnProg = gl::GlslProg(loadResource(RES_SIM_VERT_TEX_GLSL), loadResource(RES_SIM_FRAG_TEX_GLSL));
  } catch (gl::GlslProgCompileExc e) {
    std::cerr << "Error compiling GLSL program: " << e.what();
    exit(1);
  } catch (ResourceLoadExc e) {
    std::cerr << "Error loading shaders: " << e.what();
    exit(1);
  }
  
  floor.appendVertex(ci::Vec3f(-100, 0, -100));
  floor.appendNormal(ci::Vec3f(0, 1, 0));
  floor.appendTexCoord(Vec2f(-12, -12));
  floor.appendVertex(ci::Vec3f(100, 0, -100));
  floor.appendNormal(ci::Vec3f(0, 1, 0));
  floor.appendTexCoord(Vec2f(12, -12));
  floor.appendVertex(ci::Vec3f(100, 0, 100));
  floor.appendNormal(ci::Vec3f(0, 1, 0));
  floor.appendTexCoord(Vec2f(12, 12));
  floor.appendVertex(ci::Vec3f(-100, 0, 100));
  floor.appendNormal(ci::Vec3f(0, 1, 0));
  floor.appendTexCoord(Vec2f(-12, 12));
  floor.appendTriangle(0, 1, 2);
  floor.appendTriangle(0, 3, 2);
  
  ci::Surface s(4, 4, false);
  auto iter = s.getIter();
  do {
    do {
      Vec2i pos = iter.getPos();
      unsigned char val = pos.x > 0 && pos.x < 3 && pos.y > 0 && pos.y < 3 ? 100 : 150;
      iter.r() = iter.g() = iter.b() = val;
    } while (iter.pixel());
  } while (iter.line());
  floorTex = gl::Texture(s);
  floorTex.setMagFilter(GL_NEAREST);
  floorTex.setWrap(GL_REPEAT, GL_REPEAT);
  
  
  // Create Yarn Energies - Add in the order they are most likely to fail during evaluation
  
  YarnEnergy* stretch = new Stretching(*y, Implicit);
  energies.push_back(stretch);
  
  YarnEnergy* bending = new Bending(*y, Explicit);
  energies.push_back(bending);
  
  YarnEnergy* twisting = new Twisting(*y, Explicit);
  energies.push_back(twisting);
  
  YarnEnergy* gravity = new Gravity(*y, Explicit, Eigen::Vector3f(0, -9.8, 0));
  energies.push_back(gravity);
  
  mouseSpring = new MouseSpring(*y, Explicit, y->numCPs()-1, 10);
  energies.push_back(mouseSpring);
  
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
  
  // Create integrator
  integrator = new Integrator(energies);
  
}

void SimulatorApp::mouseDown( MouseEvent event )
{
  if (event.isRight()) {
    // Set targetPos to the ControlPoint we just clicked
    if (!y) return;
    Vec2i mouse = event.getPos();
    Vec2i windowSize = getWindowSize();
    Ray r = cam.generateRay((float)mouse.x/windowSize.x,
                            1.0 - (float)mouse.y/windowSize.y,
                            getWindowAspectRatio());
    float tmin = INFINITY;
    bool any = false;
    for (const CtrlPoint& p : y->cur().points) { // A bit slow, but beats keeping a KD-Tree updated
      Sphere s(toCi(p.pos), constants::radius * 1.5);
      float t;
      if (s.intersect(r, &t) && t < tmin) {
        any = true;
        tmin = t;
      }
    }
    if (!any) return;
    targetPos = r.calcPosition(tmin);
    cam.lookAt(targetPos);
  } else {
    if (!running) return;
    isMouseDown = true;
    mouseDrag(event);
  }
}

void SimulatorApp::mouseDrag( MouseEvent event )
{
  if (!running) return;
  Vec2i mouse = event.getPos();
  Vec2i windowSize = getWindowSize();
  
  Ray r = cam.generateRay((float)mouse.x/windowSize.x,
                          1.0 - (float)mouse.y/windowSize.y,
                          getWindowAspectRatio());
  
  float t;
  if(!r.calcPlaneIntersection(targetPos, targetPos-eyePos, &t)) {
    std::cerr << "Mouse ray did not intersect plane!\n";
  }
  mousePosition = r.calcPosition(t);
}

void SimulatorApp::mouseUp( MouseEvent event )
{
  if (!running) return;
  isMouseDown = false;
  
}

void SimulatorApp::mouseWheel(MouseEvent event) {
  float scroll = event.getWheelIncrement();
  eyePos += (targetPos - eyePos).normalized() * scroll;
  cam.lookAt(eyePos, targetPos);
}

void SimulatorApp::keyDown( KeyEvent event )
{
  switch (event.getCode()) {
    case event.KEY_r:
      if (event.isShiftDown()) {
        twist = 0;
      } else {
        isRotate = !isRotate;
      }
      break;
    case event.KEY_p:
      running = !running;
      break;
    case event.KEY_LEFT:
    {
      Vec2f v(eyePos.x - targetPos.x, eyePos.z - targetPos.z);
      eyePos.x = v.x*cosf(0.2) - v.y*sinf(0.2) + targetPos.x;
      eyePos.z = v.x*sinf(0.2) + v.y*cosf(0.2) + targetPos.z;
      cam.lookAt(eyePos, targetPos);
    }
      break;
    case event.KEY_RIGHT:
    {
      Vec2f v(eyePos.x - targetPos.x, eyePos.z - targetPos.z);
      eyePos.x = v.x*cosf(0.2) + v.y*sinf(0.2) + targetPos.x;
      eyePos.z = - v.x*sinf(0.2) + v.y*cosf(0.2) + targetPos.z;
      cam.lookAt(eyePos, targetPos);
    }
      break;
    case event.KEY_UP:
      eyePos.y += 1;
      cam.lookAt(eyePos, targetPos);
      break;
      case event.KEY_DOWN:
      eyePos.y -= 1;
      cam.lookAt(eyePos, targetPos);
      break;
      
      
    case event.KEY_w:
      testSpring1Clamp.z() += 1;
      testSpring2Clamp.z() += 1;
      testSpring1->setClamp(testSpring1Clamp);
      testSpring2->setClamp(testSpring2Clamp);
      break;
    case event.KEY_s:
      testSpring1Clamp.z() -= 1;
      testSpring2Clamp.z() -= 1;
      testSpring1->setClamp(testSpring1Clamp);
      testSpring2->setClamp(testSpring2Clamp);
      break;
      
    default:;
  }
}

void SimulatorApp::resize() {
  cam.setAspectRatio(getWindowAspectRatio());
}

void SimulatorApp::update()
{
  if (!running) return;
  
  Eigen::Vector3f mp;
  if (isMouseDown) mp << mousePosition.x, mousePosition.y, mousePosition.z;
  mouseSpring->setMouse(mp, isMouseDown);
  
  integrator->integrate(*y, c);
  
  /// Update Bishop frame
  for(int i=0; i<y->numSegs(); i++) {
    const Segment& prevSeg = y->cur().segments[i];
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
    twist += 2*constants::pi*c.timestep();
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
	gl::clear(Color( 0.45, 0.45, 0.5 ));
  
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
  for(int i=0; i<y->numSegs(); i++) {
    ci::Vec3f p0 = toCi(y->cur().points[i].pos);
    ci::Vec3f p1 = toCi(y->cur().points[i+1].pos);
    gl::drawLine(p0, p1);
    gl::color(1, 1, 0);
    gl::lineWidth(1);
    ci::Vec3f u = toCi(y->cur().segments[i].getU());
    gl::drawLine((p0+p1)/2, (p0+p1)/2+u);
  }
  

  m.apply();
  
  l->setDiffuse(Color::white());
  l->setAmbient(Color::white());
  l->setPosition(ci::Vec3f(0, 50, 0));
  l->enable();
  
  diffuseProg.bind();
  
  for (int i=0; i<y->numCPs(); i++) {
    gl::pushModelView();
    gl::translate(toCi(y->cur().points[i].pos));
    spheredl->draw();
    gl::popModelView();
  }
  
//  gl::draw(floor);
  diffuseProg.unbind();
  
  yarnProg.bind();

  floorTex.enableAndBind();
  gl::draw(floor);
  floorTex.disable();
  
  yarnTex.enableAndBind();
  
  float totalTwist = 0;
  for (int i=0; i<y->numSegs(); i++) {
    gl::pushModelView();
    const Segment& s = y->cur().segments[i];
    ci::Vec3f v = toCi(s.vec().normalized());
    
    gl::translate(toCi(s.getFirst().pos));
    Quatf q(ci::Vec3f(0,1,0), v);
    float angle = acosf(std::max(-1.0f, std::min(1.0f, (q*ci::Vec3f(-1, 0, 0)).dot(toCi(s.getU())))));
    if ((q*ci::Vec3f(-1, 0, 0)).dot(toCi(s.v())) > 0) angle *= -1;
    gl::rotate(Quatf(v, angle));
//    gl::rotate(Quatf(q * ci::Vec3f(-1, 0, 0), toCi(s.getU())));
    gl::rotate(q);
    gl::rotate(ci::Vec3f(0, (s.getRot() - s.getRefTwist() - totalTwist)*180/constants::pi, 0));
    gl::scale(1, s.length(), 1);
    cylinderdl->draw();
    gl::popModelView();
    totalTwist += s.getRefTwist();
  }

  yarnTex.unbind();
  yarnProg.unbind();
  
  
  for (YarnEnergy* e : energies) {
    e->draw();
  }
  
}

void SimulatorApp::loadYarnFile(std::string filename) {
  if (filename.empty()) filename = getOpenFilePath().string();
  if (filename.empty()) return;
  
  std::ifstream yarnFile;

  try {
    yarnFile.open(filename);
  } catch (Exception e) {
    std::cerr << e.what() << "\n";
    return;
  }
  
  if (!yarnFile.is_open()) {
    std::cerr << filename << " failed to open!\n";
    return;
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
  
  if (y != 0) delete y;
  y = new Yarn(yarnPoints, Eigen::Vector3f(0, 0, 1));
}

void SimulatorApp::loadDefaultYarn(int numPoints) {
  if (y != 0) delete y;
  
  std::vector<Eigen::Vector3f> yarnPoints;
  for(int i=0; i < numPoints; i++) {
    Eigen::Vector3f p(0, (numPoints-i)*20.0f/numPoints, 0);
    yarnPoints.push_back(p);
  }
  
  eyePos = ci::Vec3f(40, 10, 0);
  targetPos = ci::Vec3f(0, 10, 0);
  cam.lookAt(eyePos, targetPos, ci::Vec3f(0, 1, 0));
  
  y = new Yarn(yarnPoints, Eigen::Vector3f(0, 0, 1));
}


CINDER_APP_NATIVE( SimulatorApp, RendererGl )
