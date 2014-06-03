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
#include "YarnBuilder.h"

using namespace ci;
using namespace ci::app;
using namespace Eigen;

class SimulatorApp : public AppNative {
  public:
	void setup();
	void mouseDown(MouseEvent event);
  void mouseDrag(MouseEvent event);
  void mouseUp(MouseEvent event);
  void mouseWheel(MouseEvent event);
  void keyDown(KeyEvent event);
  void resize();
	void update();
	void draw();
  
  void loadYarnFile(std::string filename);
  void loadDefaultYarn(int numPoints);
  void loadStdEnergies();
  
  // Set to false to pause the simulation
  bool running = true;
  
  // Camera for the scene, along with its position and orientation
  CameraPersp cam;
  ci::Vec3f eyePos = ci::Vec3f(50.0f, 0.0f, 0.0f);
  ci::Vec3f targetPos = ci::Vec3f(0.0f, 0.0f, 0.0f);
  
  // Rendering stuff
  gl::GlslProg yarnProg;
  gl::GlslProg diffuseProg;
  gl::Texture yarnTex;
  gl::Texture floorTex;
  gl::DisplayList* spheredl;
  gl::DisplayList* cylinderdl;
  gl::Material m = gl::Material(Color(0.3f, 0.3f, 0.3f), Color(0.8f, 0.9f, 0.9f));
  gl::Light* l;
  TriMesh floor;
  
  Yarn* y = nullptr;
  Clock c;
  Integrator* integrator = nullptr;
  std::vector<YarnEnergy*> energies;
  MouseSpring* mouseSpring;
  
  Spring* testSpring1;
  Spring* testSpring2;
  Eigen::Vector3f testSpring1Clamp = Eigen::Vector3f(10.0f, 15.0f, 5.0f);
  Eigen::Vector3f testSpring2Clamp = Eigen::Vector3f(-10.0f, 15.0f, 5.0f);
  
  float twist = 0.0f;
  
  // Interactive stuff
  bool isMouseDown = false;
  ci::Vec3f mousePosition;
  bool isRotate = false;
};


void SimulatorApp::setup()
{
  // YarnBuilder::buildBraid();
  
  // Setup scene
  cam.setPerspective(40.0f, getWindowAspectRatio(), 0.1f, 1000.0f);
  cam.lookAt(eyePos, targetPos, ci::Vec3f(0.0f, 1.0f, 0.0f));
  
  // Setup rendering stuff
  spheredl = new gl::DisplayList(GL_COMPILE);
  spheredl->newList();
  gl::drawSphere(ci::Vec3f::zero(), constants::radius);
  spheredl->endList();
  
  cylinderdl = new gl::DisplayList(GL_COMPILE);
  cylinderdl->newList();
  gl::drawCylinder(constants::radius, constants::radius, 1.0f);
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
  
  floor.appendVertex(ci::Vec3f(-100.0f, 0.0f, -100.0f));
  floor.appendNormal(ci::Vec3f(0.0f, 1.0f, 0.0f));
  floor.appendTexCoord(ci::Vec2f(-12.0f, -12.0f));
  floor.appendVertex(ci::Vec3f(100.0f, 0.0f, -100.0f));
  floor.appendNormal(ci::Vec3f(0.0f, 1.0f, 0.0f));
  floor.appendTexCoord(ci::Vec2f(12.0f, -12.0f));
  floor.appendVertex(ci::Vec3f(100.0f, 0.0f, 100.0f));
  floor.appendNormal(ci::Vec3f(0.0f, 1.0f, 0.0f));
  floor.appendTexCoord(ci::Vec2f(12.0f, 12.0f));
  floor.appendVertex(ci::Vec3f(-100.0f, 0.0f, 100.0f));
  floor.appendNormal(ci::Vec3f(0.0f, 1.0f, 0.0f));
  floor.appendTexCoord(ci::Vec2f(-12.0f, 12.0f));
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
  
  // Load the yarn
  loadDefaultYarn(42);
  // loadYarnFile("");
  loadStdEnergies();
}

void SimulatorApp::mouseDown(MouseEvent event)
{
  if (event.isRight()) {
    // Set targetPos to the ControlPoint we just clicked
    if (!y) return;
    Vec2i mouse = event.getPos();
    Vec2i windowSize = getWindowSize();
    Ray r = cam.generateRay((float)mouse.x/windowSize.x,
                            1.0f - (float)mouse.y/windowSize.y,
                            getWindowAspectRatio());
    float tmin = INFINITY;
    bool any = false;
    for (const CtrlPoint& p : y->cur().points) { // A bit slow, but beats keeping a KD-Tree updated
      Sphere s(toCi(p.pos), constants::radius * 1.5f);
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

void SimulatorApp::mouseDrag(MouseEvent event)
{
  if (!running) return;
  Vec2i mouse = event.getPos();
  Vec2i windowSize = getWindowSize();
  
  Ray r = cam.generateRay((float)mouse.x/windowSize.x,
                          1.0f - (float)mouse.y/windowSize.y,
                          getWindowAspectRatio());
  
  float t;
  if(!r.calcPlaneIntersection(targetPos, targetPos-eyePos, &t)) {
    std::cerr << "Mouse ray did not intersect plane!\n";
  }
  mousePosition = r.calcPosition(t);
}

void SimulatorApp::mouseUp(MouseEvent event)
{
  if (!running) return;
  isMouseDown = false;
  
}

void SimulatorApp::mouseWheel(MouseEvent event) {
  float scroll = event.getWheelIncrement();
  eyePos += (targetPos - eyePos).normalized() * scroll;
  cam.lookAt(eyePos, targetPos);
}

void SimulatorApp::keyDown(KeyEvent event)
{
  switch (event.getCode()) {
    case event.KEY_r:
      if (event.isShiftDown()) {
        twist = 0.0f;
      } else {
        isRotate = !isRotate;
      }
      break;
    case event.KEY_p:
      running = !running;
      break;
    case event.KEY_y:
      running = false;
      loadYarnFile("");
      loadStdEnergies();
      break;
    case event.KEY_LEFT:
    {
      ci::Vec2f v(eyePos.x - targetPos.x, eyePos.z - targetPos.z);
      eyePos.x = v.x*cosf(0.2f) - v.y*sinf(0.2f) + targetPos.x;
      eyePos.z = v.x*sinf(0.2f) + v.y*cosf(0.2f) + targetPos.z;
      cam.lookAt(eyePos, targetPos);
    }
      break;
    case event.KEY_RIGHT:
    {
      ci::Vec2f v(eyePos.x - targetPos.x, eyePos.z - targetPos.z);
      eyePos.x = v.x*cosf(0.2f) + v.y*sinf(0.2f) + targetPos.x;
      eyePos.z = - v.x*sinf(0.2f) + v.y*cosf(0.2f) + targetPos.z;
      cam.lookAt(eyePos, targetPos);
    }
      break;
    case event.KEY_UP:
      eyePos.y += 1.0f;
      cam.lookAt(eyePos, targetPos);
      break;
      case event.KEY_DOWN:
      eyePos.y -= 1.0f;
      cam.lookAt(eyePos, targetPos);
      break;
      
      
    case event.KEY_w:
      testSpring1Clamp.z() += 1.0f;
      testSpring2Clamp.z() += 1.0f;
      testSpring1->setClamp(testSpring1Clamp);
      testSpring2->setClamp(testSpring2Clamp);
      break;
    case event.KEY_s:
      testSpring1Clamp.z() -= 1.0f;
      testSpring2Clamp.z() -= 1.0f;
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
  
  while (!integrator->integrate(*y, c)) {
    if (c.canDecreaseTimestep()) {
      c.suggestTimestep(c.timestep() / 2.0f);
    } else {
      std::cout << "Simulation Failed!\n";
      running = false;
    }
  }
  
  /// Update Bishop frame
  for(int i=0; i<y->numSegs(); i++) {
    const Segment& prevTimeSeg = y->cur().segments[i];
    Segment& seg = y->next().segments[i];
    
    if (i == 0) {
      seg.parallelTransport(prevTimeSeg);
    } else {
      const Segment& prevSpaceSeg = y->next().segments[i-1];
      seg.parallelTransport(prevTimeSeg, prevSpaceSeg);
    }
  }
  
  /// Update material frame rotation
  if (isRotate) {
    twist += 2.0f*constants::pi*c.timestep();
  }
  y->next().segments[y->numSegs()-1].setRot(twist);
  if (!integrator->setRotations(*y)) {
    std::cout << "rotations failed";
  }
  
  // Swap Yarns
  y->swapYarns();

  c.increment();
}

void SimulatorApp::draw()
{
	// Clear out the window with grey
	gl::clear(Color(0.45f, 0.45f, 0.5f));
  
  // Enable alpha blending and depth testing
  gl::enableAlphaBlending();
	gl::enableDepthRead(true);
	gl::enableDepthWrite(true);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
  // Draw framerate counter
  gl::setMatricesWindow(getWindowSize());
  std::stringstream ss;
  ss << getAverageFps();
  gl::drawStringRight(ss.str(),
                      ci::Vec2f(getWindowWidth()-toPixels(10), getWindowHeight()-toPixels(20)),
                      Color(0.0f, 0.0f, 0.0f),
                      Font("Arial", toPixels(12)));
  
  // Set projection/modelview matrices
  gl::setMatrices(cam);
  
  // Draw the rod and the normal of the bishop frame
  for(int i=0; i<y->numSegs(); i++) {
    ci::Vec3f p0 = toCi(y->cur().points[i].pos);
    ci::Vec3f p1 = toCi(y->cur().points[i+1].pos);
    gl::drawLine(p0, p1);
    gl::color(1.0f, 1.0f, 0.0f);
    gl::lineWidth(1.0f);
    ci::Vec3f u = toCi(y->cur().segments[i].getU());
    gl::drawLine((p0+p1)/2.0f, (p0+p1)/2.0f+u);
  }
  

  m.apply();
  
  l->setDiffuse(Color::white());
  l->setAmbient(Color::white());
  l->setPosition(ci::Vec3f(0.0f, 50.0f, 0.0f));
  l->enable();
  
  diffuseProg.bind();
  for (int i=0; i<y->numCPs(); i++) {
    gl::pushModelView();
    gl::translate(toCi(y->cur().points[i].pos));
    spheredl->draw();
    gl::popModelView();
  }
  diffuseProg.unbind();
  
  yarnProg.bind();

  floorTex.enableAndBind();
  gl::draw(floor);
  floorTex.disable();
  
  yarnTex.enableAndBind();
  
  for (int i=0; i<y->numSegs(); i++) {
    gl::pushModelView();
    const Segment& s = y->cur().segments[i];
    ci::Vec3f v = toCi(s.vec().normalized());
    
    gl::translate(toCi(s.getFirst().pos));
    Quatf q(ci::Vec3f(0.0f, 1.0f, 0.0f), v);
    float angle = acosf(std::max(-1.0f, std::min(1.0f, (q*ci::Vec3f(-1.0f, 0.0f, 0.0f)).dot(toCi(s.getU())))));
    if ((q*ci::Vec3f(-1.0f, 0.0f, 0.0f)).dot(toCi(s.v())) > 0.0f) angle *= -1.0f;
    gl::rotate(Quatf(v, angle));
    gl::rotate(q);
    gl::rotate(ci::Vec3f(0.0f, s.getRot()*180.0f/constants::pi, 0.0f));
    gl::scale(1.0f, s.length(), 1.0f);
    cylinderdl->draw();
    gl::popModelView();
  }

  yarnTex.unbind();
  yarnProg.unbind();
  
  
  for (YarnEnergy* e : energies) {
    e->draw();
  }
  integrator->draw();
  
}

void SimulatorApp::loadYarnFile(std::string filename) {
  if (filename.empty()) filename = getOpenFilePath().string();
  if (filename.empty()) return;
  
  std::ifstream yarnFile(filename);
  
  if (!yarnFile.is_open()) {
    std::cerr << filename << " failed to open!\n";
    return;
  }
  
  std::string line;
  std::getline(yarnFile, line);
  const int numPoints = std::stoi(line);
  
  std::vector<Eigen::Vector3f> yarnPoints;
  yarnPoints.reserve(numPoints);
  
  for (int j=0; j<numPoints; j++) {
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
  Eigen::Vector3f u;
  for (int i=0; i<3; i++) {
    std::string line;
    std::getline(yarnFile, line);
    u(i) = std::stof(line);
  }
  assert((yarnPoints[1] - yarnPoints[0]).dot(u) < 5e-6f);
  
  yarnFile.close();
  if (y) delete y;
  y = new Yarn(yarnPoints, u);
}

void SimulatorApp::loadDefaultYarn(int numPoints) {
  if (y) delete y;
  
  std::vector<Eigen::Vector3f> yarnPoints;
  for(int i=0; i < numPoints; i++) {
    Eigen::Vector3f p(0.0f, (numPoints-i)*20.0f/numPoints, 0.0f);
    yarnPoints.push_back(p);
  }
  
  eyePos = ci::Vec3f(40.0f, 10.0f, 0.0f);
  targetPos = ci::Vec3f(0.0f, 10.0f, 0.0f);
  cam.lookAt(eyePos, targetPos, ci::Vec3f(0.0f, 1.0f, 0.0f));
  
  y = new Yarn(yarnPoints, Eigen::Vector3f(0.0f, 0.0f, 1.0f));
}

void SimulatorApp::loadStdEnergies() {
  // Create Yarn Energies - Add in the order they are most likely to fail during evaluation
  assert(y && "Tried to load evergies on a null yarn");
  energies.clear();
  
  YarnEnergy* stretch = new Stretching(*y, Implicit);
  energies.push_back(stretch);
  
  YarnEnergy* bending = new Bending(*y, Explicit);
  energies.push_back(bending);
  
  YarnEnergy* twisting = new Twisting(*y, Explicit);
  energies.push_back(twisting);
  
  YarnEnergy* gravity = new Gravity(*y, Explicit, Eigen::Vector3f(0.0f, -9.8f, 0.0f));
  energies.push_back(gravity);
  
  mouseSpring = new MouseSpring(*y, Explicit, y->numCPs()-1, 100.0f);
  energies.push_back(mouseSpring);
  
  YarnEnergy* intContact = new IntContact(*y, Explicit);
  energies.push_back(intContact);
  
  Spring* clamp1 = new Spring(*y, Implicit, 0, 500.0f);
  clamp1->setClamp(y->rest().points[0].pos);
//  clamp1->setClamp(y->rest().points[0].pos + Eigen::Vector3f(0.0f, 6.0f, 2.0f));
  Spring* clamp2 = new Spring(*y, Implicit, 1, 1000.0f);
  clamp2->setClamp(y->rest().points[1].pos);
//  Spring* clamp2 = new Spring(*y, Implicit, 14, 500.0f);
//  clamp2->setClamp(y->rest().points[14].pos + Eigen::Vector3f(0.0f, -6.0f, 2.0f));
  Spring* clamp3 = new Spring(*y, Implicit, 28, 500.0f);
  clamp3->setClamp(y->rest().points[28].pos + Eigen::Vector3f(0.0f, 6.0f, -2.0f));
  Spring* clamp4 = new Spring(*y, Implicit, 42, 500.0f);
  clamp4->setClamp(y->rest().points[42].pos + Eigen::Vector3f(0.0f, -6.0f, -2.0f));
  energies.push_back(clamp1);
  energies.push_back(clamp2);
//  energies.push_back(clamp3);
//  energies.push_back(clamp4);
  
  testSpring1 = new Spring(*y, Explicit, 2*y->numCPs()/3, 50.0f);
  testSpring1->setClamp(testSpring1Clamp);
  testSpring2 = new Spring(*y, Explicit, y->numCPs()-1, 50.0f);
  testSpring2->setClamp(testSpring2Clamp);
//  energies.push_back(testSpring1);
//  energies.push_back(testSpring2);
  
  if (integrator) delete integrator;
  integrator = new Integrator(energies);
}


CINDER_APP_NATIVE(SimulatorApp, RendererGl)
