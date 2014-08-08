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
#include "IMEXIntegrator.h"
#include "ConstraintIntegrator.h"
#include "YarnBuilder.h"

using namespace ci;
using namespace ci::app;

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
  void loadStdEnergiesAndConsts();
  
  // Set to false to pause the simulation
  bool running = true;
  
  // Camera for the scene, along with its position and orientation
  CameraPersp cam;
  Vec3c eyePos = Vec3c(50.0, 0.0, 0.0);
  Vec3c targetPos = Vec3c(0.0, 0.0, 0.0);
  
  // Rendering stuff
  gl::GlslProg yarnProg;
  gl::GlslProg diffuseProg;
  gl::Texture yarnTex;
  gl::Texture floorTex;
  gl::DisplayList* spheredl;
  gl::DisplayList* cylinderdl;
  gl::Material m = gl::Material(Color(0.3, 0.3, 0.3), Color(0.8, 0.9, 0.9));
  gl::Light* l;
  TriMesh floor;
  
  Yarn* y = nullptr;
  Clock c;
  Integrator* integrator = nullptr;
  ConstraintIntegrator* cIntegrator = nullptr;
  std::vector<YarnEnergy*> energies;
  MouseSpring* mouseSpring;
  std::vector<YarnConstraint*> constraints;
  
  real twist = 0.0;
  real yarnTwist = 0.0;
  int numYarnTwists = 0;
  
  // Interactive stuff
  bool isMouseDown = false;
  Vec3c mousePosition;
  bool isRotate = false;
};

void SimulatorApp::setup()
{
  // Setup scene
  cam.setPerspective(40.0, getWindowAspectRatio(), 0.1, 1000.0);
  cam.lookAt(eyePos, targetPos, Vec3c(0.0, 1.0, 0.0));
  
  // Setup rendering stuff
  spheredl = new gl::DisplayList(GL_COMPILE);
  spheredl->newList();
  gl::drawSphere(Vec3c::zero(), constants::radius);
  spheredl->endList();
  
  cylinderdl = new gl::DisplayList(GL_COMPILE);
  cylinderdl->newList();
  gl::drawCylinder(constants::radius, constants::radius, 1.0);
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
  
  floor.appendVertex(Vec3c(-100.0, 0.0, -100.0));
  floor.appendNormal(Vec3c(0.0, 1.0, 0.0));
  floor.appendTexCoord(Vec2c(-12.0, -12.0));
  floor.appendVertex(Vec3c(100.0, 0.0, -100.0));
  floor.appendNormal(Vec3c(0.0, 1.0, 0.0));
  floor.appendTexCoord(Vec2c(12.0, -12.0));
  floor.appendVertex(Vec3c(100.0, 0.0, 100.0));
  floor.appendNormal(Vec3c(0.0, 1.0, 0.0));
  floor.appendTexCoord(Vec2c(12.0, 12.0));
  floor.appendVertex(Vec3c(-100.0, 0.0, 100.0));
  floor.appendNormal(Vec3c(0.0, 1.0, 0.0));
  floor.appendTexCoord(Vec2c(-12.0, 12.0));
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
#ifndef CONST_INTEGRATOR
  loadStdEnergies();
#else
  loadStdEnergiesAndConsts();
#endif // ifndef CONST_INTEGRATOR
}

void SimulatorApp::mouseDown(MouseEvent event)
{
  if (event.isRight()) {
    // Set targetPos to the ControlPoint we just clicked
    if (!y) return;
    Vec2i mouse = event.getPos();
    Vec2i windowSize = getWindowSize();
    Ray r = cam.generateRay((real)mouse.x/windowSize.x,
                            1.0 - (real)mouse.y/windowSize.y,
                            getWindowAspectRatio());
    real tmin = INFINITY;
    bool any = false;
    for (int i=0; i<y->numCPs(); i++) { // A bit slow, but beats keeping a KD-Tree updated
      Sphere s(EtoC(y->cur().POS(i)), constants::radius * 1.5);
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
  
  Ray r = cam.generateRay((real)mouse.x/windowSize.x,
                          1.0 - (real)mouse.y/windowSize.y,
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
  real scroll = event.getWheelIncrement();
  eyePos += (targetPos - eyePos).normalized() * scroll;
  cam.lookAt(eyePos, targetPos);
}

void SimulatorApp::keyDown(KeyEvent event)
{
  switch (event.getCode()) {
    case event.KEY_r:
      if (event.isShiftDown()) {
        twist = 0.0;
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
      Vec2c v(eyePos.x - targetPos.x, eyePos.z - targetPos.z);
      eyePos.x = v.x*cosf(0.2) - v.y*sinf(0.2) + targetPos.x;
      eyePos.z = v.x*sinf(0.2) + v.y*cosf(0.2) + targetPos.z;
      cam.lookAt(eyePos, targetPos);
    }
      break;
    case event.KEY_RIGHT:
    {
      Vec2c v(eyePos.x - targetPos.x, eyePos.z - targetPos.z);
      eyePos.x = v.x*cosf(0.2) + v.y*sinf(0.2) + targetPos.x;
      eyePos.z = - v.x*sinf(0.2) + v.y*cosf(0.2) + targetPos.z;
      cam.lookAt(eyePos, targetPos);
    }
      break;
    case event.KEY_UP:
      eyePos.y += 1.0;
      cam.lookAt(eyePos, targetPos);
      break;
      case event.KEY_DOWN:
      eyePos.y -= 1.0;
      cam.lookAt(eyePos, targetPos);
      break;
      
      
    case event.KEY_w:
      testSpring1Clamp.z() += 1.0;
      testSpring2Clamp.z() += 1.0;
      testSpring1->setClamp(testSpring1Clamp);
      testSpring2->setClamp(testSpring2Clamp);
      break;
    case event.KEY_s:
      testSpring1Clamp.z() -= 1.0;
      testSpring2Clamp.z() -= 1.0;
      testSpring1->setClamp(testSpring1Clamp);
      testSpring2->setClamp(testSpring2Clamp);
      break;
      
    case event.KEY_t:
      for(int i=1; i<y->numSegs(); i++) {
        std::cout << "seg " << i << " twist: " << (y->cur().rot(i) - y->cur().rot(i-1) + y->cur().refTwist(i))
        << " (" << y->cur().rot(i) << " + " << y->cur().refTwist(i) <<
        " = " << y->cur().rot(i) + y->cur().refTwist(i) << ")\n";
      }
      std::cout << "numYarnTwists: " << numYarnTwists << "\n";
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
  
  Vec3e mp;
  if (isMouseDown) mp << mousePosition.x, mousePosition.y, mousePosition.z;
  mouseSpring->setMouse(mp, isMouseDown);
  
#ifndef CONST_INTEGRATOR
  while (!integrator->integrate(c)) {
    if (c.canDecreaseTimestep()) {
      c.suggestTimestep(c.timestep() / 2.0);
      std::cout << "Decreasing timestep: " << c.timestep() << "\n";
    } else {
      std::cout << "Simulation Failed!\n";
      running = false;
      return;
    }
  }
#else
  while (!cIntegrator->integrate(c)) {
    std::cout << "wat\n";
    throw;
  }
#endif // ifdef CONST_INTEGRATOR
  
  /// Update Bishop frame
  y->next().updateReferenceFrames(y->cur());
  
#ifndef CONST_INTEGRATOR
  /// Update material frame rotation
  if (isRotate) {
    twist += 2.0*constants::pi*c.timestep();
  }
  
  Vec3e uRef = parallelTransport(y->next().vec(0), y->next().vec(y->numSegs()-1), y->next().u[0]);
  real cosTwist = y->next().u[y->numSegs()-1].dot(uRef.normalized());
  real oldTwist = yarnTwist;
  if (cosTwist >= 1.0) { // Avoid values like 1.0000000012 that introduce NaNs
    yarnTwist = 0.0;
  } else if (cosTwist <= -1.0) {
    yarnTwist = constants::pi;
  } else {
    yarnTwist = acos(cosTwist);
  }
  // Flip the sign if necessary
  if (y->next().v(y->numSegs()-1).dot(uRef) > 0.0) {
    yarnTwist = -yarnTwist;
  }
  real diff = yarnTwist - oldTwist;
  if (diff < -constants::pi) {
    numYarnTwists += 1;
  } else if (diff > constants::pi) {
    numYarnTwists -= 1;
  }
  y->next().rot(y->numSegs()-1) = twist - yarnTwist;
  if (!static_cast<IMEXIntegrator*>(integrator)->setRotations()) {
    std::cout << "rotations failed";
  }
#endif // ifndef CONST_INTEGRATOR
  
  // Swap Yarns
  y->swapYarns();

  c.increment();
}

void SimulatorApp::draw() {
	// Clear out the window with grey
	gl::clear(Color(0.45, 0.45, 0.5));
  
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
                      Vec2c(getWindowWidth()-toPixels(10), getWindowHeight()-toPixels(20)),
                      Color(0.0, 0.0, 0.0),
                      Font("Arial", toPixels(12)));
  
  // Set projection/modelview matrices
  gl::setMatrices(cam);
  
  // Draw the rod and the normal of the bishop frame
  for(int i=0; i<y->numSegs(); i++) {
    Vec3c p0 = EtoC(y->cur().POS(i));
    Vec3c p1 = EtoC(y->cur().POS(i+1));
    gl::drawLine(p0, p1);
    gl::color(1.0, 1.0, 0.0);
    gl::lineWidth(1.0);
    Vec3c u = EtoC(y->cur().u[i]);
    gl::drawLine((p0+p1)/2.0, (p0+p1)/2.0+u*(p1-p0).length()*2.0);
  }
  
  m.apply();
  
  l->setDiffuse(Color::white());
  l->setAmbient(Color::white());
  l->setPosition(Vec3c(0.0, 50.0, 0.0));
  l->enable();
  
  diffuseProg.bind();
  for (int i=0; i<y->numCPs(); i++) {
    gl::pushModelView();
    gl::translate(EtoC(y->cur().POS(i)));
    spheredl->draw();
    gl::popModelView();
  }
  diffuseProg.unbind();
  
  yarnProg.bind();

  floorTex.enableAndBind();
  gl::draw(floor);
  floorTex.disable();
  
  yarnProg.unbind();
  
  // Draw yarn segments
  yarnProg.bind();
  yarnTex.enableAndBind();
  for (int i=0; i<y->numSegs(); i++) {
    gl::pushModelView();
    Vec3c v = EtoC(y->cur().vec(i).normalized());
    
    gl::translate(EtoC(y->cur().POS(i)));
    Quaternion<real> q(Vec3c(0.0, 1.0, 0.0), v);
    real angle = acos(std::max((real)-1.0, std::min((real)1.0, (q*Vec3c(-1.0, 0.0, 0.0)).dot(EtoC(y->cur().u[i])))));
    if ((q*Vec3c(-1.0, 0.0, 0.0)).dot(EtoC(y->cur().v(i))) > 0.0) angle = -angle;
    gl::rotate(Quaternion<real>(v, angle));
    gl::rotate(q);
    gl::rotate(Vec3c(0.0, y->cur().rot(i)*180.0/constants::pi, 0.0));
    gl::scale(1.0, y->cur().length(i), 1.0);
    cylinderdl->draw();
    gl::popModelView();
  }
  yarnTex.unbind();
  yarnProg.unbind();
  
  for (YarnEnergy* e : energies) {
    e->draw(c.timestep());
  }
#ifndef CONST_INTEGRATOR
  integrator->draw();
#endif // ifndef CONST_INTEGRATOR
  
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
  const size_t numPoints = std::stoi(line);
  
  VecXe yarnPos(3*numPoints);
  
  for (int i=0; i<3*numPoints; i++) {
    std::string line;
    std::getline(yarnFile, line);
    if(!line.empty()) {
      yarnPos(i) = std::stof(line);
    }
  }
  Vec3e u;
  for (int i=0; i<3; i++) {
    std::string line;
    std::getline(yarnFile, line);
    u(i) = std::stof(line);
  }
  assert((yarnPos.segment<3>(3) - yarnPos.segment<3>(0)).dot(u) < 5.0e-6);
  
  yarnFile.close();
  if (y) delete y;
  y = new Yarn(yarnPos, u);
}

void SimulatorApp::loadDefaultYarn(int numPoints) {
  if (y) delete y;
  
  eyePos = Vec3c(40.0, 10.0, 0.0);
  targetPos = Vec3c(0.0, 10.0, 0.0);
  cam.lookAt(eyePos, targetPos, Vec3c(0.0, 1.0, 0.0));
  
  Vec3e start = Vec3e(0.0, 1.6069, 0.0);
  Vec3e end   = Vec3e(0.0, 1.0, 0.0);
  
  Vec3e u     = (end-start).cross(Vec3e(0.0, 0.1, 0.0)).normalized();
  if (u.hasNaN() || u.norm() < 0.95) {
    u << 1.0, 0.0, 0.0;
  }
  
  VecXe yarnPos(3*numPoints);
  for(int i=0; i < numPoints; i++) {
    real t = ((real) i) / (real) (numPoints -1);
    yarnPos.segment<3>(3*i) = (1-t)*start + t*end;
  }
  
  real massPerPoint = 0.1;
  VecXe mass = VecXe::Constant(numPoints, massPerPoint);
  
  eyePos = Vec3c(5.0, 1.5, 0.0);
  targetPos = Vec3c(0.0, 1.5, 0.0);
  cam.lookAt(eyePos, targetPos, Vec3c(0.0, 1.0, 0.0));
  
  y = new Yarn(yarnPos, Vec3e(0.0, 0.0, 1.0), &mass, 1e7, 80.0);
}

void SimulatorApp::loadStdEnergies() {
  // Create Yarn Energies - Add in the order they are most likely to fail during evaluation
  assert(y && "Tried to load energies on a null yarn");
  for (YarnEnergy* e : energies) {
    delete e;
  }
  energies.clear();
  
  
  YarnEnergy* stretch = new Stretching(*y, Implicit);
  energies.push_back(stretch);
  
  YarnEnergy* bending = new Bending(*y, Implicit);
  energies.push_back(bending);
  
  YarnEnergy* twisting = new Twisting(*y, Explicit);
  energies.push_back(twisting);
  
  YarnEnergy* gravity = new Gravity(*y, Explicit, Vec3e(0.0, -9.8, 0.0));
  energies.push_back(gravity);
  
  mouseSpring = new MouseSpring(*y, Explicit, y->numCPs()-1, 100.0);
  energies.push_back(mouseSpring);
  
  YarnEnergy* intContact = new IntContact(*y, Implicit);
  energies.push_back(intContact);
  
  Spring* clamp1 = new Spring(*y, Implicit, 0, 500.0);
  clamp1->setClamp(y->rest().POS(0));
  Spring* clamp2 = new Spring(*y, Implicit, 1, 1000.0);
  clamp2->setClamp(y->rest().POS(1));
  energies.push_back(clamp1);
  energies.push_back(clamp2);
  
  if (integrator) delete integrator;
  integrator = new IMEXIntegrator(energies, *y);
}

void SimulatorApp::loadStdEnergiesAndConsts() {
  assert(y && "Tried to load energies and constraints on a null yarn");
  for (YarnEnergy* e : energies) {
    delete e;
  }
  for (YarnConstraint* c : constraints) {
    delete c;
  }
  energies.clear();
  constraints.clear();
  
  YarnEnergy* gravity = new Gravity(*y, Explicit, Vec3e(0.0, -9.8, 0.0));
  energies.push_back(gravity);
  
  Spring* clamp1 = new Spring(*y, Implicit, 0, 500.0);
  clamp1->setClamp(y->rest().POS(0));
  Spring* clamp2 = new Spring(*y, Implicit, 1, 1000.0);
  clamp2->setClamp(y->rest().POS(1));
  energies.push_back(clamp1);
  energies.push_back(clamp2);
  
  mouseSpring = new MouseSpring(*y, Explicit, y->numCPs()-1, 10.0);
  energies.push_back(mouseSpring);
  
  YarnConstraint* length = new Length(*y);
  constraints.push_back(length);
  
  if (cIntegrator) delete cIntegrator;
  cIntegrator = new ConstraintIntegrator(*y, energies, constraints);
}


CINDER_APP_NATIVE(SimulatorApp, RendererGl)
