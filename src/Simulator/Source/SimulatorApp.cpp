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
#include "cinder/ImageIo.h"

#include "Resources.h"
#include "Util.h"
#include "Rod.h"
#include "Clock.h"
#include "Integrator.h"
#include "IMEXIntegrator.h"
#include "ConstraintIntegrator.h"
#include "YarnBuilder.h"

#include "TypeDefs.h"

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
  
  void loadRodFile(std::string filename);
  void loadDefaultRod(int numPoints);
  void loadStdEnergies();
  void loadStdEnergiesAndConsts();
  
  // Set to false to pause the simulation
  bool running = true;
  
  // Camera for the scene, along with its position and orientation
  CameraPersp cam;
  Vec3c eyePos = Vec3c(50.0, 0.0, 0.0);
  Vec3c targetPos = Vec3c(0.0, 0.0, 0.0);
  
  // Rendering stuff
  gl::GlslProg rodProg;
  gl::GlslProg diffuseProg;
  gl::Texture rodTex;
  gl::Texture floorTex;
  gl::DisplayList* spheredl;
  gl::DisplayList* cylinderdl;
  gl::Material m = gl::Material(Color(0.3, 0.3, 0.3), Color(0.8, 0.9, 0.9));
  gl::Light* l;
  TriMesh floor;
  
  Rod* r = nullptr;
  Clock c;
  Integrator* integrator = nullptr;
  ConstraintIntegrator* cIntegrator = nullptr;
  std::vector<RodEnergy*> energies;
  MouseSpring* mouseSpring;
  std::vector<RodConstraint*> constraints;
  
  real twist = 0.0;
  real rodTwist = 0.0;
  int numRodTwists = 0;
  
  // Interactive stuff
  bool isMouseDown = false;
  Vec3c mousePosition;
  bool isRotate = false;

  bool isSetup = false;
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
    rodTex = loadImage(loadResource(RES_SIM_YARN_TEX));
  } catch (ImageIoException e) {
    std::cerr << "Error loading textures: " << e.what();
    exit(1);
  }
  
  // Load and compile shaders
  try {
    diffuseProg = gl::GlslProg(loadResource(RES_SIM_VERT_GLSL), loadResource(RES_SIM_FRAG_GLSL));
    rodProg = gl::GlslProg(loadResource(RES_SIM_VERT_TEX_GLSL), loadResource(RES_SIM_FRAG_TEX_GLSL));
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
  
  // Load the rod
  loadDefaultRod(42);
  // loadRodFile("");
#ifndef CONST_INTEGRATOR
  loadStdEnergies();
#else
  loadStdEnergiesAndConsts();
#endif // ifndef CONST_INTEGRATOR

  isSetup = true;
}

void SimulatorApp::mouseDown(MouseEvent event)
{
  if (event.isRight()) {
    // Set targetPos to the ControlPoint we just clicked
    if (!r) return;
    Vec2i mouse = event.getPos();
    Vec2i windowSize = getWindowSize();
    Ray ray = cam.generateRay((real)mouse.x/windowSize.x,
                            1.0 - (real)mouse.y/windowSize.y,
                            getWindowAspectRatio());
    real tmin = INFINITY;
    bool any = false;
    for (int i=0; i<r->numCPs(); i++) { // A bit slow, but beats keeping a KD-Tree updated
      Sphere s(EtoC(r->cur().POS(i)), constants::radius * 1.5);
      float t;
      if (s.intersect(ray, &t) && t < tmin) {
        any = true;
        tmin = t;
      }
    }
    if (!any) return;
    targetPos = ray.calcPosition(tmin);
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
      loadRodFile("");
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
      
      
    case event.KEY_t:
      for(int i=1; i<r->numEdges(); i++) {
        std::cout << "edge " << i << " twist: " << (r->cur().rot(i) - r->cur().rot(i-1) + r->cur().refTwist(i))
        << " (" << r->cur().rot(i) << " + " << r->cur().refTwist(i) <<
        " = " << r->cur().rot(i) + r->cur().refTwist(i) << ")\n";
      }
      std::cout << "numRodTwists: " << numRodTwists << "\n";
      break;
      
    default:;
  }
}

void SimulatorApp::resize() {
  cam.setAspectRatio(getWindowAspectRatio());
}

void SimulatorApp::update()
{
  if (!running || !isSetup) return;
  
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
  r->next().updateReferenceFrames(r->cur());
  
#ifndef CONST_INTEGRATOR
  /// Update material frame rotation
  if (isRotate) {
    twist += 2.0*constants::pi*c.timestep();
  }
  
  Vec3e uRef = parallelTransport(r->next().edge(0), r->next().edge(r->numEdges()-1), r->next().u[0]);
  real cosTwist = r->next().u[r->numEdges()-1].dot(uRef.normalized());
  real oldTwist = rodTwist;
  if (cosTwist >= 1.0) { // Avoid values like 1.0000000012 that introduce NaNs
    rodTwist = 0.0;
  } else if (cosTwist <= -1.0) {
    rodTwist = constants::pi;
  } else {
    rodTwist = acos(cosTwist);
  }
  // Flip the sign if necessary
  if (r->next().v(r->numEdges()-1).dot(uRef) > 0.0) {
    rodTwist = -rodTwist;
  }
  real diff = rodTwist - oldTwist;
  if (diff < -constants::pi) {
    numRodTwists += 1;
  } else if (diff > constants::pi) {
    numRodTwists -= 1;
  }
  r->next().rot(r->numEdges()-1) = twist - rodTwist;
  if (!static_cast<IMEXIntegrator*>(integrator)->setRotations()) {
    std::cout << "rotations failed";
  }
#endif // ifndef CONST_INTEGRATOR
  
  // Swap Rods
  r->swapRods();

  c.increment();
}

void SimulatorApp::draw() {

	if (!isSetup) {
		return;
	}

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
  for(int i=0; i<r->numEdges(); i++) {
    Vec3c p0 = EtoC(r->cur().POS(i));
    Vec3c p1 = EtoC(r->cur().POS(i+1));
    gl::drawLine(p0, p1);
    gl::color(1.0, 1.0, 0.0);
    gl::lineWidth(1.0);
    Vec3c u = EtoC(r->cur().u[i]);
    gl::drawLine((p0+p1)/2.0, (p0+p1)/2.0+u*(p1-p0).length()*2.0);
  }
  
  m.apply();
  
  l->setDiffuse(Color::white());
  l->setAmbient(Color::white());
  l->setPosition(Vec3c(0.0, 50.0, 0.0));
  l->enable();
  
  diffuseProg.bind();
  for (int i=0; i<r->numCPs(); i++) {
    gl::pushModelView();
    gl::translate(EtoC(r->cur().POS(i)));
    spheredl->draw();
    gl::popModelView();
  }
  diffuseProg.unbind();
  
  rodProg.bind();

  floorTex.enableAndBind();
  gl::draw(floor);
  floorTex.disable();
  
  rodProg.unbind();
  
  // Draw rod edges
  rodProg.bind();
  rodTex.enableAndBind();
  for (int i=0; i<r->numEdges(); i++) {
    gl::pushModelView();
    Vec3c v = EtoC(r->cur().edge(i).normalized());
    
    gl::translate(EtoC(r->cur().POS(i)));
    Quaternion<real> q(Vec3c(0.0, 1.0, 0.0), v);
    real angle = acos(std::max((real)-1.0, std::min((real)1.0, (q*Vec3c(-1.0, 0.0, 0.0)).dot(EtoC(r->cur().u[i])))));
    if ((q*Vec3c(-1.0, 0.0, 0.0)).dot(EtoC(r->cur().v(i))) > 0.0) angle = -angle;
    gl::rotate(Quaternion<real>(v, angle));
    gl::rotate(q);
    gl::rotate(Vec3c(0.0, r->cur().rot(i)*180.0/constants::pi, 0.0));
    gl::scale(1.0, r->cur().edgeLength(i), 1.0);
    cylinderdl->draw();
    gl::popModelView();
  }
  rodTex.unbind();
  rodProg.unbind();
  
  for (RodEnergy* e : energies) {
    e->draw(c.timestep());
  }
#ifndef CONST_INTEGRATOR
  integrator->draw();
#endif // ifndef CONST_INTEGRATOR
  
}

void SimulatorApp::loadRodFile(std::string filename) {
  if (filename.empty()) filename = getOpenFilePath().string();
  if (filename.empty()) return;
  
  std::ifstream rodFile(filename);
  
  if (!rodFile.is_open()) {
    std::cerr << filename << " failed to open!\n";
    return;
  }
  
  std::string line;
  std::getline(rodFile, line);
  const std::size_t numPoints = std::stoi(line);
  
  VecXe rodPos(3*numPoints);
  
  for (int i=0; i<3*numPoints; i++) {
    std::string line;
    std::getline(rodFile, line);
    if(!line.empty()) {
      rodPos(i) = std::stof(line);
    }
  }
  Vec3e u;
  for (int i=0; i<3; i++) {
    std::string line;
    std::getline(rodFile, line);
    u(i) = std::stof(line);
  }
  //assert((rodPos.segment<3>(3) - rodPos.segment<3>(0)).dot(u) < 5.0e-6);
  
  rodFile.close();
  if (r) delete r;
  r = new Rod(rodPos, u);
}

void SimulatorApp::loadDefaultRod(int numPoints) {
  if (r) delete r;
  
  eyePos = Vec3c(40.0, 10.0, 0.0);
  targetPos = Vec3c(0.0, 10.0, 0.0);
  cam.lookAt(eyePos, targetPos, Vec3c(0.0, 1.0, 0.0));
  
  Vec3e start = Vec3e(0.0, 1.6069, 0.0);
  Vec3e end   = Vec3e(0.0, 1.0, 0.0);
  
  Vec3e u     = (end-start).cross(Vec3e(0.0, 0.1, 0.0)).normalized();
  if (u.hasNaN() || u.norm() < 0.95) {
    u << 1.0, 0.0, 0.0;
  }
  
  VecXe rodPos(3*numPoints);
  for(int i=0; i < numPoints; i++) {
    real t = ((real) i) / (real) (numPoints -1);
    rodPos.segment<3>(3*i) = (1-t)*start + t*end;
  }
  
  real massPerPoint = 0.1;
  VecXe mass = VecXe::Constant(numPoints, massPerPoint);
  
  eyePos = Vec3c(5.0, 1.5, 0.0);
  targetPos = Vec3c(0.0, 1.5, 0.0);
  cam.lookAt(eyePos, targetPos, Vec3c(0.0, 1.0, 0.0));
  
  r = new Rod(rodPos, Vec3e(0.0, 0.0, 1.0), &mass, 1e7, 80.0);
}

void SimulatorApp::loadStdEnergies() {
  // Create Rod Energies - Add in the order they are most likely to fail during evaluation
  assert(r && "Tried to load energies on a null rod");
  for (RodEnergy* e : energies) {
    delete e;
  }
  energies.clear();
  
  
  RodEnergy* stretch = new Stretching(*r, Implicit);
  energies.push_back(stretch);
  
  RodEnergy* bending = new Bending(*r, Implicit);
  energies.push_back(bending);
  
  RodEnergy* twisting = new Twisting(*r, Explicit);
  energies.push_back(twisting);
  
  RodEnergy* gravity = new Gravity(*r, Explicit, Vec3e(0.0, -9.8, 0.0));
  energies.push_back(gravity);
  
  mouseSpring = new MouseSpring(*r, Explicit, r->numCPs()-1, 100.0);
  energies.push_back(mouseSpring);
  
  RodEnergy* intContact = new IntContact(*r, Implicit);
  energies.push_back(intContact);
  
  Spring* clamp1 = new Spring(*r, Implicit, 0, 500.0);
  clamp1->setClamp(r->rest().POS(0));
  Spring* clamp2 = new Spring(*r, Implicit, 1, 1000.0);
  clamp2->setClamp(r->rest().POS(1));
  energies.push_back(clamp1);
  energies.push_back(clamp2);
  
  if (integrator) delete integrator;
  integrator = new IMEXIntegrator(energies, *r);
}

void SimulatorApp::loadStdEnergiesAndConsts() {
  assert(r && "Tried to load energies and constraints on a null rod");
  for (RodEnergy* e : energies) {
    delete e;
  }
  for (RodConstraint* c : constraints) {
    delete c;
  }
  energies.clear();
  constraints.clear();
  
  RodEnergy* gravity = new Gravity(*r, Explicit, Vec3e(0.0, -9.8, 0.0));
  energies.push_back(gravity);
  
  Spring* clamp1 = new Spring(*r, Implicit, 0, 500.0);
  clamp1->setClamp(r->rest().POS(0));
  Spring* clamp2 = new Spring(*r, Implicit, 1, 1000.0);
  clamp2->setClamp(r->rest().POS(1));
  energies.push_back(clamp1);
  energies.push_back(clamp2);
  
  mouseSpring = new MouseSpring(*r, Explicit, r->numCPs()-1, 10.0);
  energies.push_back(mouseSpring);
  
  RodConstraint* length = new Length(*r);
  constraints.push_back(length);
  
  if (cIntegrator) delete cIntegrator;
  cIntegrator = new ConstraintIntegrator(*r, energies, constraints);
}


CINDER_APP_NATIVE(SimulatorApp, RendererGl)
