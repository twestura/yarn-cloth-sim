//
//  RodSoundApp.cpp
//  Visualizer
//
//  Created by eschweickart on 6/9/14.
//
//

#include "RodSoundApp.h"
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
#include "Rod.h"
#include "FrameExporter.h"
#include "ExIntegrator.h"
#include "YarnBuilder.h"
#include "BEMSolver.h"

using namespace ci;
using namespace ci::app;

class RodSoundApp : public AppNative {
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
  
  // Set to false to pause the simulation
  bool running = true;
  
  // Camera for the scene, along with its position and orientation
  CameraPersp cam;
  Vec3c eyePos   = Vec3c(50.0, 0.0, 0.0);
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
  std::vector<RodEnergy*> energies;
  std::vector<RodConstraint*> constraints;
  MouseSpring* mouseSpring;
  
  real twist = 0.0;
  real rodTwist = 0.0;
  int numRodTwists = 0;
  
  // Interactive stuff
  bool isMouseDown = false;
  Vec3c mousePosition;
  bool isRotate = false;
  
  // Sound stuff
#define SimulationLength 3.0 // in seconds
  const static std::size_t BufferSize = (std::size_t)(SampleRate * SimulationLength);
  double sampleBuffer[BufferSize];
  
  real tAtLastDraw = 0.0;
  bool stopNow = false;
  Vec3e ear2Pos = Vec3e(28.0, 10.0, 28.0);
  double sampleBuffer2[BufferSize];
  std::size_t curSample = 0;
  double sampleBuffer3[BufferSize];
  
  std::size_t multiSample = 13;
  FrameExporter fe;
};

void RodSoundApp::setup()
{
//  std::cout << solveBEM(constants::radius) << "\n\n";
//  std::cout << "Expected:\n" << -constants::pi * constants::radius * constants::radius * Mat2e::Identity() << "\n\n";
  
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
  loadDefaultRod(50);
  // loadRodFile("");
  loadStdEnergies();
  
  PROFILER_START("Total");
}

void RodSoundApp::mouseDown(MouseEvent event)
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

void RodSoundApp::mouseDrag(MouseEvent event)
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

void RodSoundApp::mouseUp(MouseEvent event)
{
  if (!running) return;
  isMouseDown = false;
  
}

void RodSoundApp::mouseWheel(MouseEvent event) {
  real scroll = event.getWheelIncrement();
  eyePos += (targetPos - eyePos).normalized() * scroll;
  cam.lookAt(eyePos, targetPos);
}

void RodSoundApp::keyDown(KeyEvent event)
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
      
    case event.KEY_s:
      stopNow = true;
      break;
      
    default:;
  }
}

void RodSoundApp::resize() {
  cam.setAspectRatio(getWindowAspectRatio());
}

void RodSoundApp::update()
{
  if (!running) return;
  
  if (curSample % 5000 == 0 && curSample != 0 && c.getTicks() % multiSample == 0) {
    std::cout << curSample << " / " << BufferSize << " (" << (curSample*100.0)/BufferSize << "%)\n";
    PROFILER_PRINT_ELAPSED();
    PROFILER_RESET_ALL();
    std::cout << "\n";
  }
  
  if (curSample >= BufferSize || stopNow) { // We're done!
    sampleBuffer[0] = 0.0; // To prevent the click of forces suddenly being applied
    double max = 0;
    for (int i=0; i<BufferSize; i++) {
      max = std::max(max, std::fabs(sampleBuffer[i]));
    }
    std::cout << "Max1: " << max << "\n";
    uint16_t buffer[BufferSize];
    for (int i=0; i<BufferSize; i++) {
      buffer[i] = toSample(sampleBuffer[i], max);
    }
    writeWAVData((constants::ResultPath+"result.wav").data(), buffer,
                 curSample * sizeof(uint16_t), SampleRate, 1);
    
    sampleBuffer2[0] = 0.0;
    max = 0;
    for (int i=0; i<BufferSize; i++) {
      max = std::max(max, std::fabs(sampleBuffer2[i]));
    }
    std::cout << "Max2: " << max << "\n";
    for (int i=0; i<BufferSize; i++) {
      buffer[i] = toSample(sampleBuffer2[i], max);
    }
    writeWAVData((constants::ResultPath+"result2.wav").data(), buffer,
                 curSample * sizeof(uint16_t), SampleRate, 1);
    
    sampleBuffer3[0] = 0.0;
    max = 0;
    for (int i=0; i<BufferSize; i++) {
      max = std::max(max, std::fabs(sampleBuffer3[i]));
    }
    std::cout << "Max3: " << max << "\n";
    for (int i=0; i<BufferSize; i++) {
      buffer[i] = toSample(sampleBuffer3[i], max);
    }
    writeWAVData((constants::ResultPath+"result3.wav").data(), buffer,
                 curSample * sizeof(uint16_t), SampleRate, 1);
    
    fe.writeMPEG("result");
    std::cout << "Total simulation time: " << app::getElapsedSeconds() << "\n"; // FIXME: This is inaccurate
    
    running = false;
    return;
  }
  
  PROFILER_START("Update");
  
  c.suggestTimestep(1.0 / (real) SampleRate / multiSample);
  // FIXME: Normally the frame exporter would suggest a timestep, but this interferes with the audio
  // recording, as it assumes all timesteps are 1/SampleRate. However, any error the frame exporter
  // experiences is small since 1/60 >> 1/SampleRate.
  // fe.suggestTimestep(c);
  
  Vec3e mp;
  if (isMouseDown) mp << mousePosition.x, mousePosition.y, mousePosition.z;
  mouseSpring->setMouse(mp, isMouseDown);
  
  if (!integrator->integrate(c)) throw;
  
  /// Update Bishop frame
  r->next().updateReferenceFrames(r->cur());
  
  // Sound Calculations
  if (c.getTicks() % multiSample == 0) {
    real sample = 0;
    real sample2 = 0;
    real avgX = 0;
    VecXe jerkVec = r->next().dVel - r->cur().dVel;
    for (int i=1; i<r->numCPs()-1; i++) {
      avgX += r->next().VEL(i).x();
      
      // Calculate jerk
      Vec3e jerk = jerkVec.segment<3>(3*i);
      // Project jerk to transverse plane
      Vec3e tPlaneNormal = (r->next().edge(i-1) + r->next().edge(i)).normalized();
      jerk = jerk - jerk.dot(tPlaneNormal) * tPlaneNormal; // Vector rejection of jerk from tPlaneNormal
      
      /*
      real m0 = r->restVoronoiLength(i)*constants::pi*r->radius()*r->radius()*constants::rhoAir;
      // Rotation to align system so that the cylinder is coaxial with the z-axis
      Eigen::Quaternion<real> q = Eigen::Quaternion<real>::FromTwoVectors(tPlaneNormal, Vec3e(0, 0, 1));
      Vec3e rotJerk = q * jerk;
      rotJerk = rotJerk.cwiseProduct(Vec3e(2.0*m0, 2.0*m0, m0));
      
      // Calculate sample contribution
      Vec3e earVec = CtoE(eyePos) - r->next().points[i].pos;
      sample +=  (q * earVec).dot(rotJerk) / (4.0 * constants::pi * constants::cAir * earVec.dot(earVec));
      
      earVec = ear2Pos - r->next().points[i].pos;
      sample2 +=  (q * earVec).dot(rotJerk) / (4.0 * constants::pi * constants::cAir * earVec.dot(earVec));
      */
       
      
      Vec3e earVec = CtoE(eyePos) - r->next().POS(i);
      // Calculate sample contribution
      sample += r->getCS()[i].calcSample(earVec, jerk);
    
      earVec = ear2Pos - r->next().POS(i);
      sample2 += r->getCS()[i].calcSample(earVec, jerk);
    }
    avgX = avgX/(r->numCPs()-2);
    sampleBuffer[curSample] = sample;
    sampleBuffer2[curSample] = sample2;
    
    sampleBuffer3[curSample] = r->next().VEL(r->numCPs()/2).x() - avgX;
    
    curSample++;
  }
  
  // Swap Rods
  r->swapRods();

  c.increment();
  PROFILER_STOP("Update");
}

void RodSoundApp::draw() {
  while (running &&
//         app::getElapsedSeconds() - tAtLastDraw < 1.0/app::getFrameRate() &&
         fe.nextTimestep(c) > 1.0 / (real) SampleRate) {
    update();
  }
  tAtLastDraw = app::getElapsedSeconds();
  
  PROFILER_START("Draw");
  
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
  integrator->draw();
 
  fe.record(c);
  
  PROFILER_STOP("Draw");
}

void RodSoundApp::loadRodFile(std::string filename) {
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
    if (!line.empty()) {
      rodPos(i) = std::stof(line);
    }
  }
  Vec3e u;
  for (int i=0; i<3; i++) {
    std::string line;
    std::getline(rodFile, line);
    u(i) = std::stof(line);
  }
  assert((rodPos.segment<3>(3) - rodPos.segment<3>(0)).dot(u) < 5.0e-6);
  
  rodFile.close();
  if (r) delete r;
  r = new Rod(rodPos, u);
}

void RodSoundApp::loadDefaultRod(int numPoints) {
  if (r) delete r;
  
  Vec3e start = Vec3e(0.0, 1.9144, 0.0); // Vec3e(0.5, 1.0, 0.5);
  Vec3e end   = Vec3e(0.0, 1.0, 0.0); // start + (Vec3e(-1.0, -1.0, -1.0).normalized() * 0.6069);
  // 1ft. = 0.3048m
  // 2ft. = 0.6069m
  // 3ft. = 0.9144m

  Vec3e u     = (end-start).cross(Vec3e(0.0, 0.1, 0.0)).normalized();
  if (u.hasNaN() || u.norm() < 0.95) {
    u << 1.0, 0.0, 0.0;
  }
  
  VecXe rodPos(3*numPoints);
  for(int i=0; i < numPoints; i++) {
    real t = ((real) i) / (real) (numPoints -1);
    rodPos.segment<3>(3*i) = (1-t)*start + t*end;
  }
  
  eyePos = Vec3c(5.0, 1.5, 0.0);
  targetPos = Vec3c(0.0, 1.5, 0.0);
  cam.lookAt(eyePos, targetPos, Vec3c(0.0, 1.0, 0.0));
  
  // FIXME: assumes circular cylindrical rod of default radius
  real totalMass = constants::radius * constants::radius * constants::pi * (start - end).norm() * constants::rhoRod;
  real massPerPoint = totalMass / numPoints;
  
  VecXe mass = VecXe::Constant(numPoints, massPerPoint);
  r = new Rod(rodPos, u, &mass);
  
  real l = (start - end).norm();
  real kappa = sqrt(r->youngsModulus * r->getCS()[0].areaMoment()(0, 0) /
                    (constants::rhoRod * r->getCS()[0].area() * l * l * l * l));
  std::cout << "kappa: " << kappa << "\n";
  real h = l / (numPoints-1);
  real k = 1.0 / (44100*multiSample);
  real mu = kappa * k / (h * h);
  std::cout << "mu: " << mu << "\n";
  
  real fmax = asin(2.0 * mu) / (constants::pi * k);
  std::cout << "fmax: " << fmax << "\n";
  /*
  real cb = std::pow(constants::youngsModulus * r->getCS()[0].areaMoment()(0, 0) / constants::rhoRod / r->getCS()[0].area(), 0.25);
  real cbHigh = cb * 141.4;
  real cbLow = cb * 4.47;
  real dx = (start - end).norm() / (numPoints - 1.0);
  std::cout << "Min timestep: " << dx / cbHigh << "\n";
   */
}

void RodSoundApp::loadStdEnergies() {
  // Create Rod Energies - Add in the order they are most likely to fail during evaluation
  assert(r && "Tried to load energies on a null rod");
  for (RodEnergy* e : energies) {
    delete e;
  }
  energies.clear();
  
  
  RodEnergy* stretch = new Stretching(*r, Explicit);
  energies.push_back(stretch);
  // OR
  //RodConstraint* length = new Length(*r);
  //constraints.push_back(length);
  
  RodEnergy* bending = new Bending(*r, Explicit);
  energies.push_back(bending);
  
  RodEnergy* fembending = new FEMBending(*r, Explicit);
//  energies.push_back(fembending);
  
  RodEnergy* twisting = new Twisting(*r, Explicit);
//  energies.push_back(twisting);
  
  Vec3e dir(0.0, -9.8, 0.0);
  RodEnergy* gravity = new Gravity(*r, Explicit, dir);
//  energies.push_back(gravity);
  
  mouseSpring = new MouseSpring(*r, Explicit, r->numCPs()-1, 100.0);
  energies.push_back(mouseSpring);
  
  Vec3e floorNormal(0.0, 1.0, 0.0);
  Vec3e floorOrigin = Vec3e::Zero();
  RodEnergy* floor = new PlaneContact(*r, Explicit, floorNormal, floorOrigin, 5000.0);
//  energies.push_back(floor);
  
  
  Vec3e imp1dir(1.0e-5, 0.0, 0.0);
  Vec3e imp2dir(-1.0e-5, 0.0, 0.0);
  Vec3e imp3dir(0.0, 0.0, 1.0e-10);
  RodEnergy* imp1 = new Impulse(*r, Explicit, c, 0.2, 0.201, imp1dir, 3);
  RodEnergy* imp2 = new Impulse(*r, Explicit, c, 0.2, 0.201, imp2dir, r->numCPs()-2);
  RodEnergy* imp3 = new Impulse(*r, Explicit, c, 0.2, 0.201, imp3dir, r->numCPs()-1);
  energies.push_back(imp1);  // energies.push_back(imp2); // energies.push_back(imp3);
  
  RodEnergy* spr1 = new Spring(*r, Explicit, 0, 1e5);
  Vec3e spr1clamp = r->rest().POS(0) + Vec3e(0.0, 0.1, 0.0);
  static_cast<Spring*>(spr1)->setClamp(spr1clamp);
  RodEnergy* spr2 = new Spring(*r, Explicit, r->numCPs()-1, 1e5);
  Vec3e spr2clamp = r->rest().POS(r->numCPs()-1);
  static_cast<Spring*>(spr2)->setClamp(spr2clamp);
  energies.push_back(spr1); energies.push_back(spr2);
  
  /*
  
  RodEnergy* intContact = new IntContact(*r, Explicit);
  energies.push_back(intContact);
  
   */
  
  if (integrator) delete integrator;
  integrator = new ExIntegrator(*r, energies, &constraints);
}


CINDER_APP_NATIVE(RodSoundApp, RendererGl)
