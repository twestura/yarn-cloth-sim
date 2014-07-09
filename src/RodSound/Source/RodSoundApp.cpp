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
#include "Yarn.h"
#include "FrameExporter.h"
#include "ExIntegrator.h"
#include "YarnBuilder.h"

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
  
  void loadYarnFile(std::string filename);
  void loadDefaultYarn(int numPoints);
  void loadStdEnergies();
  
  // Set to false to pause the simulation
  bool running = true;
  
  // Camera for the scene, along with its position and orientation
  CameraPersp cam;
  Vec3c eyePos   = Vec3c(50.0, 0.0, 0.0);
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
  std::vector<YarnEnergy*> energies;
  MouseSpring* mouseSpring;
  
  Spring* testSpring1;
  Spring* testSpring2;
  Vec3e testSpring1Clamp = Vec3e(10.0, 15.0, 5.0);
  Vec3e testSpring2Clamp = Vec3e(-10.0, 15.0, 5.0);
  
  real twist = 0.0;
  real yarnTwist = 0.0;
  int numYarnTwists = 0;
  
  // Interactive stuff
  bool isMouseDown = false;
  Vec3c mousePosition;
  bool isRotate = false;
  
  // Sound stuff
  constexpr static real SimulationLength = 3.0; // in seconds
  constexpr static size_t BufferSize = (size_t)(SampleRate * SimulationLength);
  double sampleBuffer[BufferSize];
  
  real tAtLastDraw = 0.0;
  bool stopNow = false;
  Vec3e ear2Pos = Vec3e(28.0, 10.0, 28.0);
  double sampleBuffer2[BufferSize];
  size_t curSample = 0;
  
  FrameExporter fe;
};

void RodSoundApp::setup()
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
  loadStdEnergies();
  
  PROFILER_START("Total");
}

void RodSoundApp::mouseDown(MouseEvent event)
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
    for (const CtrlPoint& p : y->cur().points) { // A bit slow, but beats keeping a KD-Tree updated
      Sphere s(EtoC(p.pos), constants::radius * 1.5);
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
      /*
      testSpring1Clamp.z() -= 1.0;
      testSpring2Clamp.z() -= 1.0;
      testSpring1->setClamp(testSpring1Clamp);
      testSpring2->setClamp(testSpring2Clamp);
       */
      stopNow = true;
      break;
      
    case event.KEY_t:
      for(int i=1; i<y->numSegs(); i++) {
        const Segment& s = y->cur().segments[i];
        const Segment& sPrev = y->cur().segments[i-1];
        std::cout << "seg " << i << " twist: " << (s.getRot() - sPrev.getRot() + s.getRefTwist())
        << " (" << y->cur().segments[i].getRot() << " + " << y->cur().segments[i].getRefTwist() <<
        " = " << y->cur().segments[i].getRot() + y->cur().segments[i].getRefTwist() << ")\n";
      }
      std::cout << "numYarnTwists: " << numYarnTwists << "\n";
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
  
  if (curSample % 5000 == 0 && curSample != 0) {
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
    std::cout << "Max: " << max << "\n";
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
    for (int i=0; i<BufferSize; i++) {
      buffer[i] = toSample(sampleBuffer2[i], max);
    }
    writeWAVData((constants::ResultPath+"result2.wav").data(), buffer,
                 curSample * sizeof(uint16_t), SampleRate, 1);
    
    fe.writeMPEG("result");
    std::cout << "Total simulation time: " << app::getElapsedSeconds() << "\n"; // FIXME: This is inaccurate
    
    running = false;
    return;
  }
  
  PROFILER_START("Update");
  
  c.suggestTimestep(1.0 / (real) SampleRate);
  // FIXME: Normally the frame exporter would suggest a timestep, but this interferes with the audio
  // recording, as it assumes all timesteps are 1/SampleRate. However, any error the frame exporter
  // experiences is small since 1/60 >> 1/SampleRate.
  // fe.suggestTimestep(c);
  
  Vec3e mp;
  if (isMouseDown) mp << mousePosition.x, mousePosition.y, mousePosition.z;
  mouseSpring->setMouse(mp, isMouseDown);
  
  if (!integrator->integrate(c)) throw;
  
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
  
  /*
  // Update material frame rotation
  if (isRotate) {
    twist += 2.0*constants::pi*c.timestep();
  }
  const Segment& sFirst = y->next().segments[0];
  Segment& sLast = y->next().segments[y->numSegs()-1];
  Vec3e uRef = Segment::parallelTransport(sFirst.vec(), sLast.vec(), sFirst.getU());
  real cosTwist = sLast.getU().normalized().dot(uRef.normalized());
  real oldTwist = yarnTwist;
  if (cosTwist >= 1.0) { // Avoid values like 1.0000000012 that introduce NaNs
    yarnTwist = 0.0;
  } else if (cosTwist <= -1.0) {
    yarnTwist = constants::pi;
  } else {
    yarnTwist = acos(cosTwist);
  }
  // Flip the sign if necessary
  if (sLast.v().normalized().dot(uRef) > 0.0) {
    yarnTwist = -yarnTwist;
  }
  real diff = yarnTwist - oldTwist;
  if (diff < -constants::pi) {
    numYarnTwists += 1;
  } else if (diff > constants::pi) {
    numYarnTwists -= 1;
  }
  sLast.setRot(twist - (yarnTwist)); // + 2*constants::pi*numYarnTwists));
  if (!integrator->setRotations()) {
    std::cout << "rotations failed";
  }
   */
  
  // Sound Calculations
  if (c.getTicks() % 1 == 0) {
    real sample = 0;
    real sample2 = 0;
    for (int i=1; i<y->numCPs()-1; i++) {
      // Calculate jerk
      Vec3e jerk = y->next().points[i].accel - y->cur().points[i].accel;
      // Project it to transverse plane
      Vec3e tPlaneNormal = (y->next().segments[i-1].vec() + y->next().segments[i].vec()).normalized();
      jerk = jerk - jerk.dot(tPlaneNormal) * tPlaneNormal; // Vector rejection of jerk from tPlaneNormal
      
      /*
      real m0 = y->restVoronoiLength(i)*constants::pi*y->radius()*y->radius()*constants::rhoAir;
      // Rotation to align system so that the cylinder is coaxial with the z-axis
      Eigen::Quaternion<real> q = Eigen::Quaternion<real>::FromTwoVectors(tPlaneNormal, Vec3e(0, 0, 1));
      Vec3e rotJerk = q * jerk;
      rotJerk = rotJerk.cwiseProduct(Vec3e(2.0*m0, 2.0*m0, m0));
      
      // Calculate sample contribution
      Vec3e earVec = CtoE(eyePos) - y->next().points[i].pos;
      sample +=  (q * earVec).dot(rotJerk) / (4.0 * constants::pi * constants::cAir * earVec.dot(earVec));
      
      earVec = ear2Pos - y->next().points[i].pos;
      sample2 +=  (q * earVec).dot(rotJerk) / (4.0 * constants::pi * constants::cAir * earVec.dot(earVec));
      */
       
      
      Vec3e earVec = CtoE(eyePos) - y->next().points[i].pos;
      // Calculate sample contribution
      sample += (constants::rhoAir*y->radius()*y->radius()*y->radius() /
                 (2.0*constants::cAir*earVec.norm()*earVec.norm())) * (earVec.dot(jerk));
    
      earVec = ear2Pos - y->next().points[i].pos;
      sample2 += (constants::rhoAir*y->radius()*y->radius()*y->radius() /
                  (2.0*constants::cAir*earVec.norm()*earVec.norm())) * (earVec.dot(jerk));
    }
    sampleBuffer[curSample] = sample;
    sampleBuffer2[curSample] = sample2;
    curSample++;
  }
  
  // Swap Yarns
  y->swapYarns();
  
  /*
  // Update twists in new yarn
  for (int i=0; i<y->numSegs(); i++) {
    y->next().segments[i].updateTwists(y->cur().segments[i]);
  }
   */

  c.increment();
  PROFILER_STOP("Update");
}

void RodSoundApp::draw() {
  while (running &&
         // app::getElapsedSeconds() - tAtLastDraw < 1.0/app::getFrameRate() &&
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
  for(int i=0; i<y->numSegs(); i++) {
    Vec3c p0 = EtoC(y->cur().points[i].pos);
    Vec3c p1 = EtoC(y->cur().points[i+1].pos);
    gl::drawLine(p0, p1);
    gl::color(1.0, 1.0, 0.0);
    gl::lineWidth(1.0);
    Vec3c u = EtoC(y->cur().segments[i].getU());
    gl::drawLine((p0+p1)/2.0, (p0+p1)/2.0+u);
  }
  
  m.apply();
  
  l->setDiffuse(Color::white());
  l->setAmbient(Color::white());
  l->setPosition(Vec3c(0.0, 50.0, 0.0));
  l->enable();
  
  diffuseProg.bind();
  for (int i=0; i<y->numCPs(); i++) {
    gl::pushModelView();
    gl::translate(EtoC(y->cur().points[i].pos));
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
    const Segment& s = y->cur().segments[i];
    Vec3c v = EtoC(s.vec().normalized());
    
    gl::translate(EtoC(s.getFirst().pos));
    Quaternion<real> q(Vec3c(0.0, 1.0, 0.0), v);
    real angle = acosf(std::max((real)-1.0,
                                std::min((real)1.0,(q*Vec3c(-1.0, 0.0, 0.0)).dot(EtoC(s.getU())))));
    if ((q*Vec3c(-1.0, 0.0, 0.0)).dot(EtoC(s.v())) > 0.0) angle = -angle;
    gl::rotate(Quatf(v, angle));
    gl::rotate(q);
    gl::rotate(Vec3c(0.0, s.getRot()*180.0/constants::pi, 0.0));
    gl::scale(1.0, s.length(), 1.0);
    cylinderdl->draw();
    gl::popModelView();
  }
  yarnTex.unbind();
  yarnProg.unbind();

  for (YarnEnergy* e : energies) {
    e->draw(c.timestep());
  }
  integrator->draw();
 
  fe.record(c);
  
  PROFILER_STOP("Draw");
}

void RodSoundApp::loadYarnFile(std::string filename) {
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
  
  std::vector<Vec3e> yarnPoints;
  yarnPoints.reserve(numPoints);
  
  for (int j=0; j<numPoints; j++) {
    Vec3e p;
    for (int i=0; i<3; i++) {
      std::string line;
      std::getline(yarnFile, line);
      if(!line.empty()) {
        p(i) = std::stof(line);
      }
    }
    yarnPoints.push_back(p);
  }
  Vec3e u;
  for (int i=0; i<3; i++) {
    std::string line;
    std::getline(yarnFile, line);
    u(i) = std::stof(line);
  }
  assert((yarnPoints[1] - yarnPoints[0]).dot(u) < 5.0e-6);
  
  yarnFile.close();
  if (y) delete y;
  y = new Yarn(yarnPoints, u);
}

void RodSoundApp::loadDefaultYarn(int numPoints) {
  if (y) delete y;
  
  Vec3e start = Vec3e(0.0, 20.0, 0.0); // Vec3e(-5.0, 4.0, 3.0);
  Vec3e end   = Vec3e(0.0, 1.0, 0.0); // Vec3e(5.0, 3.0, -3.0);

  Vec3e u     = (end-start).cross(Vec3e(0.0, 0.1, 0.0)).normalized();
  if (u.hasNaN() || u.norm() < 0.95) {
    u << 1.0, 0.0, 0.0;
  }
  
  std::vector<Vec3e> yarnPoints;
  for(int i=0; i < numPoints; i++) {
    real t = ((real) i) / (real) (numPoints -1);
    Vec3e p = (1-t)*start + t*end;
    yarnPoints.push_back(p);
  }
  
  eyePos = Vec3c(40.0, 10.0, 0.0);
  targetPos = Vec3c(0.0, 10.0, 0.0);
  cam.lookAt(eyePos, targetPos, Vec3c(0.0, 1.0, 0.0));
  
//  VecXe mass = VecXe::Constant(numPoints, 0.2);
  y = new Yarn(yarnPoints, u); //, &mass);
}

void RodSoundApp::loadStdEnergies() {
  // Create Yarn Energies - Add in the order they are most likely to fail during evaluation
  assert(y && "Tried to load evergies on a null yarn");
  for (YarnEnergy* e : energies) {
    delete e;
  }
  energies.clear();
  
  
  YarnEnergy* stretch = new Stretching(*y, Explicit);
  energies.push_back(stretch);
  
  YarnEnergy* bending = new Bending(*y, Explicit);
  energies.push_back(bending);
  
  YarnEnergy* twisting = new Twisting(*y, Explicit);
//  energies.push_back(twisting);
  
  YarnEnergy* gravity = new Gravity(*y, Explicit, Vec3e(0.0, -9.8, 0.0));
//  energies.push_back(gravity);
  
  mouseSpring = new MouseSpring(*y, Explicit, y->numCPs()-1, 100.0);
  energies.push_back(mouseSpring);
  
  YarnEnergy* floor = new PlaneContact(*y, Explicit, Vec3e(0.0, 1.0, 0.0), Vec3e::Zero(), 5000.0);
//  energies.push_back(floor);
  
  
  YarnEnergy* imp1 = new Impulse(*y, Explicit, c, 0.2, 0.21, Vec3e(0.0, 0.0, -500.0), 0);
  YarnEnergy* imp2 = new Impulse(*y, Explicit, c, 0.2, 0.21, Vec3e(0.0, 0.0, 500.0), y->numCPs()/2);
  YarnEnergy* imp3 = new Impulse(*y, Explicit, c, 0.2, 0.21, Vec3e(0.0, 0.0, -500.0), y->numCPs()-1);
  energies.push_back(imp1); energies.push_back(imp2); energies.push_back(imp3);
  
  /*
  
  YarnEnergy* intContact = new IntContact(*y, Explicit);
  energies.push_back(intContact);
  
  Spring* clamp1 = new Spring(*y, Implicit, 0, 500.0);
  clamp1->setClamp(y->rest().points[0].pos);
//  clamp1->setClamp(y->rest().points[0].pos + Vec3e(0.0, 6.0, 2.0));
  Spring* clamp2 = new Spring(*y, Implicit, 1, 1000.0);
  clamp2->setClamp(y->rest().points[1].pos);
//  Spring* clamp2 = new Spring(*y, Implicit, 14, 500.0);
//  clamp2->setClamp(y->rest().points[14].pos + Vec3e(0.0, -6.0, 2.0));
  Spring* clamp3 = new Spring(*y, Implicit, 28, 500.0);
  clamp3->setClamp(y->rest().points[28].pos + Vec3e(0.0, 6.0, -2.0));
  Spring* clamp4 = new Spring(*y, Implicit, 42, 500.0);
  clamp4->setClamp(y->rest().points[42].pos + Vec3e(0.0, -6.0, -2.0));
  energies.push_back(clamp1);
  energies.push_back(clamp2);
//  energies.push_back(clamp3);
//  energies.push_back(clamp4);
  
  testSpring1 = new Spring(*y, Explicit, 2*y->numCPs()/3, 50.0);
  testSpring1->setClamp(testSpring1Clamp);
  testSpring2 = new Spring(*y, Explicit, y->numCPs()-1, 50.0);
  testSpring2->setClamp(testSpring2Clamp);
//  energies.push_back(testSpring1);
//  energies.push_back(testSpring2);
   */
  
  if (integrator) delete integrator;
  integrator = new ExIntegrator(*y, energies);
}


CINDER_APP_NATIVE(RodSoundApp, RendererGl)
