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
#include "Eigen/Dense"

#include "Util.h"
#include "Yarn.h"
#include "Clock.h"
#include "Integrator.h"

using namespace ci;
using namespace ci::app;
using namespace Eigen;

#define N 40
#define NUM_EDGES (N+1)
#define NUM_VERTICES (N+2)
#define H 0.01f // Timestep.
#define GRAVITY ci::Vec3f(0, -40, 0) // The direction/intensity of gravity.
#define ANG_VEL (10*3.14159) // Constant angular velocity of the material frame.
#define NUM_CONSTRAINTS (NUM_EDGES + 6) // Number of constraints to solve

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
  
  // Set to false to pause the simulation
  bool running = true;
  
  // Camera for the scene, along with its position and orientation
  CameraOrtho cam = CameraOrtho();
  Vec4f cameraBounds = Vec4f(-15, 15, -15, 15);
  ci::Vec3f eyePos = ci::Vec3f( 50, 0, 0 );
  ci::Vec3f targetPos = ci::Vec3f( 0, 0, 0 );
  
  Yarn* y;
  Clock c;
  Integrator* integrator;
  std::vector<YarnEnergy*> energies;
  MouseSpring* mouseSpring;

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
  
  std::vector<Eigen::Vector3f> yarnPoints;
  // Define Yarn
  for(int i=0; i < NUM_VERTICES; i++) {
    Eigen::Vector3f p(0, (NUM_VERTICES-i)*20.0f/NUM_VERTICES-10, 0);
    yarnPoints.push_back(p);
  }
  
  y = new Yarn(yarnPoints, Eigen::Vector3f(0, 0, 1));
  
  // Create Yarn Energies
  YarnEnergy* gravity = new Gravity(*y, Explicit, Eigen::Vector3f(0, -9.8, 0));
  energies.push_back(gravity);
  
  mouseSpring = new MouseSpring(*y, Explicit, NUM_VERTICES-1, 10);
  energies.push_back(mouseSpring);
  
  YarnEnergy* bending = new Bending(*y, Implicit);
  energies.push_back(bending);
  
  YarnEnergy* twisting = new Twisting(*y, Explicit);
  energies.push_back(twisting);
  
  YarnEnergy* intContact = new IntContact(*y, Explicit);
  energies.push_back(intContact);
  
  Spring* clamp1 = new Spring(*y, Explicit, 0, 1000);
  clamp1->setClamp(y->rest().points[0].pos);
  Spring* clamp2 = new Spring(*y, Explicit, 1, 1000);
  clamp2->setClamp(y->rest().points[1].pos);
  energies.push_back(clamp1);
  energies.push_back(clamp2);
  
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
      p0.y() += 1;
      break;
      case 's':
      p0.y() -= 1;
      break;
      
    default:;
  }
}

void ThreadTestApp::update()
{
  if (!running) return;

  // TODO: adaptive timesteps
  c.suggestTimestep(1.0f/30.0f);
  
  Eigen::Vector3f mp;
  if (isMouseDown) {
    mp << mousePosition.x, mousePosition.y, mousePosition.z;
  } else {
    mp = y->cur().points[y->numCPs()-1].pos;
  }
  mouseSpring->setMouse(mp, isMouseDown);
  
  integrator->integrate(*y, c);
  
  /// 3. Enforce inextensibility/clamped edge as a velocity filter
#ifdef MANIFOLD_PROJECTION
  
  Matrix<float, NUM_CONSTRAINTS, 3*NUM_VERTICES> Cgrad
    = Matrix<float, NUM_CONSTRAINTS, 3*NUM_VERTICES>::Zero();
  Matrix<float, NUM_CONSTRAINTS, 1> C = Matrix<float, NUM_CONSTRAINTS, 1>::Zero();
  float Cerror = 0;
  
  
  for(int i=0; i<NUM_EDGES; i++){
    Segment& si = y->next().segments[i];
    const Segment& restsi = y->rest().segments[i];
    float restLength = restsi.length();
    Eigen::Vector3f siv = si.vec();
    
    // The Kirsch paper recommends this constraint, which is supposedly more rubust
    // that the one proposed in the Bergou paper.
    C(i) = (siv.dot(siv)/(2*restLength)) - (restLength / 2);
    
    Cgrad(i, 3*i) = -siv.x()/restLength;
    Cgrad(i, 3*i+1) = -siv.y()/restLength;
    Cgrad(i, 3*i+2) = -siv.z()/restLength;
    Cgrad(i, 3*(i+1)) = siv.x()/restLength;
    Cgrad(i, 3*(i+1)+1) = siv.y()/restLength;
    Cgrad(i, 3*(i+1)+2) = siv.z()/restLength;
    
  }
  
  // Boundary edge constraints
  Eigen::Vector3f& p0 = y->next().points[0].pos;
  Eigen::Vector3f& p1 = y->next().points[1].pos;
  const Eigen::Vector3f& p0rest = y->rest().points[0].pos;
  const Eigen::Vector3f& p1rest = y->rest().points[1].pos;
  
  C(NUM_CONSTRAINTS-6) = p0rest.x() - p0.x();
  C(NUM_CONSTRAINTS-5) = p0rest.y() - p0.y();
  C(NUM_CONSTRAINTS-4) = p0rest.z() - p0.z();
  C(NUM_CONSTRAINTS-3) = p1rest.x() - p1.x();
  C(NUM_CONSTRAINTS-2) = p1rest.y() - p1.y();
  C(NUM_CONSTRAINTS-1) = p1rest.z() - p1.z();
  
  for (int i=0; i<6; i++) {
    Cgrad(NUM_CONSTRAINTS-6+i, i) = -1;
  }
  
  Cerror = fmaxf(std::abs(C.minCoeff()), std::abs(C.maxCoeff()));
  
  int iter = 0;
  while (Cerror > 0.0001) {
    iter++;
    
    //Timer timer = Timer(true);
    
    // Solving delta = -M^-1 * G^T * [G * M^-1 * G^T]^-1 * lambda
    // where G=Cgrad, M is the identity, and lambda is the vector of Lagrange multipliers
    Matrix<float, 3*NUM_VERTICES, NUM_CONSTRAINTS> Cgradtrans = Cgrad.transpose();
    Matrix<float, NUM_CONSTRAINTS, NUM_CONSTRAINTS> Cgradsym = Cgrad * Cgradtrans;
    Matrix<float, NUM_CONSTRAINTS, 1> lambda = Cgradsym.ldlt().solve(-C);
    Matrix<float, 3*NUM_VERTICES, 1> delta = Cgradtrans * lambda;
    
    //cout << timer.getSeconds() << "\n";
    
    // Update where are control points should be
    for(int i=0; i<NUM_VERTICES; i++) {
      y->next().points[i].pos.x() += delta(3*i);
      y->next().points[i].pos.y() += delta(3*i+1);
      y->next().points[i].pos.z() += delta(3*i+2);
    }
    
    
    // Did we fix it?
    Cerror=0;
    Cgrad.setZero();
    
    for(int i=0; i<NUM_EDGES; i++){
      Segment& si = y->next().segments[i];
      const Segment& restsi = y->rest().segments[i];
      float restLength = restsi.length();
      Eigen::Vector3f siv = si.vec();
      
      // The Kirsch paper recommends this constraint, which is supposedly more rubust
      // that the one proposed in the Bergou paper.
      C(i) = (siv.dot(siv)/(2*restLength)) - (restLength / 2);
      
      Cgrad(i, 3*i) = -siv.x()/restLength;
      Cgrad(i, 3*i+1) = -siv.y()/restLength;
      Cgrad(i, 3*i+2) = -siv.z()/restLength;
      Cgrad(i, 3*(i+1)) = siv.x()/restLength;
      Cgrad(i, 3*(i+1)+1) = siv.y()/restLength;
      Cgrad(i, 3*(i+1)+2) = siv.z()/restLength;
      
    }
    
    // Boundary edge constraints
    Eigen::Vector3f& p0 = y->next().points[0].pos;
    Eigen::Vector3f& p1 = y->next().points[1].pos;
    const Eigen::Vector3f& p0rest = y->rest().points[0].pos;
    const Eigen::Vector3f& p1rest = y->rest().points[1].pos;
    
    C(NUM_CONSTRAINTS-6) = p0rest.x() - p0.x();
    C(NUM_CONSTRAINTS-5) = p0rest.y() - p0.y();
    C(NUM_CONSTRAINTS-4) = p0rest.z() - p0.z();
    C(NUM_CONSTRAINTS-3) = p1rest.x() - p1.x();
    C(NUM_CONSTRAINTS-2) = p1rest.y() - p1.y();
    C(NUM_CONSTRAINTS-1) = p1rest.z() - p1.z();
    
    for (int i=0; i<6; i++) {
      Cgrad(NUM_CONSTRAINTS-6+i, i) = -1;
    }
    
    Cerror = fmaxf(std::abs(C.minCoeff()), std::abs(C.maxCoeff()));
    
    
    if (iter >= 100) { // Things will probably just go downhill from here
      std::cout << "Iterator failed to converge; aborting...\n";
      running = false;
      break;
    }
  }
  
  for(int i=0; i<NUM_VERTICES; i++) {
    // Update velocities once we've converged
    y->next().points[i].vel = (y->next().points[i].pos - y->cur().points[i].pos) / c.timestep();
  }
  
#endif // ifdef MANIFOLD_PROJECTION
  
  /// 4. (Collision detection and response would go here)
  
  /// 5. Update Bishop frame
  for(int i=1; i<NUM_EDGES; i++) {
    Segment& prevSeg = y->cur().segments[i];
    Segment& seg = y->next().segments[i];
    
    if (i < NUM_EDGES-1) {
      Segment& refSeg = y->next().segments[i+1];
      seg.parallelTransport(prevSeg, refSeg);
    } else {
      seg.parallelTransport(prevSeg);
    }
  }
  
  /// 6. Update material frame rotation
  if (isRotate) {
    twist += ANG_VEL*H;
  }
  
  float delta = twist/(N); // Twist is constant and instant
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
  
  // Set projection/modelview matrices
  gl::setMatrices(cam);
  
  // Draw the rod and the normal of the bishop frame
  for(int i=1; i<NUM_EDGES; i++) {
    ci::Color c(0.4, -cos(y->cur().segments[i].getRot()), cos(y->cur().segments[i].getRot()));
    c *= (y->cur().segments[i].getFirst().pos.z() + 5) / 10;
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
    ci::Vec3f p(y->cur().points[NUM_VERTICES-1].pos.x(), y->cur().points[NUM_VERTICES-1].pos.y(), y->cur().points[NUM_VERTICES-1].pos.z());
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
    p = s.eval(t, true);
    ci::Vec3f v1(p.x(), p.y(), p.z());
    plTest.push_back(v1);
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




CINDER_APP_NATIVE( ThreadTestApp, RendererGl )
