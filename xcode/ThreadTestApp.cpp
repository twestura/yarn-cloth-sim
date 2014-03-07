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

#define N 20
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
  
  Yarn restYarn;
  Yarn* curYarn;
  Yarn* nextYarn;
  Workspace ws;
  Clock c;

  float twist = 0;
  
  // Interactive stuff
  bool isMouseDown = false;
  ci::Vec3f mousePosition = ci::Vec3f(0, 0, 0);
  bool isRotate = false;
  
};


void ThreadTestApp::setup()
{
  // Setup scene
  cam.setOrtho(cameraBounds[0], cameraBounds[1], cameraBounds[2], cameraBounds[3], 1, 100);
  cam.lookAt( eyePos, targetPos, ci::Vec3f( 0, 1, 0 ) );
  
  // Define Yarn
  for( int i = 0; i < NUM_VERTICES; i++ ) {
    CtrlPoint p;
    p.pos.x() = 0;
    p.pos.y() = (NUM_VERTICES-i)*20.0f/NUM_VERTICES-10;
    p.pos.z() = 0;
    p.vel.setZero();
    p.force.setZero();
    restYarn.points.push_back(p);
  }
  
  curYarn = new Yarn(restYarn);
  nextYarn = new Yarn(restYarn);
  
  for(int i=0; i<NUM_EDGES; i++) {
    Segment rest(restYarn.points[i], restYarn.points[i+1], Eigen::Vector3f(0, 0, 1));
    restYarn.segments.push_back(rest);
    
    Segment cur(curYarn->points[i], curYarn->points[i+1], Eigen::Vector3f(0, 0, 1));
    curYarn->segments.push_back(cur);
    
    Segment next(nextYarn->points[i], nextYarn->points[i+1], Eigen::Vector3f(0, 0, 1));
    nextYarn->segments.push_back(next);
  }
  
  for (int i=1; i<=restYarn.numIntCPs(); i++) {
    Segment& ePrev = restYarn.segments[i-1];
    Segment& eNext = restYarn.segments[i];
    
    Eigen::Vector3f curveBinorm = 2*ePrev.vec().cross(eNext.vec()) /
    (ePrev.length()*eNext.length() + ePrev.vec().dot(eNext.vec()));
    
    offVec<Eigen::Vector2f> matCurvature(-(i-1));
    matCurvature.push_back(Eigen::Vector2f(curveBinorm.dot(ePrev.m2()), -(curveBinorm.dot(ePrev.m1()))));
    matCurvature.push_back(Eigen::Vector2f(curveBinorm.dot(eNext.m2()), -(curveBinorm.dot(eNext.m1()))));
    ws.restMatCurvature.push_back(matCurvature);
  }
  
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
    default:;
  }
}

void ThreadTestApp::update()
{
  if (!running) return;
  
  /// 1. Compute forces on centerline
  
  c.suggestTimestep(1.0f/30.0f);
  Eigen::Vector3f mp;
  if (isMouseDown) {
    mp << mousePosition.x, mousePosition.y, mousePosition.z;
  } else {
    mp = curYarn->points[curYarn->numCPs()-1].pos;
  }
  Integrator::integrate(*curYarn, *nextYarn, restYarn, ws, c, mp);
  
  
  /// 3. Enforce inextensibility/clamped edge as a velocity filter
  
  Matrix<float, NUM_CONSTRAINTS, 3*NUM_VERTICES> Cgrad
    = Matrix<float, NUM_CONSTRAINTS, 3*NUM_VERTICES>::Zero();
  Matrix<float, NUM_CONSTRAINTS, 1> C = Matrix<float, NUM_CONSTRAINTS, 1>::Zero();
  float Cerror = 0;
  
  
  for(int i=0; i<NUM_EDGES; i++){
    Segment& si = nextYarn->segments[i];
    Segment& restsi = restYarn.segments[i];
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
  Eigen::Vector3f p0 = nextYarn->points[0].pos;
  Eigen::Vector3f p1 = nextYarn->points[1].pos;
  Eigen::Vector3f p0rest = restYarn.points[0].pos;
  Eigen::Vector3f p1rest = restYarn.points[1].pos;
  
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
      nextYarn->points[i].pos.x() += delta(3*i);
      nextYarn->points[i].pos.y() += delta(3*i+1);
      nextYarn->points[i].pos.z() += delta(3*i+2);
    }
    
    
    // Did we fix it?
    Cerror=0;
    Cgrad.setZero();
    
    for(int i=0; i<NUM_EDGES; i++){
      Segment& si = nextYarn->segments[i];
      Segment& restsi = restYarn.segments[i];
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
    Eigen::Vector3f& p0 = nextYarn->points[0].pos;
    Eigen::Vector3f& p1 = nextYarn->points[1].pos;
    Eigen::Vector3f& p0rest = restYarn.points[0].pos;
    Eigen::Vector3f& p1rest = restYarn.points[1].pos;
    
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
    nextYarn->points[i].vel = (nextYarn->points[i].pos - curYarn->points[i].pos) / c.timestep();
  }
  
  // Clamp top edge
  nextYarn->points[0].vel.setZero();
  nextYarn->points[1].vel.setZero();
  
  /// 4. (Collision detection and response would go here)
  
  /// 5. Update Bishop frame
  for(int i=1; i<NUM_EDGES; i++) {
    Segment& prevSeg = curYarn->segments[i];
    Segment& seg = nextYarn->segments[i];
    
    seg.parallelTransport(prevSeg);
  }
  
  /// 6. Update material frame rotation
  if (isRotate) {
    twist += ANG_VEL*H;
  }
  
  float delta = twist/(N); // Twist is constant and instant
  for(int i=1; i<nextYarn->numSegs(); i++) {
    nextYarn->segments[i].setRot(nextYarn->segments[i-1].getRot() + delta);
  }
  
  // Swap Yarns
  Yarn* temp = curYarn;
  curYarn = nextYarn;
  nextYarn = temp;

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
    gl::color( 0.4, -cos(curYarn->segments[i].getRot()), cos(curYarn->segments[i].getRot()) );
    gl::lineWidth(2);
    ci::Vec3f p0(curYarn->points[i].pos.x(), curYarn->points[i].pos.y(), curYarn->points[i].pos.z());
    ci::Vec3f p1(curYarn->points[i+1].pos.x(), curYarn->points[i+1].pos.y(), curYarn->points[i+1].pos.z());
    gl::drawLine(p0, p1);
    gl::color(1, 1, 0);
    gl::lineWidth(1);
    Eigen::Vector3f eu = curYarn->segments[i].getU();
    ci::Vec3f u(eu.x(), eu.y(), eu.z());
    gl::drawLine((p0+p1)/2, (p0+p1)/2+u);
  }
  
  if (isMouseDown) {
    gl::color(1, 0, 0);
    ci::Vec3f p(curYarn->points[NUM_VERTICES-1].pos.x(), curYarn->points[NUM_VERTICES-1].pos.y(), curYarn->points[NUM_VERTICES-1].pos.z());
    gl::drawLine(p, mousePosition);
  }
  
}




CINDER_APP_NATIVE( ThreadTestApp, RendererGl )
