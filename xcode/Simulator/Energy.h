//
//  Energy.h
//  Visualizer
//
//  Created by eschweickart on 3/7/14.
//
//

#ifndef Visualizer_Energy_h
#define Visualizer_Energy_h

#include "Eigen/Sparse"
#include "autodiff.h"
#include "Yarn.h"
#include "Constants.h"
#include <boost/timer.hpp>

typedef Eigen::Triplet<float> Triplet;
typedef Eigen::VectorXf VecXf;

enum EvalType {
  Implicit,
  Explicit,
  // Finite difference?
};

class YarnEnergy {
protected:
  const Yarn& y;
  EvalType et;
public:
  YarnEnergy(const Yarn& y, EvalType et) : y(y), et(et) { }
  virtual void eval(VecXf&, std::vector<Triplet>&, const VecXf&, Clock&) =0;
};

class Gravity : public YarnEnergy {
private:
  Vec3f dir;
public:
  Gravity(const Yarn& y, EvalType et, Vec3f dir) : YarnEnergy(y, et), dir(dir) {}
  void eval(VecXf& Fx, std::vector<Triplet>& GradFx, const VecXf& dqdot, Clock& c) {
    assert(et == Explicit && "Unsupported EvalType");
    for (int i=0; i<y.numCPs(); i++) {
      Fx(i*3)   -= dir.x() * c.timestep();
      Fx(i*3+1) -= dir.y() * c.timestep();
      Fx(i*3+2) -= dir.z() * c.timestep();
    }
  }
};

class Spring : public YarnEnergy {
protected:
  Vec3f clamp;
  size_t index;
  float stiffness;
public:
  Spring(const Yarn& y, EvalType et, size_t index, float stiffness) :
    YarnEnergy(y, et), index(index), stiffness(stiffness) {}
  
  void eval(VecXf& Fx, std::vector<Triplet>& GradFx, const VecXf& dqdot, Clock& c) {
    assert(et == Explicit && "Unsupported EvalType");
    size_t i = 3*index;
    Fx(i)   -= c.timestep() * stiffness * (clamp.x() - y.cur().points[index].pos.x());
    Fx(i+1) -= c.timestep() * stiffness * (clamp.y() - y.cur().points[index].pos.y());
    Fx(i+2) -= c.timestep() * stiffness * (clamp.z() - y.cur().points[index].pos.z());
  }
  
  void setClamp(Vec3f newClamp) { clamp = newClamp; }
};

class MouseSpring : public YarnEnergy {
private:
  bool mouseDown = false;
  bool mouseSet = false;
  Vec3f mouse;
  size_t index;
  float stiffness;
public:
  MouseSpring(const Yarn& y, EvalType et, size_t index, float stiffness) :
    YarnEnergy(y, et), index(index), stiffness(stiffness) {}
  
  void eval(VecXf& Fx, std::vector<Triplet>& GradFx, const VecXf& dqdot, Clock& c) {
    if (!mouseDown) return;
    assert(et == Explicit && "Unsupported EvalType");
    assert(mouseSet && "Set the mouse position each time you call eval()!");
    size_t i = 3*index;
    Fx(i)   -= c.timestep() * stiffness * (mouse.x() - y.cur().points[index].pos.x());
    Fx(i+1) -= c.timestep() * stiffness * (mouse.y() - y.cur().points[index].pos.y());
    Fx(i+2) -= c.timestep() * stiffness * (mouse.z() - y.cur().points[index].pos.z());
  }
  
  void setMouse(Vec3f newMouse, bool newDown) {
    mouse = newMouse;
    mouseDown = newDown;
    mouseSet = true;
  }
};


class Bending : public YarnEnergy {
private:
#define NUM_VARS 9
  typedef Eigen::Matrix<float, NUM_VARS, 1> Gradient;
  typedef Eigen::Matrix<float, NUM_VARS, NUM_VARS> Hessian;
  typedef DScalar2<float, NUM_VARS, Gradient, Hessian> DScalar;
  typedef DScalar::DVector3 DVector3;
  typedef DScalar::DVector2 DVector2;
  
  std::vector<Vec2f> restCurve;
  bool init = false;
  
public:
  Bending(const Yarn& y, EvalType et) : YarnEnergy(y, et) {
    // Init rest curve
    for (int i=1; i<y.numCPs()-1; i++) {
      const Segment& ePrev = y.rest().segments[i-1];
      const Segment& eNext = y.rest().segments[i];
      
      Vec3f curveBinorm = 2*ePrev.vec().cross(eNext.vec()) /
      (ePrev.length()*eNext.length() + ePrev.vec().dot(eNext.vec()));
      
      Vec2f restMatCurvePrev(curveBinorm.dot(ePrev.m2()), -(curveBinorm.dot(ePrev.m1())));
      Vec2f restMatCurveNext(curveBinorm.dot(eNext.m2()), -(curveBinorm.dot(eNext.m1())));
      Vec2f restMatCurve = 0.5*(restMatCurvePrev + restMatCurveNext);
      
      restCurve.push_back(restMatCurve);
    }
  }
  
  void eval(VecXf& Fx, std::vector<Triplet>& GradFx, const VecXf& dqdot, Clock& c) {
    DiffScalarBase::setVariableCount(NUM_VARS);
    float h = c.timestep();
    
    for (int i=1; i<y.numCPs()-1; i++) {
      const CtrlPoint& curPoint  = y.cur().points[i];
      const CtrlPoint& prevPoint = y.cur().points[i-1];
      const CtrlPoint& nextPoint = y.cur().points[i+1];
      const Segment&   prevSeg   = y.cur().segments[i-1];
      const Segment&   nextSeg   = y.cur().segments[i];
      
      // Redefine NUM_VARS if you change these
      DVector3 dPrevPoint(DScalar(0, prevPoint.pos.x() + h*(dqdot(3*(i-1))   + prevPoint.vel.x())),
                          DScalar(1, prevPoint.pos.y() + h*(dqdot(3*(i-1)+1) + prevPoint.vel.y())),
                          DScalar(2, prevPoint.pos.z() + h*(dqdot(3*(i-1)+2) + prevPoint.vel.z())));
      
      DVector3 dCurPoint(DScalar(3, curPoint.pos.x() + h*(dqdot(3*i)   + curPoint.vel.x())),
                         DScalar(4, curPoint.pos.y() + h*(dqdot(3*i+1) + curPoint.vel.y())),
                         DScalar(5, curPoint.pos.z() + h*(dqdot(3*i+2) + curPoint.vel.z())));
      
      DVector3 dNextPoint(DScalar(6, nextPoint.pos.x() + h*(dqdot(3*(i+1))   + nextPoint.vel.x())),
                          DScalar(7, nextPoint.pos.y() + h*(dqdot(3*(i+1)+1) + nextPoint.vel.y())),
                          DScalar(8, nextPoint.pos.z() + h*(dqdot(3*(i+1)+2) + nextPoint.vel.z())));
      
      DVector3 dPrevSeg = dCurPoint - dPrevPoint;
      DVector3 dNextSeg = dNextPoint - dCurPoint;
      assert(dPrevSeg.norm() != 0 && dNextSeg.norm() != 0 && "Edge length is 0");
      DVector3 dPrevSegN = dPrevSeg.normalized();
      DVector3 dNextSegN = dNextSeg.normalized();
      DScalar dotProd = dPrevSegN.dot(dNextSegN);
      assert(dotProd != -1 && "Segments are pointing in exactly opposite directions");
      
      DVector3 curveBinorm = (DScalar(2)*dPrevSegN.cross(dNextSegN))/(1+dotProd);
      
      Vec3f prevm1 = prevSeg.m1();
      Vec3f prevm2 = prevSeg.m2();
      Vec3f nextm1 = nextSeg.m1();
      Vec3f nextm2 = nextSeg.m2();
      
      DVector3 d1prev(DScalar(prevm1.x()), DScalar(prevm1.y()), DScalar(prevm1.z()));
      DVector3 d2prev(DScalar(prevm2.x()), DScalar(prevm2.y()), DScalar(prevm2.z()));
      DVector3 d1next(DScalar(nextm1.x()), DScalar(nextm1.y()), DScalar(nextm1.z()));
      DVector3 d2next(DScalar(nextm2.x()), DScalar(nextm2.y()), DScalar(nextm2.z()));
      
      DVector2 matCurvePrev(curveBinorm.dot(d2prev), -curveBinorm.dot(d1prev));
      DVector2 matCurveNext(curveBinorm.dot(d2next), -curveBinorm.dot(d1next));
      DVector2 matCurve = DScalar(0.5)*(matCurvePrev + matCurveNext);
      
      DVector2 restMatCurve(DScalar(restCurve[i-1].x()), DScalar(restCurve[i-1].y()));
      
      // TODO: bending matrix may not be I
      
      DScalar voronoiCell = 0.5*(dPrevSeg.norm()+dNextSeg.norm());
      DVector2 curveDiff = matCurve - restMatCurve;
      
      DScalar bendEnergy = 0.5*(1/voronoiCell)*curveDiff.dot(curveDiff);
      
      Gradient grad = bendEnergy.getGradient();
      Hessian hess = bendEnergy.getHessian();
      
      assert(et == Implicit && "Unsupported EvalType");
      for (int j=0; j<NUM_VARS; j++) {
        // TODO: mass matrix may not be I
        Fx(3*(i-1)+j) += h*grad(j);
        for (int k=0; k<NUM_VARS; k++) {
          float val = h*h*hess(j,k);
          CHECK_NAN(val);
          if (val != 0) {
            GradFx.push_back(Triplet(3*(i-1)+j, 3*(i-1)+k, val));
          }
        }
      }
    }
    
  }
#undef NUM_VARS
};

class Stretching : public YarnEnergy {
private:
#define NUM_VARS 6
  typedef Eigen::Matrix<float, NUM_VARS, 1> Gradient;
  typedef Eigen::Matrix<float, NUM_VARS, NUM_VARS> Hessian;
  typedef DScalar2<float, NUM_VARS, Gradient, Hessian> DScalar;
  typedef DScalar::DVector3 DVector3;
  
  // TODO: find a good value for this
  float youngsModulus = 1e7;
  float xArea = constants::pi * constants::radius * constants::radius;
  float stretchScalar = xArea * youngsModulus;
  
public:
  Stretching(const Yarn& y, EvalType et) : YarnEnergy(y, et) { }
  void eval(VecXf& Fx, std::vector<Triplet>& GradFx, const VecXf& dqdot, Clock& c) {
    DiffScalarBase::setVariableCount(NUM_VARS);
    float h = c.timestep();
    
    for (int i=0; i<y.numSegs(); i++) {
      const Segment& seg = y.cur().segments[i];
      float restSegLength = y.rest().segments[i].length();
      const CtrlPoint& prevPoint = seg.getFirst();
      const CtrlPoint& nextPoint = seg.getSecond();
      
      // Redefine NUM_VARS if you change these
      DVector3 dPrevPoint(DScalar(0, prevPoint.pos.x() + h*(dqdot(3*i)   + prevPoint.vel.x())),
                          DScalar(1, prevPoint.pos.y() + h*(dqdot(3*i+1) + prevPoint.vel.y())),
                          DScalar(2, prevPoint.pos.z() + h*(dqdot(3*i+2) + prevPoint.vel.z())));
      
      DVector3 dNextPoint(DScalar(3, nextPoint.pos.x() + h*(dqdot(3*(i+1))   + nextPoint.vel.x())),
                          DScalar(4, nextPoint.pos.y() + h*(dqdot(3*(i+1)+1) + nextPoint.vel.y())),
                          DScalar(5, nextPoint.pos.z() + h*(dqdot(3*(i+1)+2) + nextPoint.vel.z())));
      
      DScalar axialStrain = (dNextPoint - dPrevPoint).norm()/restSegLength - 1;
      
      DScalar stretchEnergy = (0.5 * stretchScalar * restSegLength) * axialStrain * axialStrain;
      
      Gradient grad = stretchEnergy.getGradient();
      Hessian hess = stretchEnergy.getHessian();
      
      assert(et == Implicit && "Unsuported EvalType");
      for (int j=0; j<NUM_VARS; j++) {
        Fx(3*i+j) += h*grad(j);
        for (int k=0; k<NUM_VARS; k++) {
          float val = h*h*hess(j,k);
          if (val != 0) {
            GradFx.push_back(Triplet(3*i+j, 3*i+k, val));
          }
        }
      }
      
    }
  }
#undef NUM_VARS
};

#endif
