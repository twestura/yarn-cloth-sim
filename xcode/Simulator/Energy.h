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
    for (int i=2; i<y.numCPs(); i++) { // TODO: that "2" should not be hard-coded
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
public:
  Spring(const Yarn& y, EvalType et, size_t index) : YarnEnergy(y, et), index(index) {}
  
  void eval(VecXf& Fx, std::vector<Triplet>& GradFx, const VecXf& dqdot, Clock& c) {
    assert(et == Explicit && "Unsupported EvalType");
    size_t i = 3*index;
    Fx(i)   -= c.timestep() * (clamp.x() - y.cur().points[index].pos.x());
    Fx(i+1) -= c.timestep() * (clamp.y() - y.cur().points[index].pos.y());
    Fx(i+2) -= c.timestep() * (clamp.z() - y.cur().points[index].pos.z());
  }
  
  void setClamp(Vec3f newClamp) { clamp = newClamp; }
};

class MouseSpring : public YarnEnergy {
private:
  bool mouseDown = false;
  bool mouseSet = false;
  Vec3f mouse;
  size_t index;
public:
  MouseSpring(const Yarn& y, EvalType et, size_t index) : YarnEnergy(y, et), index(index) {}
  
  void eval(VecXf& Fx, std::vector<Triplet>& GradFx, const VecXf& dqdot, Clock& c) {
    if (!mouseDown) return;
    assert(et == Explicit && "Unsupported EvalType");
    assert(mouseSet && "Set the mouse position each time you call eval()!");
    size_t i = 3*index;
    Fx(i)   -= c.timestep() * (mouse.x() - y.cur().points[index].pos.x());
    Fx(i+1) -= c.timestep() * (mouse.y() - y.cur().points[index].pos.y());
    Fx(i+2) -= c.timestep() * (mouse.z() - y.cur().points[index].pos.z());
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
  typedef DScalar2<float, Gradient, Hessian> DScalar;
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
      DVector3 dPrevPoint(DScalar(0, prevPoint.pos.x() + h*dqdot(3*(i-1))   + h*prevPoint.vel.x()),
                          DScalar(1, prevPoint.pos.y() + h*dqdot(3*(i-1)+1) + h*prevPoint.vel.y()),
                          DScalar(2, prevPoint.pos.z() + h*dqdot(3*(i-1)+2) + h*prevPoint.vel.z()));
      
      DVector3 dCurPoint(DScalar(3, curPoint.pos.x() + h*dqdot(3*i)   + h*curPoint.vel.x()),
                         DScalar(4, curPoint.pos.y() + h*dqdot(3*i+1) + h*curPoint.vel.y()),
                         DScalar(5, curPoint.pos.z() + h*dqdot(3*i+2) + h*curPoint.vel.z()));
      
      DVector3 dNextPoint(DScalar(6, nextPoint.pos.x() + h*dqdot(3*(i+1))   + h*nextPoint.vel.x()),
                          DScalar(7, nextPoint.pos.y() + h*dqdot(3*(i+1)+1) + h*nextPoint.vel.y()),
                          DScalar(8, nextPoint.pos.z() + h*dqdot(3*(i+1)+2) + h*nextPoint.vel.z()));
      
      DVector3 dPrevSeg = dCurPoint - dPrevPoint;
      DVector3 dNextSeg = dNextPoint - dCurPoint;
      assert(dPrevSeg.norm() != 0 && dNextSeg.norm() != 0 && "Edge length is 0");
      DVector3 dPrevSegN = dPrevSeg.normalized();
      DVector3 dNextSegN = dNextSeg.normalized();
      DScalar dotProd = dPrevSegN.dot(dNextSegN);
      assert(dotProd != -1 && "Segments are pointing in exactly opposite directions");
      
      DVector3 curveBinorm = (DScalar(2)*dPrevSegN.cross(dNextSegN))/(1+dotProd);
      
      DVector3 d1prev(DScalar(prevSeg.m1().x()),
                      DScalar(prevSeg.m1().y()),
                      DScalar(prevSeg.m1().z()));
      
      DVector3 d2prev(DScalar(prevSeg.m2().x()),
                      DScalar(prevSeg.m2().y()),
                      DScalar(prevSeg.m2().z()));
      
      DVector3 d1next(DScalar(nextSeg.m1().x()),
                      DScalar(nextSeg.m1().y()),
                      DScalar(nextSeg.m1().z()));
      
      DVector3 d2next(DScalar(nextSeg.m2().x()),
                      DScalar(nextSeg.m2().y()),
                      DScalar(nextSeg.m2().z()));
      
      DVector2 matCurvePrev(curveBinorm.dot(d2prev), -curveBinorm.dot(d1prev));
      DVector2 matCurveNext(curveBinorm.dot(d2next), -curveBinorm.dot(d1next));
      DVector2 matCurve = DScalar(0.5)*(matCurvePrev + matCurveNext);
      
      DVector2 restMatCurve(DScalar(restCurve[i-1].x()), DScalar(restCurve[i-1].y()));
      
      // TODO: bending matrix may not be I
      
      DScalar voronoiCell = 0.5*(dPrevSeg.norm()+dNextSeg.norm());
      
      DScalar bendEnergy = 0.5*(1/voronoiCell)*
      (matCurve - restMatCurve).dot(matCurve - restMatCurve);
      
      Gradient grad = bendEnergy.getGradient();
      Hessian hess = bendEnergy.getHessian();
      
      assert(et == Implicit && "Unsupported EvalType");
      for (int j=0; j<NUM_VARS; j++) {
        // TODO: mass matrix may not be I
        Fx(3*(i-1)+j) -= h*grad(j);
        for (int k=0; k<NUM_VARS; k++) {
          float val = -h*h*hess(j,k);
          CHECK_NAN(val);
          if (val != 0) {
            GradFx.push_back(Triplet(3*(i-1)+j, 3*(i-1)+k, val));
          }
        }
      }
    }
  }
};

class Stretching : public YarnEnergy {
public:
  void eval(VecXf& Fx, std::vector<Triplet>& GradFx, const VecXf& dqdot, Clock& c) {
    // TODO
  }
};

#endif
