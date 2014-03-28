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

#define ENABLE_AUTODIFF

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
      Fx.block<3,1>(3*i, 0) -= dir * c.timestep();
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
    Fx.block<3,1>(3*index, 0) -= (clamp - y.cur().points[index].pos) * stiffness * c.timestep();
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
    Fx.block<3,1>(3*index, 0) -= (mouse - y.cur().points[index].pos) * stiffness * c.timestep();
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
  typedef Eigen::Vector2f Vec2f;
  typedef Eigen::Matrix<float, NUM_VARS, 1> Gradient;
  typedef Eigen::Matrix<float, NUM_VARS, NUM_VARS> Hessian;
  typedef DScalar2<float, NUM_VARS, Gradient, Hessian> DScalar;
  typedef DScalar::DVector3 DVector3;
  typedef DScalar::DVector2 DVector2;
  
  std::vector<Vec2f> restCurve;
  std::vector<float> voronoiCell;
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
      
      voronoiCell.push_back(0.5*(ePrev.length()+eNext.length()));
    }
  }
  
  void eval(VecXf& Fx, std::vector<Triplet>& GradFx, const VecXf& dqdot, Clock& c) {
    Profiler::start("Bend Eval");
    DiffScalarBase::setVariableCount(NUM_VARS);
    float h = c.timestep();
    
    for (int i=1; i<y.numCPs()-1; i++) {
      const CtrlPoint& curPoint  = y.cur().points[i];
      const CtrlPoint& prevPoint = y.cur().points[i-1];
      const CtrlPoint& nextPoint = y.cur().points[i+1];
      const Segment&   prevSeg   = y.cur().segments[i-1];
      const Segment&   nextSeg   = y.cur().segments[i];
      
#ifdef ENABLE_AUTODIFF
      // Redefine NUM_VARS if you change these
      DVector3 dvPrevPoint(DScalar(0, prevPoint.pos.x() + h*(dqdot(3*(i-1))   + prevPoint.vel.x())),
                           DScalar(1, prevPoint.pos.y() + h*(dqdot(3*(i-1)+1) + prevPoint.vel.y())),
                           DScalar(2, prevPoint.pos.z() + h*(dqdot(3*(i-1)+2) + prevPoint.vel.z())));
      
      DVector3 dvCurPoint(DScalar(3, curPoint.pos.x() + h*(dqdot(3*i)   + curPoint.vel.x())),
                          DScalar(4, curPoint.pos.y() + h*(dqdot(3*i+1) + curPoint.vel.y())),
                          DScalar(5, curPoint.pos.z() + h*(dqdot(3*i+2) + curPoint.vel.z())));
      
      DVector3 dvNextPoint(DScalar(6, nextPoint.pos.x() + h*(dqdot(3*(i+1))   + nextPoint.vel.x())),
                           DScalar(7, nextPoint.pos.y() + h*(dqdot(3*(i+1)+1) + nextPoint.vel.y())),
                           DScalar(8, nextPoint.pos.z() + h*(dqdot(3*(i+1)+2) + nextPoint.vel.z())));
            
      DVector3 dvPrevSeg = dvCurPoint - dvPrevPoint;
      DVector3 dvNextSeg = dvNextPoint - dvCurPoint;
      assert(dvPrevSeg.norm() != 0 && dvNextSeg.norm() != 0 && "Edge length is 0");
      
      DVector3 dvPrevSegN = dvPrevSeg.normalized();
      DVector3 dvNextSegN = dvNextSeg.normalized();
      DScalar dotProd = dvPrevSegN.dot(dvNextSegN);
      assert(dotProd != -1 && "Segments are pointing in exactly opposite directions");
      
      DVector3 dvcurveBinorm = (DScalar(2)*dvPrevSegN.cross(dvNextSegN))/(1+dotProd);
      
      Vec3f prevm1 = prevSeg.m1();
      Vec3f prevm2 = prevSeg.m2();
      Vec3f nextm1 = nextSeg.m1();
      Vec3f nextm2 = nextSeg.m2();
      
      DVector3 d1prev(DScalar(prevm1.x()), DScalar(prevm1.y()), DScalar(prevm1.z()));
      DVector3 d2prev(DScalar(prevm2.x()), DScalar(prevm2.y()), DScalar(prevm2.z()));
      DVector3 d1next(DScalar(nextm1.x()), DScalar(nextm1.y()), DScalar(nextm1.z()));
      DVector3 d2next(DScalar(nextm2.x()), DScalar(nextm2.y()), DScalar(nextm2.z()));
      
      DVector2 matCurvePrev(dvcurveBinorm.dot(d2prev), -dvcurveBinorm.dot(d1prev));
      DVector2 matCurveNext(dvcurveBinorm.dot(d2next), -dvcurveBinorm.dot(d1next));
      DVector2 dvmatCurve = DScalar(0.5)*(matCurvePrev + matCurveNext);
      
      DVector2 restMatCurve(DScalar(restCurve[i-1].x()), DScalar(restCurve[i-1].y()));
      
      // TODO: bending matrix may not be I
      
      DVector2 curveDiff = dvmatCurve - restMatCurve;
      
      DScalar bendEnergy = (0.5/voronoiCell[i-1])*curveDiff.dot(curveDiff);
      
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

      
      /*
      if (i == y.numCPs()-2) {
        std::cout << hess << "\n\n";
      }
       */
      
#else // ifdef ENABLE_AUTODIFF
    
      Vec3f dPrevPoint = prevPoint.pos + h*(dqdot.block<3, 1>(3*(i-1), 0) + prevPoint.vel);
      Vec3f dCurPoint  = curPoint.pos  + h*(dqdot.block<3, 1>(3*i,     0) + curPoint.vel);
      Vec3f dNextPoint = nextPoint.pos + h*(dqdot.block<3, 1>(3*(i+1), 0) + nextPoint.vel);
      Vec3f dPrevSeg   = dCurPoint - dPrevPoint;
      Vec3f dNextSeg   = dNextPoint - dCurPoint;
    
      // WARNING: assumes that twist in the material curvature changes minimally
      // between Newton iterations. This may not be the case.
    
      Vec3f tPrev = dPrevSeg.normalized();
      Vec3f tNext = dNextSeg.normalized();
      float chi = 1 + (tPrev.dot(tNext));
      Vec3f tTilde = (tPrev + tNext)/chi;
      Vec3f d1 = prevSeg.m1() + nextSeg.m1();
      Vec3f d1tilde = d1/chi;
      Vec3f d2 = prevSeg.m2() + nextSeg.m2();
      Vec3f d2tilde = d2/chi;
      Vec3f curveBinorm = (2*tPrev.cross(tNext))/chi;
      Vec2f matCurve = 0.5*Vec2f(d2.dot(curveBinorm), -d1.dot(curveBinorm));
      
      Vec3f gradK1ePrev = (-matCurve.x()*tTilde + tNext.cross(d2tilde)) / dPrevSeg.norm();
      Vec3f gradK1eNext = (-matCurve.x()*tTilde - tPrev.cross(d2tilde)) / dNextSeg.norm();
      Vec3f gradK2ePrev = (-matCurve.y()*tTilde + tNext.cross(d1tilde)) / dPrevSeg.norm();
      Vec3f gradK2eNext = (-matCurve.y()*tTilde - tPrev.cross(d1tilde)) / dNextSeg.norm();
      
      // WARNING: assumes that the bending matrix is the identity.
      
      Vec2f& dRestCurve = restCurve[i-1];
      // b11*2*(k1-restk1) + (b21+b12)(k2-restk2)
      float k1coeff = 2*(matCurve.x()-dRestCurve.x());
      // b22*2*(k2-restk2) + (b21+b12)(k1-restk1)
      float k2coeff = 2*(matCurve.y()-dRestCurve.y());
      float totalcoeff = 0.5/voronoiCell[i-1];
      
      Vec3f gradePrev = totalcoeff * (gradK1ePrev * k1coeff + gradK2ePrev * k2coeff);
      Vec3f gradeNext = totalcoeff * (gradK1eNext * k1coeff + gradK2eNext * k2coeff);
      
      typedef Eigen::Matrix3f Mat3f;
      
      Mat3f tTilde2 = tTilde*tTilde.transpose();
      
      Mat3f tNextxd2TildextTilde = (tNext.cross(d2tilde))*tTilde.transpose();
      Mat3f tPrevxd2TildextTilde = (tPrev.cross(d2tilde))*tTilde.transpose();
      Mat3f tNextxd1TildextTilde = (tNext.cross(d1tilde))*tTilde.transpose();
      Mat3f tPrevxd1TildextTilde = (tPrev.cross(d1tilde))*tTilde.transpose();
      
      Mat3f d2TildeCross, d1TildeCross;
      d2TildeCross << 0, -d2tilde.z(), d2tilde.y(),
                      d2tilde.z(), 0, -d2tilde.x(),
                      -d2tilde.y(), d2tilde.x(), 0;
      d1TildeCross << 0, -d1tilde.z(), d1tilde.y(),
                      d1tilde.z(), 0, -d1tilde.x(),
                      -d1tilde.y(), d1tilde.x(), 0;
      
      Mat3f hessK1ePrev2 = 2*matCurve.x()*tTilde2-tNextxd2TildextTilde-tNextxd2TildextTilde.transpose();
      hessK1ePrev2 -= (matCurve.x()/chi)*(Mat3f::Identity() - (tPrev*tPrev.transpose()));
      hessK1ePrev2 += 0.25*(curveBinorm*prevSeg.m2().transpose()+prevSeg.m2()*curveBinorm.transpose());
      hessK1ePrev2 /= dPrevSeg.dot(dPrevSeg);
      
      Mat3f hessK2ePrev2 = 2*matCurve.y()*tTilde2-tNextxd1TildextTilde-tNextxd1TildextTilde.transpose();
      hessK2ePrev2 -= (matCurve.y()/chi)*(Mat3f::Identity() - (tPrev*tPrev.transpose()));
      hessK2ePrev2 += 0.25*(curveBinorm*prevSeg.m1().transpose()+prevSeg.m1()*curveBinorm.transpose());
      hessK2ePrev2 /= dPrevSeg.dot(dPrevSeg);
      
      Mat3f hessK1eNext2 = 2*matCurve.x()*tTilde2+tPrevxd2TildextTilde+tPrevxd2TildextTilde.transpose();
      hessK1eNext2 -= (matCurve.x()/chi)*(Mat3f::Identity() - (tNext*tNext.transpose()));
      hessK1eNext2 += 0.25*(curveBinorm*nextSeg.m2().transpose()+nextSeg.m2()*curveBinorm.transpose());
      hessK1eNext2 /= dNextSeg.dot(dNextSeg);
      
      Mat3f hessK2eNext2 = 2*matCurve.y()*tTilde2+tPrevxd1TildextTilde+tPrevxd1TildextTilde.transpose();
      hessK2eNext2 -= (matCurve.y()/chi)*(Mat3f::Identity() - (tNext*tNext.transpose()));
      hessK2eNext2 += 0.25*(curveBinorm*nextSeg.m1().transpose()+nextSeg.m1()*curveBinorm.transpose());
      hessK2eNext2 /= dNextSeg.dot(dNextSeg);
      
      Mat3f hessK1ePreveNext = (-matCurve.x()/chi)*(Mat3f::Identity() + (tPrev * tNext.transpose()));
      hessK1ePreveNext += (2*matCurve.x()*tTilde2) - tNextxd2TildextTilde + tPrevxd2TildextTilde.transpose();
      hessK1ePreveNext -= d2TildeCross;
      hessK1ePreveNext /= (dNextSeg.norm() * dPrevSeg.norm());
      
      Mat3f hessK2ePreveNext = (-matCurve.y()/chi)*(Mat3f::Identity() + (tPrev * tNext.transpose()));
      hessK2ePreveNext += (2*matCurve.y()*tTilde2) - tNextxd1TildextTilde + tPrevxd1TildextTilde.transpose();
      hessK2ePreveNext -= d1TildeCross;
      hessK2ePreveNext /= (dNextSeg.norm() * dPrevSeg.norm());
      
      
      Mat3f hessePrev2 = 2*gradK1ePrev*gradK1ePrev.transpose() + k1coeff * hessK1ePrev2;
      hessePrev2 += 2*gradK2ePrev*gradK2ePrev.transpose() + k2coeff * hessK2ePrev2;
      hessePrev2 *= totalcoeff*h*h;
      
      Mat3f hesseNext2 = 2*gradK1eNext*gradK1eNext.transpose() + k1coeff * hessK1eNext2;
      hesseNext2 += 2*gradK2eNext*gradK2eNext.transpose() + k2coeff * hessK2eNext2;
      hesseNext2 *= totalcoeff*h*h;
      
      Mat3f  hessePreveNext = 2*gradK1ePrev*gradK1eNext.transpose() + k1coeff * hessK1ePreveNext;
      hessePreveNext += 2*gradK2ePrev*gradK2eNext.transpose() + k2coeff * hessK2ePreveNext;
      hessePreveNext *= totalcoeff*h*h;
      
      Fx.block<3,1>(3*(i-1), 0) += h*gradePrev;
      Fx.block<3,1>(3*i,     0) += h*(gradePrev + gradeNext);
      Fx.block<3,1>(3*(i+1), 0) += h*gradeNext;
      
      for(int j=0; j<3; j++) {
        for(int k=0; k<3; k++) {
          GradFx.push_back(Triplet(3*(i-1)+j, 3*(i-1)+k, hessePrev2(j, k)));
          GradFx.push_back(Triplet(3*i+j,     3*(i-1)+k, hessePreveNext(k, j) - hessePrev2(k, j)));
          GradFx.push_back(Triplet(3*(i-1)+j, 3*i+k,     hessePreveNext(j, k) - hessePrev2(j, k)));
          GradFx.push_back(Triplet(3*i+j,     3*i+k,     hessePrev2(j, k) + hesseNext2(j, k) - hessePreveNext(j, k) - hessePreveNext(k, j)));
          GradFx.push_back(Triplet(3*(i+1)+j, 3*i+k,     hessePreveNext(k, j) - hesseNext2(k, j)));
          GradFx.push_back(Triplet(3*i+j,     3*(i+1)+k, hessePreveNext(j, k) - hesseNext2(j, k)));
          GradFx.push_back(Triplet(3*(i+1)+j, 3*(i+1)+k, hesseNext2(j, k)));
          
          GradFx.push_back(Triplet(3*(i-1)+j, 3*(i+1)+k, -hessePreveNext(j, k)));
          GradFx.push_back(Triplet(3*(i+1)+j, 3*(i-1)+k, -hessePreveNext(k, j)));
        }
      }
      
/*
      Eigen::Matrix<float, 9, 9> hess;
      hess.block<3,3>(0, 0) = hessePrev2;
      hess.block<3,3>(3, 0) = hessePreveNext.transpose() - hessePrev2;
      hess.block<3,3>(0, 3) = hessePreveNext - hessePrev2;
      hess.block<3,3>(3, 3) = hessePrev2 + hesseNext2 - hessePreveNext - hessePreveNext.transpose();
      hess.block<3,3>(6, 3) = hessePreveNext.transpose() - hesseNext2;
      hess.block<3,3>(3, 6) = hessePreveNext - hesseNext2;
      hess.block<3,3>(6, 6) = hesseNext2;
      
      hess.block<3,3>(0, 6) = -hessePreveNext;
      hess.block<3,3>(6, 0) = -hessePreveNext.transpose();
      if (i == y.numCPs()-2) {
        std::cout << hess << "\n\n";
      }
 */
      
      
#endif //ifdef ENABLE_AUTODIFF
      
    }
    
    Profiler::stop("Bend Eval");
    
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
  float youngsModulus = 2e7;
  float xArea = constants::pi * constants::radius * constants::radius;
  float stretchScalar = xArea * youngsModulus;
  
public:
  Stretching(const Yarn& y, EvalType et) : YarnEnergy(y, et) { }
  
  void eval(VecXf& Fx, std::vector<Triplet>& GradFx, const VecXf& dqdot, Clock& c) {
    Profiler::start("Stretch Eval");
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
    
    Profiler::stop("Stretch Eval");
  }
#undef NUM_VARS
};


class Twisting : public YarnEnergy {
private:
  std::vector<float> voronoiCell;
  
  float shearModulus = 1e-5;
  float xArea = constants::pi * constants::radius * constants::radius;
  float twistMod = xArea * shearModulus * constants::radius * constants::radius / 2;

public:
  Twisting(const Yarn& y, EvalType et) : YarnEnergy(y, et) {
    for (int i=1; i<y.numCPs()-1; i++) {
      const Segment& ePrev = y.rest().segments[i-1];
      const Segment& eNext = y.rest().segments[i];
      voronoiCell.push_back(0.5*(ePrev.length()+eNext.length()));
    }
  }
  
  void eval(VecXf& Fx, std::vector<Triplet>& GradFx, const VecXf& dqdot, Clock& c) {
    for (int i=1; i<y.numCPs()-1; i++) {
      const Segment& segPrev = y.cur().segments[i-1];
      const Segment& segNext = y.cur().segments[i];

      Vec3f ePrev = segPrev.vec();
      Vec3f eNext = segNext.vec();
      Vec3f tPrev = ePrev.normalized();
      Vec3f tNext = eNext.normalized();
      
      Vec3f curveBinorm = (2*tPrev.cross(tNext))/(1+tPrev.dot(tNext));
      float dThetaHat = twistMod * (segNext.getRefTwist() - (segNext.getRot() - segPrev.getRot())) / voronoiCell[i-1];
      Vec3f dxi = curveBinorm / y.rest().segments[i].length() - curveBinorm / y.rest().segments[i-1].length();
      Fx.block<3,1>(3*i, 0) += c.timestep() * dThetaHat * dxi;
    }

  }
};

class IntContact : public YarnEnergy {
private:
  inline float f(float d) {
    return 1/d/d + d*d - 2;
  }
  
  inline float df(float d) {
    return -2/d/d/d + 2*d;
  }
  
  const float contactMod = 1e-3;

public:
  IntContact(const Yarn& y, EvalType et) : YarnEnergy(y, et) {}
  
  void eval(VecXf& Fx, std::vector<Triplet>& GradFx, const VecXf& dqdot, Clock& c) {
    for (int i=0; i<y.numSegs(); i++) {
      for (int j=0; j<y.numSegs(); j++) {
        if (i-j <= 1 && j-i <= 1) continue; // Don't evaluate adjacent edges.
        float h = c.timestep();
        const Segment& e1 = y.cur().segments[i];
        const Segment& e2 = y.cur().segments[j];
        
        Vec3f e1p1 = e1.getFirst().pos;  // + h*(dqdot.block<3,1>(3*i, 0) + e1.getFirst().vel);
        Vec3f e1p2 = e1.getSecond().pos; // + h*(dqdot.block<3,1>(3*i, 0) + e1.getSecond().vel);
        Vec3f e2p1 = e2.getFirst().pos;  // + h*(dqdot.block<3,1>(3*j, 0) + e2.getFirst().vel);
        Vec3f e2p2 = e2.getSecond().pos; // + h*(dqdot.block<3,1>(3*j, 0) + e2.getSecond().vel);
        
        Vec3f e1mid = (e1p1 + e1p2) / 2;
        Vec3f e2mid = (e2p1 + e2p2) / 2;
        
        if ((e1mid - e2mid).norm() > fmaxf(e1.length(), e2.length())/2) continue;
        
        // Splines are close, evaluate them.
        
        // FIXME: evaluate ends!!!
        if (i == 0 || j == 0 || i == y.numSegs()-1 || j == y.numSegs()-1) continue;
        
        /*
        if  (i == 1) {
          const Segment& e1next = y.cur().segments[i+1];
          s1(e1.getFirst(), e1.getSecond(), e1next.getSecond(), false);
        } else if (i == y.numSegs()-1) {
          const Segment& e1prev = y.cur().segments[i-1];
          s1(e1prev.getFirst(), e1.getFirst(), e1.getSecond(), true);
        } else {
          const Segment& e1prev = y.cur().segments[i-1];
          const Segment& e1next = y.cur().segments[i+1];
          s1(e1prev.getFirst(), e1.getFirst(), e1.getSecond(), e1next.getSecond());
        }
        */
        
        const Segment& e1prev = y.cur().segments[i-1];
        const Segment& e1next = y.cur().segments[i+1];
        Spline s1(e1prev.getFirst(), e1.getFirst(), e1.getSecond(), e1next.getSecond());

        const Segment& e2prev = y.cur().segments[j-1];
        const Segment& e2next = y.cur().segments[j+1];
        Spline s2(e2prev.getFirst(), e2.getFirst(), e2.getSecond(), e2next.getSecond());
        
        // FIXME: hack to get an estimate of spline length.
        float l1 = e1.length();
        float l2 = e2.length();
        
        for (int n=0; n<=12; n++) { // TODO: remove hardcoded 12
          for (int m=0; m<=12; m++) {
            float t1 = ((float) n) / 12;
            float t2 = ((float) m) / 12;
            Vec3f p1 = s1.eval(t1, false);
            Vec3f p2 = s2.eval(t2, false);
            float dist = (p1 - p2).norm();
            float coeff = contactMod*l1*l2*df(dist / 2 / constants::radius)/2/constants::radius;
            
            for (int k=0; k<4; k++) {
              Fx.block<3,1>(3*(i-1+k), 0) += h*coeff*s1.distGrad(k, t1, p2);
              Fx.block<3,1>(3*(j-1+k), 0) += h*coeff*s2.distGrad(k, t2, p1);
            }
          }
        }
        
      }
    }
  }
};


 
#endif
