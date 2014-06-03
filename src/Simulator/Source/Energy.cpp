//
//  Energy.cpp
//  Visualizer
//
//  Created by eschweickart on 4/14/14.
//
//

#include "Energy.h"

void pushBackIfNotZero(std::vector<Triplet>& GradFx, Triplet value) {
  if (value.value() != 0.0f) {
    GradFx.push_back(value);
  }
}

void const YarnEnergy::draw() {
  if (frames.empty()) return;
  for (std::function<void(void)> f : frames) {
    f();
  }
}

// GRAVITY

Gravity::Gravity(const Yarn& y, EvalType et, Vec3f dir) : YarnEnergy(y, et), dir(dir) {}

bool Gravity::eval(const VecXf& dqdot, Clock& c, VecXf& Fx, std::vector<Triplet>* GradFx) {
  assert(et == Explicit && "Unsupported EvalType");
  for (int i=0; i<y.numCPs(); i++) {
    Fx.block<3,1>(3*i, 0) -= dir * c.timestep();
  }
  return true;
}

// SPRING

Spring::Spring(const Yarn& y, EvalType et, size_t index, float stiffness) :
  YarnEnergy(y, et), index(index), stiffness(stiffness) {}

bool Spring::eval(const VecXf& dqdot, Clock& c, VecXf& Fx, std::vector<Triplet>* GradFx) {
  // Drawing ops
  frames.clear();
  ci::Vec3f ciclamp = toCi(clamp);
  ci::Vec3f ciindex = toCi(y.cur().points[index].pos);
  frames.push_back([ciclamp, ciindex] () {
    ci::gl::color(1.0f, 0.5f, 0.0f);
    ci::gl::drawLine(ciclamp, ciindex);
  });
  
  const float h = c.timestep();
  if (et == Explicit) {
    Fx.block<3,1>(3*index, 0) -= (clamp - y.cur().points[index].pos) * stiffness * h;
  } else {
    Fx.block<3,1>(3*index, 0) -= (clamp - (y.cur().points[index].pos +
                                           h*(y.cur().points[index].vel
                                              + dqdot.block<3,1>(3*index, 0)))) * stiffness * h;
    
    if (GradFx) {
      GradFx->push_back(Triplet(3*index,   3*index,   h*h*stiffness));
      GradFx->push_back(Triplet(3*index+1, 3*index+1, h*h*stiffness));
      GradFx->push_back(Triplet(3*index+2, 3*index+2, h*h*stiffness));
    }
  }
  
  return true;
}


void Spring::setClamp(Vec3f newClamp) { clamp = newClamp; }


// MOUSE SPRING

MouseSpring::MouseSpring(const Yarn& y, EvalType et, size_t index, float stiffness) :
  YarnEnergy(y, et), index(index), stiffness(stiffness) {}

bool MouseSpring::eval(const VecXf& dqdot, Clock& c, VecXf& Fx, std::vector<Triplet>* GradFx) {
  frames.clear();
  if (!mouseDown) return true;
  assert(mouseSet && "Set the mouse position each time you call eval()!");
  const float h = c.timestep();
  if (et == Explicit) {
    Fx.block<3,1>(3*index, 0) -= (mouse - y.cur().points[index].pos) * stiffness * h;
  } else {
    Fx.block<3,1>(3*index, 0) -= (mouse - (y.cur().points[index].pos +
                                           h*(y.cur().points[index].vel
                                              + dqdot.block<3, 1>(3*index, 0)))) * stiffness * h;
    if (GradFx) {
      GradFx->push_back(Triplet(3*index,   3*index,   h*h*stiffness));
      GradFx->push_back(Triplet(3*index+1, 3*index+1, h*h*stiffness));
      GradFx->push_back(Triplet(3*index+2, 3*index+2, h*h*stiffness));
    }
  }
  
  // Drawing Stuff
  ci::Vec3f cimouse = toCi(mouse);
  ci::Vec3f ciindex = toCi(y.cur().points[index].pos);
  frames.push_back([cimouse, ciindex] () {
    ci::gl::color(1.0f, 0.0f, 0.0f);
    ci::gl::drawLine(cimouse, ciindex);
  });
  
  return true;
}

void MouseSpring::setMouse(Vec3f newMouse, bool newDown) {
  mouse = newMouse;
  mouseDown = newDown;
  mouseSet = true;
}


// BENDING

Bending::Bending(const Yarn& y, EvalType et) : YarnEnergy(y, et) { }

// #define ENABLE_BEND_AUTODIFF
bool Bending::eval(const VecXf& dqdot, Clock& c, VecXf& Fx, std::vector<Triplet>* GradFx) {
  Profiler::start("Bend Eval");
  const float h = c.timestep();
  
  //VecXf test1 = VecXf::Zero(Fx.rows());
  //VecXf test2 = VecXf::Zero(Fx.rows());
  
  std::vector<Triplet> test1;
  std::vector<Triplet> test2;
  
  for (int i=1; i<y.numCPs()-1; i++) {
    const Segment&   prevSeg   = y.cur().segments[i-1];
    const Segment&   nextSeg   = y.cur().segments[i];
    Vec3f prevPoint = y.cur().points[i-1].pos;
    Vec3f curPoint  = y.cur().points[i].pos;
    Vec3f nextPoint = y.cur().points[i+1].pos;
    
    if (et == Implicit) {
      prevPoint += h*(dqdot.block<3,1>(3*(i-1), 0) + y.cur().points[i-1].vel);
      curPoint  += h*(dqdot.block<3,1>(3*i,     0) + y.cur().points[i].vel);
      nextPoint += h*(dqdot.block<3,1>(3*(i+1), 0) + y.cur().points[i+1].vel);
    }
    
#ifdef ENABLE_BEND_AUTODIFF
#define NUM_VARS 9
    typedef Eigen::Matrix<float, NUM_VARS, 1> Gradient;
    typedef Eigen::Matrix<float, NUM_VARS, NUM_VARS> Hessian;
    typedef DScalar2<float, NUM_VARS, Gradient, Hessian> DScalar;
    typedef DScalar::DVector3 DVector3;
    typedef DScalar::DVector2 DVector2;
    
    DVector3 dvPrevPoint(DScalar(0, prevPoint.x()),
                         DScalar(1, prevPoint.y()),
                         DScalar(2, prevPoint.z()));
    
    DVector3 dvCurPoint(DScalar(3, curPoint.x()),
                        DScalar(4, curPoint.y()),
                        DScalar(5, curPoint.z()));
    
    DVector3 dvNextPoint(DScalar(6, nextPoint.x()),
                         DScalar(7, nextPoint.y()),
                         DScalar(8, nextPoint.z()));
    
    
    DVector3 dvPrevSeg = dvCurPoint - dvPrevPoint;
    DVector3 dvNextSeg = dvNextPoint - dvCurPoint;
    assert(dvPrevSeg.norm() != 0.0f && dvNextSeg.norm() != 0.0f && "Edge length is 0");
    
    DVector3 dvPrevSegN = dvPrevSeg.normalized();
    DVector3 dvNextSegN = dvNextSeg.normalized();
    DScalar dotProd = dvPrevSegN.dot(dvNextSegN);
    assert(dotProd > -1.0f && "Segments are pointing in exactly opposite directions");
    
    DVector3 dvcurveBinorm = (DScalar(2.0f)*dvPrevSegN.cross(dvNextSegN))/(1.0f+dotProd);
    
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
    DVector2 dvmatCurve = DScalar(0.5f)*(matCurvePrev + matCurveNext);
    
    DVector2 restMatCurve(DScalar(y.restCurve(i).x()), DScalar(y.restCurve(i).y()));
    
    // TODO: bending matrix may not be I
    DVector2 curveDiff = dvmatCurve - restMatCurve;
    DScalar bendEnergy = (0.5f/y.restVoronoiLength(i))*curveDiff.dot(curveDiff);
    
    Gradient grad = bendEnergy.getGradient();
    Fx.block<NUM_VARS, 1>(3*(i-1), 0) += h*y.bendCoeff()*grad;
    //test1.block<9,1>(3*(i-1), 0) += grad;
    
    if (et == Implicit && GradFx) {
      Hessian hess = bendEnergy.getHessian();
      for (int j=0; j<NUM_VARS; j++) {
        for (int k=0; k<NUM_VARS; k++) {
          float val = h*h*y.bendCoeff()*hess(j,k);
          CHECK_NAN(val);
          pushBackIfNotZero(*GradFx, Triplet(3*(i-1)+j, 3*(i-1)+k, val));
          pushBackIfNotZero(test1, Triplet(3*(i-1)+j, 3*(i-1)+k, hess(j,k)));
        }
      }
    }
    
    
    /*
    if (i == y.numCPs()-2) {
      std::cout << hess << "\n\n";
    }
     */
    
#undef NUM_VARS
#else // ifdef ENABLE_BEND_AUTODIFF
    
    Vec3f prevVec   = curPoint - prevPoint;
    Vec3f nextVec   = nextPoint - curPoint;
    assert(prevVec.norm() > 0.0f && nextVec.norm() > 0.0f && "Segment length is zero!");
    
    // WARNING: assumes that twist in the material curvature changes minimally
    // between Newton iterations. This may not be the case.
    
    Vec3f tPrev = prevVec.normalized();
    Vec3f tNext = nextVec.normalized();
    float chi = 1.0f + (tPrev.dot(tNext));
    assert(chi > 0.0f && "Segments are pointing in exactly opposite directions!");
    Vec3f tTilde = (tPrev + tNext)/chi;
    Vec3f d1 = prevSeg.m1() + nextSeg.m1();
    Vec3f d1tilde = d1/chi;
    Vec3f d2 = prevSeg.m2() + nextSeg.m2();
    Vec3f d2tilde = d2/chi;
    Vec3f curveBinorm = (2.0f*tPrev.cross(tNext))/chi; // Verified
    Vec2f matCurve = 0.5f*Vec2f(d2.dot(curveBinorm), -d1.dot(curveBinorm)); // Verified
    
    Vec3f gradK1ePrev = (-matCurve.x()*tTilde + tNext.cross(d2tilde)) / prevVec.norm(); // Verified
    Vec3f gradK1eNext = (-matCurve.x()*tTilde - tPrev.cross(d2tilde)) / nextVec.norm(); // Verified
    Vec3f gradK2ePrev = (-matCurve.y()*tTilde - tNext.cross(d1tilde)) / prevVec.norm(); // Verified
    Vec3f gradK2eNext = (-matCurve.y()*tTilde + tPrev.cross(d1tilde)) / nextVec.norm(); // Verified
    
    
    // WARNING: assumes that the bending matrix is the identity.
    
    Vec2f restCurveVec = y.restCurve(i);
    // b11*2*(k1-restk1) + (b21+b12)(k2-restk2)
    float k1coeff = 2.0f * (matCurve.x()-restCurveVec.x()); // Verified
    // b22*2*(k2-restk2) + (b21+b12)(k1-restk1)
    float k2coeff = 2.0f * (matCurve.y()-restCurveVec.y()); // Verified
    float totalcoeff = 0.5f * y.bendCoeff()/y.restVoronoiLength(i);
    
    Vec3f gradePrev = totalcoeff * (gradK1ePrev * k1coeff + gradK2ePrev * k2coeff); // Verified
    Vec3f gradeNext = totalcoeff * (gradK1eNext * k1coeff + gradK2eNext * k2coeff); // Verified
    
    if (et == Implicit && GradFx) {
      typedef Eigen::Matrix3f Mat3f;
      
      Mat3f tTilde2 = tTilde*tTilde.transpose();
      
      Mat3f tNextxd2TildextTilde = (tNext.cross(d2tilde))*tTilde.transpose();
      Mat3f tPrevxd2TildextTilde = (tPrev.cross(d2tilde))*tTilde.transpose();
      Mat3f tNextxd1TildextTilde = (tNext.cross(d1tilde))*tTilde.transpose();
      Mat3f tPrevxd1TildextTilde = (tPrev.cross(d1tilde))*tTilde.transpose();
      
      Mat3f d2TildeCross, d1TildeCross;
      d2TildeCross << 0.0f, -d2tilde.z(), d2tilde.y(),
      d2tilde.z(), 0.0f, -d2tilde.x(),
      -d2tilde.y(), d2tilde.x(), 0.0f;
      d1TildeCross << 0.0f, -d1tilde.z(), d1tilde.y(),
      d1tilde.z(), 0.0f, -d1tilde.x(),
      -d1tilde.y(), d1tilde.x(), 0.0f;
      
      Mat3f hessK1ePrev2 = 2.0f*matCurve.x()*tTilde2-tNextxd2TildextTilde-tNextxd2TildextTilde.transpose();
      hessK1ePrev2 -= (matCurve.x()/chi)*(Mat3f::Identity() - (tPrev*tPrev.transpose()));
      hessK1ePrev2 += 0.25f*(curveBinorm*prevSeg.m2().transpose()+prevSeg.m2()*curveBinorm.transpose());
      hessK1ePrev2 /= prevVec.dot(prevVec);
      
      Mat3f hessK2ePrev2 = 2.0f*matCurve.y()*tTilde2+tNextxd1TildextTilde+tNextxd1TildextTilde.transpose();
      hessK2ePrev2 -= (matCurve.y()/chi)*(Mat3f::Identity() - (tPrev*tPrev.transpose()));
      hessK2ePrev2 += 0.25f*(curveBinorm*prevSeg.m1().transpose()+prevSeg.m1()*curveBinorm.transpose());
      hessK2ePrev2 /= prevVec.dot(prevVec);
      
      Mat3f hessK1eNext2 = 2.0f*matCurve.x()*tTilde2+tPrevxd2TildextTilde+tPrevxd2TildextTilde.transpose();
      hessK1eNext2 -= (matCurve.x()/chi)*(Mat3f::Identity() - (tNext*tNext.transpose()));
      hessK1eNext2 += 0.25f*(curveBinorm*nextSeg.m2().transpose()+nextSeg.m2()*curveBinorm.transpose());
      hessK1eNext2 /= nextVec.dot(nextVec);
      
      Mat3f hessK2eNext2 = 2.0f*matCurve.y()*tTilde2-tPrevxd1TildextTilde-tPrevxd1TildextTilde.transpose();
      hessK2eNext2 -= (matCurve.y()/chi)*(Mat3f::Identity() - (tNext*tNext.transpose()));
      hessK2eNext2 += 0.25f*(curveBinorm*nextSeg.m1().transpose()+nextSeg.m1()*curveBinorm.transpose());
      hessK2eNext2 /= nextVec.dot(nextVec);
      
      Mat3f hessK1ePreveNext = (-matCurve.x()/chi)*(Mat3f::Identity() + (tPrev * tNext.transpose()));
      hessK1ePreveNext += (2.0f*matCurve.x()*tTilde2) - tNextxd2TildextTilde + tPrevxd2TildextTilde.transpose();
      hessK1ePreveNext -= d2TildeCross;
      hessK1ePreveNext /= (nextVec.norm() * prevVec.norm());
      
      Mat3f hessK2ePreveNext = (-matCurve.y()/chi)*(Mat3f::Identity() + (tPrev * tNext.transpose()));
      hessK2ePreveNext += (2.0f*matCurve.y()*tTilde2) - tNextxd1TildextTilde + tPrevxd1TildextTilde.transpose();
      hessK2ePreveNext -= d1TildeCross;
      hessK2ePreveNext /= (nextVec.norm() * prevVec.norm());
      
      
      /*
       typedef DScalar2<float, 3, Vec3f, Mat3f> TDS;
       typedef TDS::DVector3 TVec3;
       typedef TDS::DVector2 TVec2;
       
       TVec3 teprev = TVec3(TDS(0, dPrevSeg.x()), TDS(1, dPrevSeg.y()), TDS(2, dPrevSeg.z()));
       //    TVec3 teprev = TVec3(TDS(dPrevSeg.x()), TDS(dPrevSeg.y()), TDS(dPrevSeg.z()));
       TVec3 tenext = TVec3(TDS(dNextSeg.x()), TDS(dNextSeg.y()), TDS(dNextSeg.z()));
       //    TVec3 tenext = TVec3(TDS(0, dNextSeg.x()), TDS(1, dNextSeg.y()), TDS(3, dNextSeg.z()));
       TVec3 ttprev = teprev.normalized();
       TVec3 ttnext = tenext.normalized();
       TDS tchi = ttprev.dot(ttnext) + 1.0f;
       TVec3 tkb = TDS(2.0f)*(ttprev.cross(ttnext)) / tchi;
       TDS tk1 = TDS(0.5)*TVec3(TDS(d2.x()), TDS(d2.y()), TDS(d2.z())).dot(tkb);
       TDS tk2 = TDS(-.5)*TVec3(TDS(d1.x()), TDS(d1.y()), TDS(d1.z())).dot(tkb);
       //    Eigen::Matrix<float, 3, 2> t;
       //    t << tk2.getGradient(), gradK2eNext;
       //    std::cout << t << "\n\n";
       Mat3f tdiff = tk1.getHessian() - hessK1ePrev2;
       if (tdiff.cwiseAbs().maxCoeff() >= 1.0f)
       std::cout << tdiff << "\n\n";
       // std::cout << tdiff.norm() << " / " << tk1.getHessian().cwiseAbs().maxCoeff() << "\n";
       */
      
      
      /// Populate the 3x3 block matrix:
      /// 0 1 2
      /// 3 4 5
      /// 6 7 8
      
      Mat3f block[9];
      
      block[0] = k1coeff*hessK1ePrev2 + k2coeff*hessK2ePrev2 + 2.0f*(gradK1ePrev*gradK1ePrev.transpose() + gradK2ePrev*gradK2ePrev.transpose());
      
      block[1] = k1coeff*(hessK1ePreveNext-hessK1ePrev2) + k2coeff*(hessK2ePreveNext-hessK2ePrev2)
      + 2.0f*(-gradK1ePrev*(gradK1ePrev-gradK1eNext).transpose() + -gradK2ePrev*(gradK2ePrev-gradK2eNext).transpose());
      
      block[2] = -k1coeff*hessK1ePreveNext - k2coeff*hessK2ePreveNext + 2.0f*(-gradK1ePrev*gradK1eNext.transpose() + -gradK2ePrev*gradK2eNext.transpose());
      
      block[3] = k1coeff*(hessK1ePreveNext.transpose() - hessK1ePrev2) + k2coeff*(hessK2ePreveNext.transpose() - hessK2ePrev2) + 2.0f*((gradK1ePrev-gradK1eNext)*-gradK1ePrev.transpose() + (gradK2ePrev-gradK2eNext)*-gradK2ePrev.transpose());
      
      block[4] = k1coeff*(hessK1ePrev2 + hessK1eNext2 - hessK1ePreveNext - hessK1ePreveNext.transpose()) +
      k2coeff*(hessK2ePrev2 + hessK2eNext2 - hessK2ePreveNext - hessK2ePreveNext.transpose()) +
      2.0f*((gradK1ePrev-gradK1eNext)*(gradK1ePrev-gradK1eNext).transpose() + (gradK2ePrev-gradK2eNext)*(gradK2ePrev-gradK2eNext).transpose());
      
      block[5] = k1coeff*(hessK1ePreveNext - hessK1eNext2) + k2coeff*(hessK2ePreveNext - hessK2eNext2)
      + 2.0f*((gradK1ePrev-gradK1eNext)*gradK1eNext.transpose() + (gradK2ePrev-gradK2eNext)*gradK2eNext.transpose());
      
      block[6] = -k1coeff*hessK1ePreveNext.transpose() - k2coeff*hessK2ePreveNext.transpose()
      + 2.0f*(gradK1eNext*-gradK1ePrev.transpose() + gradK2eNext*-gradK2ePrev.transpose());
      
      block[7] = k1coeff*(hessK1ePreveNext.transpose()-hessK1eNext2) + k2coeff*(hessK2ePreveNext.transpose()-hessK2eNext2) + 2.0f*(gradK1eNext*(gradK1ePrev-gradK1eNext).transpose() + gradK2eNext*(gradK2ePrev-gradK2eNext).transpose());
      
      block[8] = k1coeff*hessK1eNext2 + k2coeff*hessK2eNext2 + 2.0f*(gradK1eNext*gradK1eNext.transpose() + gradK2eNext*gradK2eNext.transpose());
      
      
      for(int j=0; j<3; j++) {
        for(int k=0; k<3; k++) {
          
          pushBackIfNotZero(*GradFx, Triplet(3*(i-1)+j, 3*(i-1)+k, h*h*totalcoeff*block[0](j, k)));
          pushBackIfNotZero(*GradFx, Triplet(3*(i-1)+j, 3*i+k,     h*h*totalcoeff*block[1](j, k)));
          pushBackIfNotZero(*GradFx, Triplet(3*(i-1)+j, 3*(i+1)+k, h*h*totalcoeff*block[2](j, k)));
          
          pushBackIfNotZero(*GradFx, Triplet(3*i+j,     3*(i-1)+k, h*h*totalcoeff*block[3](j, k)));
          pushBackIfNotZero(*GradFx, Triplet(3*i+j,     3*i+k,     h*h*totalcoeff*block[4](j, k)));
          pushBackIfNotZero(*GradFx, Triplet(3*i+j,     3*(i+1)+k, h*h*totalcoeff*block[5](j, k)));
          
          pushBackIfNotZero(*GradFx, Triplet(3*(i+1)+j, 3*(i-1)+k, h*h*totalcoeff*block[6](j, k)));
          pushBackIfNotZero(*GradFx, Triplet(3*(i+1)+j, 3*i+k,     h*h*totalcoeff*block[7](j, k)));
          pushBackIfNotZero(*GradFx, Triplet(3*(i+1)+j, 3*(i+1)+k, h*h*totalcoeff*block[8](j, k)));
          
          /*
           test2.push_back(Triplet(3*(i-1)+j, 3*(i-1)+k, totalcoeff*block[0](j, k)));
           test2.push_back(Triplet(3*(i-1)+j, 3*i+k,     totalcoeff*block[1](j, k)));
           test2.push_back(Triplet(3*(i-1)+j, 3*(i+1)+k, totalcoeff*block[2](j, k)));
           
           test2.push_back(Triplet(3*i+j,     3*(i-1)+k, totalcoeff*block[3](j, k)));
           test2.push_back(Triplet(3*i+j,     3*i+k,     totalcoeff*block[4](j, k)));
           test2.push_back(Triplet(3*i+j,     3*(i+1)+k, totalcoeff*block[5](j, k)));
           
           test2.push_back(Triplet(3*(i+1)+j, 3*(i-1)+k, totalcoeff*block[6](j, k)));
           test2.push_back(Triplet(3*(i+1)+j, 3*i+k,     totalcoeff*block[7](j, k)));
           test2.push_back(Triplet(3*(i+1)+j, 3*(i+1)+k, totalcoeff*block[8](j, k)));
           */
        }
      }
    }
    
    Fx.block<3,1>(3*(i-1), 0) += h*-gradePrev;
    Fx.block<3,1>(3*i,     0) += h*(gradePrev - gradeNext);
    Fx.block<3,1>(3*(i+1), 0) += h*gradeNext;
    
    /*
    test2.block<3,1>(3*(i-1), 0) += -gradePrev;
    test2.block<3,1>(3*i,     0) += (gradePrev - gradeNext);
    test2.block<3,1>(3*(i+1), 0) += gradeNext;
    */
    
    /*
    if (i == y.numCPs()-2) {
      Eigen::Matrix<float, 9, 9> hess2;
      
      hess2.block<3,3>(0, 0) = totalcoeff*block[0];
      hess2.block<3,3>(0, 3) = totalcoeff*block[1];
      hess2.block<3,3>(0, 6) = totalcoeff*block[2];
      
      hess2.block<3,3>(3, 0) = totalcoeff*block[3];
      hess2.block<3,3>(3, 3) = totalcoeff*block[4];
      hess2.block<3,3>(3, 6) = totalcoeff*block[5];
      
      hess2.block<3,3>(6, 0) = totalcoeff*block[6];
      hess2.block<3,3>(6, 3) = totalcoeff*block[7];
      hess2.block<3,3>(6, 6) = totalcoeff*block[8];
      
      std::cout << (hess - hess2).cwiseAbs().maxCoeff() << " / " << hess.cwiseAbs().maxCoeff() << "\n";
    }
    */
    
    
#endif //ifdef ENABLE_BEND_AUTODIFF
  }
  
  /*
  Eigen::SparseMatrix<float> test1m(Fx.rows(), Fx.rows());
  Eigen::SparseMatrix<float> test2m(Fx.rows(), Fx.rows());
  test1m.setFromTriplets(test1.begin(), test1.end());
  test2m.setFromTriplets(test2.begin(), test2.end());
  
  
   std::cout << (test1m-test2m).norm() << "\n";
   */
  
  //std::cout << test1.block<3, 1>(3*(y.numCPs()-2), 0) << "\n";
  //std::cout << test2.block<3, 1>(3*(y.numCPs()-2), 0) << "\n\n";
  
  Profiler::stop("Bend Eval");
  
  return true;
}


// STRETCHING

Stretching::Stretching(const Yarn& y, EvalType et) : YarnEnergy(y, et) { }

void Stretching::suggestTimestep(Clock& c) {
  YarnEnergy::suggestTimestep(c);
  /*
  float maxStrain = 1;
  float minStrain = 1;
  for (int i=0; i<y.numSegs(); i++) {
    float ratio = y.cur().segments[i].length() / y.rest().segments[i].length();
    if (ratio > maxStrain) maxStrain = ratio;
    if (ratio < minStrain) minStrain = ratio;
  }
  if (maxStrain >= 2.0f || minStrain <= 0.5f) {
    c.suggestTimestep(c.minTimestep);
  } else if (maxStrain <= 1.3f && minStrain >= 0.77f) {
    c.suggestTimestep(c.maxTimestep);
  } else { // lerp between min and max timesteps
    float alpha = (maxStrain - 1.3f) / (2.0f - 1.3f);
    float beta  = 1.0f - ((minStrain - 0.5f) / (0.77f - 0.5f));
    float gamma = fmax(alpha, beta);
    c.suggestTimestep(c.minTimestep*gamma+c.maxTimestep*(1.0f-gamma));
  }
   */
}

//#define ENABLE_STRETCH_AUTODIFF
bool Stretching::eval(const VecXf& dqdot, Clock& c, VecXf& Fx, std::vector<Triplet>* GradFx) {
  Profiler::start("Stretch Eval");
  float h = c.timestep();
  
  for (int i=0; i<y.numSegs(); i++) {
    const Segment& seg = y.cur().segments[i];
    float restSegLength = y.rest().segments[i].length();
    Vec3f prevPoint = seg.getFirst().pos;
    Vec3f nextPoint = seg.getSecond().pos;
    
    if (et == Implicit) {
      prevPoint += h*(dqdot.block<3,1>(3*i,     0) + seg.getFirst().vel);
      nextPoint += h*(dqdot.block<3,1>(3*(i+1), 0) + seg.getSecond().vel);
    }
    
#ifdef ENABLE_STRETCH_AUTODIFF
#define NUM_VARS 6
    typedef Eigen::Matrix<float, NUM_VARS, 1> Gradient;
    typedef Eigen::Matrix<float, NUM_VARS, NUM_VARS> Hessian;
    typedef DScalar2<float, NUM_VARS, Gradient, Hessian> DScalar;
    typedef DScalar::DVector3 DVector3;
    
    DVector3 dPrevPoint(DScalar(0, prevPoint.x()),
                        DScalar(1, prevPoint.y()),
                        DScalar(2, prevPoint.z()));
    
    DVector3 dNextPoint(DScalar(3, nextPoint.x()),
                        DScalar(4, nextPoint.y()),
                        DScalar(5, nextPoint.z()));
    
    DScalar axialStrain = (dNextPoint - dPrevPoint).norm()/restSegLength - 1.0f;
    DScalar stretchEnergy = (0.5f * restSegLength) * axialStrain * axialStrain;
    
    Gradient grad = stretchEnergy.getGradient();
    
    Fx.block<NUM_VARS, 1>(3*i, 0) += h * stretchScalar * grad;
    
    if (et == Implicit && GradFx) {
      Hessian hess = stretchEnergy.getHessian();
    
      for (int j=0; j<NUM_VARS; j++) {
        for (int k=0; k<NUM_VARS; k++) {
          float val = h*h*stretchScalar*hess(j,k);
          pushBackIfNotZero(*GradFx, Triplet(3*i+j, 3*i+k, val));
        }
      }
    }
#undef NUM_VARS
#else // ifdef ENABLE_STRETCH_AUTODIFF

    Vec3f ej = nextPoint - prevPoint;
    float l = ej.norm();
    Vec3f myGrad = ((1.0f/restSegLength) - (1.0f/l)) * ej; // Verified
     
    /*
    typedef DScalar2<float, 3, Vec3f, Mat3f> TDS;
    typedef TDS::DVector3 TVec3;
    typedef TDS::DVector2 TVec2;
    TVec3 te(TDS(0, ej.x()), TDS(1, ej.y()), TDS(2, ej.z()));
    TDS tas = te.norm() / restSegLength - 1.0f;
    Vec3f diff = tas.getGradient() - ej/ej.norm()/restSegLength;
    if (diff.norm() > 1e-5f) {
      std::cout << "axial strain badness: " << diff << "\n";
    }
    */
    
    /*
    Gradient graddiff = grad;
    graddiff.block<3,1>(0, 0) += myGrad;
    graddiff.block<3,1>(3, 0) -= myGrad;
    if (graddiff.norm() > 1e-7f) std::cout << "stretch error: " << graddiff << "\n\n";
     */
    
    /*
    Hessian diff = hess;
    diff.block<3,3>(0,0) -= myHess;
    diff.block<3,3>(3,0) += myHess;
    diff.block<3,3>(0,3) += myHess;
    diff.block<3,3>(3,3) -= myHess;
    if (diff.norm() > .005f) {
      std::cout << "stretch error:\n" << diff << "\n\n";
    } else {
      std::cout << "success!\n";
    }
     */
    
    
    Fx.block<3,1>(3*i, 0) -= h * y.stretchCoeff() * myGrad;
    Fx.block<3,1>(3*(i+1), 0) += h * y.stretchCoeff() * myGrad;
    if (et == Implicit && GradFx) {
      typedef Eigen::Matrix3f Mat3f;
      Mat3f myHess = (Mat3f::Identity()/restSegLength -
                      (Mat3f::Identity()/l - ej * ej.transpose() / (l * l * l))); // Verified
      
      for (int j=0; j<3; j++) {
        for (int k=0; k<3; k++) {
          pushBackIfNotZero(*GradFx, Triplet(3*i+j,     3*i+k,     h*h*y.stretchCoeff()*myHess(j,k)));
          pushBackIfNotZero(*GradFx, Triplet(3*(i+1)+j, 3*i+k,     -h*h*y.stretchCoeff()*myHess(j,k)));
          pushBackIfNotZero(*GradFx, Triplet(3*i+j,     3*(i+1)+k, -h*h*y.stretchCoeff()*myHess(j,k)));
          pushBackIfNotZero(*GradFx, Triplet(3*(i+1)+j, 3*(i+1)+k, h*h*y.stretchCoeff()*myHess(j,k)));
        }
      }
    }
    
    
#endif // ifdef ENABLE_STRETCH_AUTODIFF
  }
  
  Profiler::stop("Stretch Eval");
  
  return true;
}



// TWISTING

Twisting::Twisting(const Yarn& y, EvalType et) : YarnEnergy(y, et) { }

#define DRAW_TWIST
bool Twisting::eval(const VecXf& dqdot, Clock& c, VecXf& Fx, std::vector<Triplet>* GradFx) {
#ifdef DRAW_TWIST
  frames.clear();
#endif // ifdef DRAW_TWIST
  assert(et == Explicit && "EvalType unsupported");
  for (int i=1; i<y.numCPs()-1; i++) {
    const Segment& segPrev = y.cur().segments[i-1];
    const Segment& segNext = y.cur().segments[i];
    
    Vec3f tPrev = segPrev.vec().normalized();
    Vec3f tNext = segNext.vec().normalized();
    float chi = 1.0f+tPrev.dot(tNext);
    assert(chi > 0.0f && "Segments are pointing in exactly opposite directions!");
    Vec3f curveBinorm = (2.0f*tPrev.cross(tNext))/chi;
    float dThetaHat = y.twistCoeff() * (segNext.getRefTwist() - (segNext.getRot() - segPrev.getRot())) / y.restVoronoiLength(i);
    Vec3f dxi = curveBinorm / y.rest().segments[i].length() - curveBinorm / y.rest().segments[i-1].length();
    Fx.block<3,1>(3*i, 0) += c.timestep() * dThetaHat * dxi;
    
#ifdef DRAW_TWIST
    ci::Vec3f delta = toCi(c.timestep() * dThetaHat * dxi);
    const Yarn* yp = &y;
    if (delta.length() > 0.0f) frames.push_back([i, delta, yp] () {
      ci::Vec3f s = toCi(yp->cur().points[i].pos);
      ci::Vec3f e = s - delta;
      ci::gl::color(ci::Color(1.0f, 0.5f, 0.5f));
      ci::gl::drawVector(s, e);
    });
#endif // ifdef DRAW_TWIST
  }
  return true;
}

// INTERNAL CONTACT

IntContact::IntContact(const Yarn& y, EvalType et) : YarnEnergy(y, et) {}

// FIXME: This is hacky and slow
void IntContact::suggestTimestep(Clock& c) {
  /*
  float minDist = 1e6f;
  const float r = constants::radius;
  for (int i=0; i<y.numSegs(); i++) {
    for (int j=i+2; j<y.numSegs(); j++) {
      if (i == 0 || j == 0 || i == y.numSegs()-1 || j == y.numSegs()-1) continue; // FIXME: evaluate ends
      if (i-j <= 3 && j-i <= 3) continue; // Don't evaluate close edges. FIXME: should be 1 ideally
      
      const Segment& e1 = y.cur().segments[i];
      const Segment& e2 = y.cur().segments[j];
      
      Vec3f e1p1 = e1.getFirst().pos;
      Vec3f e1p2 = e1.getSecond().pos;
      Vec3f e2p1 = e2.getFirst().pos;
      Vec3f e2p2 = e2.getSecond().pos;
      
      Vec3f e1mid = (e1p1 + e1p2) / 2.0f;
      Vec3f e2mid = (e2p1 + e2p2) / 2.0f;
      
      if ((e1mid - e2mid).norm() > fmaxf(e1.length(), e2.length())) continue;
      
      const Segment& e1prev = y.cur().segments[i-1];
      const Segment& e1next = y.cur().segments[i+1];
      Spline s1(e1prev.getFirst(), e1.getFirst(), e1.getSecond(), e1next.getSecond());
      
      const Segment& e2prev = y.cur().segments[j-1];
      const Segment& e2next = y.cur().segments[j+1];
      Spline s2(e2prev.getFirst(), e2.getFirst(), e2.getSecond(), e2next.getSecond());
     
      for (int n=0; n<=nb; n++) {
        for (int m=0; m<=nb; m++) {
          float t1 = ((float) n) / nb;
          float t2 = ((float) m) / nb;
          Vec3f p1 = s1.eval(t1);
          Vec3f p2 = s2.eval(t2);
          Vec3f v = p2 - p1;
          float norm = v.norm();
          if (norm < minDist) minDist = norm;
        }
      }
    }
  }
  c.suggestTimestep(minDist / (4.0f * r) * constants::INITIAL_TIMESTEP);
   */
  YarnEnergy::suggestTimestep(c);
}


#define DRAW_INT_CONTACT
bool IntContact::eval(const VecXf& dqdot, Clock& c, VecXf& Fx, std::vector<Triplet>* GradFx) {
#ifdef DRAW_INT_CONTACT
  frames.clear();
#endif // ifdef DRAW_INT_CONTACT
  const float h = c.timestep();
  const float r = y.radius();
  
  for (int i=0; i<y.numSegs(); i++) {
    for (int j=i+2; j<y.numSegs(); j++) {
      if (i == 0 || j == 0 || i == y.numSegs()-1 || j == y.numSegs()-1) continue; // FIXME: evaluate ends
      if (i-j <= 3 && j-i <= 3) continue; // Don't evaluate close edges. FIXME: should be 1 ideally

      const Segment& e1 = y.cur().segments[i];
      const Segment& e2 = y.cur().segments[j];
      
      Vec3f e1p1 = e1.getFirst().pos;
      Vec3f e1p2 = e1.getSecond().pos;
      Vec3f e2p1 = e2.getFirst().pos;
      Vec3f e2p2 = e2.getSecond().pos;
      
      if (et == Implicit) {
        e1p1 += h*(dqdot.block<3,1>(3*i,     0) + e1.getFirst().vel);
        e1p2 += h*(dqdot.block<3,1>(3*(i+1), 0) + e1.getSecond().vel);
        e2p1 += h*(dqdot.block<3,1>(3*j,     0) + e2.getFirst().vel);
        e2p2 += h*(dqdot.block<3,1>(3*(j+1), 0) + e2.getSecond().vel);
      }
      
      Vec3f e1mid = (e1p1 + e1p2) / 2.0f;
      Vec3f e2mid = (e2p1 + e2p2) / 2.0f;
      
      if ((e1mid - e2mid).norm() > fmaxf(e1.length(), e2.length())) {
        if (ptd) {
          std::pair<int, int> id = ptd->id();
          if (id.first == i && id.second == j) {
            delete ptd;
            ptd = nullptr;
          }
        }
        continue;
      }
      
      // Splines are close, evaluate them.
      
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
      float coeff = contactMod * l1 * l2;
      
      Vec3f ref = e1mid - e2mid;
      float s1dot = e1.getU().dot(ref);
      float s2dot = e2.getU().dot(ref);
      
      if (ptd) {
        std::pair<int, int> id = ptd->id();
        if (id.first == i && id.second == j) {
          if (ptd->pass(s1dot, s2dot)) {
            std::cout << "Pullthrough detected: " << i << " " << j << "\n";
          }
        }
      } else {
        ptd = new PTDetector(i, j, s1dot, s2dot);
      }
      
      typedef Eigen::Matrix3f Mat3f;
      Mat3f hess[8][8];
      if (et == Implicit && GradFx) {
        for (int k=0; k<8; k++) {
          for (int l=0; l<8; l++) {
            hess[k][l] = Mat3f::Zero();
          }
        }
      }
      
#ifdef DRAW_INT_CONTACT
      Vec3f s1draw[4];
      Vec3f s2draw[4];
      for (int k=0; k<4; k++) {
        s1draw[k] = Vec3f::Zero();
        s2draw[k] = Vec3f::Zero();
      }
#endif // ifdef DRAW_INT_CONTACT
      
      for (int n=0; n<=nb; n++) {
        for (int m=0; m<=nb; m++) {
          float t1 = ((float) n) / nb;
          float t2 = ((float) m) / nb;
          Vec3f p1 = s1.eval(t1);
          Vec3f p2 = s2.eval(t2);
          Vec3f v = p2 - p1; // FIXME: shouldn't this be p1 - p2??
          float norm = v.norm();
          if (norm > 2*r) continue; // Quadratures are not touching, will eval to 0
          Vec4f u1(t1*t1*t1, t1*t1, t1, 1.0f);
          Vec4f u2(t2*t2*t2, t2*t2, t2, 1.0f);
          
          float dfnorm = df(norm / 2.0f / r) / 2.0f / r;
          Vec3f gradBase = v/norm * dfnorm;
          
          float d2fnorm = d2f(norm / 2.0f / r) / 4.0f / r / r;
          Mat3f grad2Base = v*v.transpose()/norm/norm;
          Mat3f hessBase = (Mat3f::Identity() - grad2Base) * dfnorm / norm + grad2Base * d2fnorm;
          
          /// TESTING
          /*
          typedef Eigen::Matrix<float, 6, 1> Grad6;
          typedef Eigen::Matrix<float, 6, 6> Hess6;
          typedef DScalar2<float, 6, Grad6, Hess6> TDS;
          typedef TDS::DVector3 TVec3;
          
          Vec3f x10 = e1prev.getFirst().pos;
          Vec3f x11 = e1.getFirst().pos;
          Vec3f x12 = e1.getSecond().pos;
          Vec3f x13 = e1next.getSecond().pos;
          Vec3f x20 = e2prev.getFirst().pos;
          Vec3f x21 = e2.getFirst().pos;
          Vec3f x22 = e2.getSecond().pos;
          Vec3f x23 = e2next.getSecond().pos;
          
          TVec3 tx1(TDS(0, x10.x()), TDS(1, x10.y()), TDS(2, x10.z()));
          TVec3 tx2(TDS(3, x20.x()), TDS(4, x20.y()), TDS(5, x20.z()));
          
          Vec4f basis1(constants::basis[0].dot(u1), constants::basis[1].dot(u1),
                       constants::basis[2].dot(u1), constants::basis[3].dot(u1));
          Vec4f basis2(constants::basis[0].dot(u2), constants::basis[1].dot(u2),
                       constants::basis[2].dot(u2), constants::basis[3].dot(u2));
          
          TVec3 tp1(basis1[0]*tx1.x() + (basis1[1]*x11.x() + basis1[2]*x12.x() + basis1[3]*x13.x()),
                    basis1[0]*tx1.y() + (basis1[1]*x11.y() + basis1[2]*x12.y() + basis1[3]*x13.y()),
                    basis1[0]*tx1.z() + (basis1[1]*x11.z() + basis1[2]*x12.z() + basis1[3]*x13.z()));
          TVec3 tp2(basis2[0]*tx2.x() + (basis2[1]*x21.x() + basis2[2]*x22.x() + basis2[3]*x23.x()),
                    basis2[0]*tx2.y() + (basis2[1]*x21.y() + basis2[2]*x22.y() + basis2[3]*x23.y()),
                    basis2[0]*tx2.z() + (basis2[1]*x21.z() + basis2[2]*x22.z() + basis2[3]*x23.z()));
          
          assert(is_approx(tp1.x().getValue(), p1.x()) && is_approx(tp1.y().getValue(), p1.y()) && is_approx(tp1.z().getValue(), p1.z()));
          assert(is_approx(tp2.x().getValue(), p2.x()) && is_approx(tp2.y().getValue(), p2.y()) && is_approx(tp2.z().getValue(), p2.z()));
          
          TDS tnorm = (tp2 - tp1).norm();
          Grad6 tgrad = tnorm.getGradient();
          Hess6 thess = tnorm.getHessian();
          
          // tgrad.block<3,1>(0, 0) += u1.dot(constants::basis[0]) * gradBase;
          // tgrad.block<3,1>(3, 0) -= u2.dot(constants::basis[0]) * gradBase;
          
          Mat3f tgrad2 = (Mat3f::Identity() - grad2Base) / norm;
          thess.block<3,3>(0, 0) -= u1.dot(constants::basis[0]) * u1.dot(constants::basis[0]) * tgrad2;
          thess.block<3,3>(3, 0) += u1.dot(constants::basis[0]) * u2.dot(constants::basis[0]) * tgrad2;
          thess.block<3,3>(0, 3) += u1.dot(constants::basis[0]) * u2.dot(constants::basis[0]) * tgrad2;
          thess.block<3,3>(3, 3) -= u2.dot(constants::basis[0]) * u2.dot(constants::basis[0]) * tgrad2;
        
          if (thess.norm() > 1e-5f) {
            std::cout << "contact badness:\n" << thess << "\n\n";
          } else {
            std::cout << "contact success\n";
          }
          */
          ///
          
          for (int k=0; k<4; k++) {
            float c = u1.dot(constants::basis[k]);
            float d = u2.dot(constants::basis[k]);
            
            Fx.block<3,1>(3*(i-1+k), 0) -= h * coeff * c * gradBase; // Verified
            Fx.block<3,1>(3*(j-1+k), 0) += h * coeff * d * gradBase; // Verified
#ifdef DRAW_INT_CONTACT
            s1draw[k] += h * coeff * c * gradBase;
            s2draw[k] -= h * coeff * d * gradBase;
#endif // ifdef DRAW_INT_CONTACT
          }
          
          if (et == Implicit && GradFx) {
            for (int k=0; k<4; k++) {
              for (int l=0; l<4; l++) {
                float ck = u1.dot(constants::basis[k]);
                float cl = u1.dot(constants::basis[l]);
                float dk = u2.dot(constants::basis[k]);
                float dl = u2.dot(constants::basis[l]);
                
                hess[k][l]     += ck * cl * hessBase;
                hess[k+4][l]   -= ck * dl * hessBase;
                hess[k][l+4]   -= dk * cl * hessBase;
                hess[k+4][l+4] += dk * dl * hessBase;
              }
            }
          }
        }
      }
      
      if (et == Implicit && GradFx) {
        for (int k=0; k<8; k++) {
          for (int l=0; l<8; l++) {
            hess[k][l] *= h * h * coeff;
          }
        }

        for (int k=0; k<4; k++) {
          for (int l=0; l<4; l++) {
            for (int p=0; p<3; p++) {
              for (int q=0; q<3; q++) {
                pushBackIfNotZero(*GradFx, Triplet(3*(i-1+k)+p, 3*(i-1+l)+q, hess[k][l](p, q)));
                pushBackIfNotZero(*GradFx, Triplet(3*(i-1+k)+p, 3*(j-1+l)+q, hess[k+4][l](p, q)));
                pushBackIfNotZero(*GradFx, Triplet(3*(j-1+k)+p, 3*(i-1+l)+q, hess[k][l+4](p, q)));
                pushBackIfNotZero(*GradFx, Triplet(3*(j-1+k)+p, 3*(j-1+l)+q, hess[k+4][l+4](p, q)));
              }
            }
          }
        }
      }

#ifdef DRAW_INT_CONTACT
      for (int k=0; k<4; k++) {
        ci::Vec3f s1start = toCi(y.cur().points[i-1+k].pos);
        ci::Vec3f s1end = toCi(y.cur().points[i-1+k].pos + s1draw[k]);
        ci::Vec3f s2start = toCi(y.cur().points[j-1+k].pos);
        ci::Vec3f s2end = toCi(y.cur().points[j-1+k].pos + s2draw[k]);
        frames.push_back([s1start, s1end, s2start, s2end] () {
          ci::gl::color(ci::Color(0.0f, 1.0f, 0.0f));
          ci::gl::drawVector(s1start, s1end);
          ci::gl::drawVector(s2start, s2end);
        });
      }
#endif // ifdef DRAW_INT_CONTACT
      
    }
  }
  return true;
}
