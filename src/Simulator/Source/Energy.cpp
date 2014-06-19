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

void const YarnEnergy::draw(float scale) {
  if (frames.empty()) return;
  for (std::function<void(float)> f : frames) {
    f(scale);
  }
}

// GRAVITY

Gravity::Gravity(const Yarn& y, EvalType et, Vec3f dir) : YarnEnergy(y, et), dir(dir) { }

bool Gravity::eval(VecXf* Fx, std::vector<Triplet>* GradFx, const VecXf* offset) {
  if (Fx) {
    for (int i=0; i<y.numCPs(); i++) {
      Fx->block<3,1>(3*i, 0) += dir;
    }
  }
  return true;
}

// SPRING

Spring::Spring(const Yarn& y, EvalType et, size_t index, float stiffness) :
  YarnEnergy(y, et), index(index), stiffness(stiffness) {}

bool Spring::eval(VecXf* Fx, std::vector<Triplet>* GradFx, const VecXf* offset) {
  // Drawing ops
#ifdef DRAW_SPRING
  frames.clear();
  ci::Vec3f ciclamp = toCi(clamp);
  ci::Vec3f ciindex = toCi(y.cur().points[index].pos);
  frames.push_back([ciclamp, ciindex] (float scale) {
    ci::gl::color(1.0f, 0.5f, 0.0f);
    ci::gl::drawLine(ciclamp, ciindex);
  });
#endif // ifdef DRAW_SPRING
  
  Vec3f yPoint = y.cur().points[index].pos;
  if (offset) {
    yPoint += offset->block<3,1>(3*index, 0);
  }
  if (Fx) {
    Fx->block<3,1>(3*index, 0) += stiffness * (clamp - yPoint);
  }
  if (GradFx) {
    GradFx->push_back(Triplet(3*index,   3*index,   -stiffness));
    GradFx->push_back(Triplet(3*index+1, 3*index+1, -stiffness));
    GradFx->push_back(Triplet(3*index+2, 3*index+2, -stiffness));
  }
  
  return true;
}

void Spring::setClamp(Vec3f newClamp) { clamp = newClamp; }


// MOUSE SPRING

MouseSpring::MouseSpring(const Yarn& y, EvalType et, size_t index, float stiffness) :
  YarnEnergy(y, et), index(index), stiffness(stiffness) {}

bool MouseSpring::eval(VecXf* Fx, std::vector<Triplet>* GradFx, const VecXf* offset) {
  frames.clear();
  if (!mouseDown) return true;
  assert(mouseSet && "Set the mouse position each time you call eval()!");

  Vec3f yPoint = y.cur().points[index].pos;
  if (offset) {
    yPoint += offset->block<3,1>(3*index, 0);
  }
  if (Fx) {
    Fx->block<3,1>(3*index, 0) += stiffness * (mouse - yPoint);
  }
  if (GradFx) {
    GradFx->push_back(Triplet(3*index,   3*index,   -stiffness));
    GradFx->push_back(Triplet(3*index+1, 3*index+1, -stiffness));
    GradFx->push_back(Triplet(3*index+2, 3*index+2, -stiffness));
  }
#ifdef DRAW_MOUSE_SPRING
  // Drawing Stuff
  ci::Vec3f cimouse = toCi(mouse);
  ci::Vec3f ciindex = toCi(yPoint);
  frames.push_back([cimouse, ciindex] (float scale) {
    ci::gl::color(1.0f, 0.0f, 0.0f);
    ci::gl::drawVector(ciindex, cimouse);
  });
#endif // ifdef DRAW_MOUSE_SPRING
  
  return true;
}

void MouseSpring::setMouse(Vec3f newMouse, bool newDown) {
  mouse = newMouse;
  mouseDown = newDown;
  mouseSet = true;
}


// BENDING

Bending::Bending(const Yarn& y, EvalType et) : YarnEnergy(y, et) { }

bool Bending::eval(VecXf* Fx, std::vector<Triplet>* GradFx, const VecXf* offset) {
  Profiler::start("Bend Eval");
  
#ifdef DRAW_BENDING
  frames.clear();
  VecXf forces = VecXf::Zero(3*y.numCPs());
#endif //ifdef DRAW_BENDING
  
  for (int i=1; i<y.numCPs()-1; i++) {
    const Segment&   prevSeg   = y.cur().segments[i-1];
    const Segment&   nextSeg   = y.cur().segments[i];
    Vec3f prevPoint = y.cur().points[i-1].pos;
    Vec3f curPoint  = y.cur().points[i].pos;
    Vec3f nextPoint = y.cur().points[i+1].pos;
    
    if (offset) {
      prevPoint += offset->block<3,1>(3*(i-1), 0);
      curPoint  += offset->block<3,1>(3*i,     0);
      nextPoint += offset->block<3,1>(3*(i+1), 0);
    }
    
#ifdef ENABLE_BEND_AUTODIFF
    typedef Eigen::Matrix<float, 9, 1> Gradient;
    typedef Eigen::Matrix<float, 9, 9> Hessian;
    typedef DScalar2<float, 9, Gradient, Hessian> DScalar;
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
    DScalar bendEnergy = (-0.5f/y.restVoronoiLength(i))*curveDiff.dot(curveDiff);
    
    if (Fx) {
      Gradient grad = bendEnergy.getGradient();
      Fx->block<9,1>(3*(i-1), 0) += y.bendCoeff()*grad;
    }
    
    if (GradFx) {
      Hessian hess = bendEnergy.getHessian();
      for (int j=0; j<9; j++) {
        for (int k=0; k<9; k++) {
          float val = y.bendCoeff()*hess(j,k);
          CHECK_NAN(val);
          pushBackIfNotZero(*GradFx, Triplet(3*(i-1)+j, 3*(i-1)+k, val));
        }
      }
    }
    
#else // ifdef ENABLE_BEND_AUTODIFF
    
    Vec3f prevVec = curPoint - prevPoint;
    Vec3f nextVec = nextPoint - curPoint;
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
    
    if (Fx) {
      Vec3f gradePrev = totalcoeff * (gradK1ePrev * k1coeff + gradK2ePrev * k2coeff); // Verified
      Vec3f gradeNext = totalcoeff * (gradK1eNext * k1coeff + gradK2eNext * k2coeff); // Verified
    
      Fx->block<3,1>(3*(i-1), 0) += gradePrev;
      Fx->block<3,1>(3*i,     0) += (gradeNext - gradePrev);
      Fx->block<3,1>(3*(i+1), 0) += -gradeNext;
      
#ifdef DRAW_BENDING
      forces.block<3,1>(3*(i-1), 0) += gradePrev;
      forces.block<3,1>(3*i,     0) += (gradeNext - gradePrev);
      forces.block<3,1>(3*(i+1), 0) += -gradeNext;
#endif //ifdef DRAW_BENDING
    }
    
    if (GradFx) {
      typedef Eigen::Matrix<float, 9, 1> Gradient;
      typedef Eigen::Matrix<float, 9, 9> Hessian;
      typedef DScalar2<float, 9, Gradient, Hessian> DScalar;
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
      DScalar bendEnergy = (-0.5f/y.restVoronoiLength(i))*curveDiff.dot(curveDiff);
      
      Hessian hess = bendEnergy.getHessian();
      for (int j=0; j<9; j++) {
        for (int k=0; k<9; k++) {
          float val = y.bendCoeff()*hess(j,k);
          CHECK_NAN(val);
          pushBackIfNotZero(*GradFx, Triplet(3*(i-1)+j, 3*(i-1)+k, val));
        }
      }
      
#ifdef NOT_DEFINED
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
      
      /// TESTING
      /*
       typedef DScalar2<float, 3, Vec3f, Mat3f> TDS;
       typedef TDS::DVector3 TVec3;
       typedef TDS::DVector2 TVec2;
       
       TVec3 teprev = TVec3(TDS(0, prevVec.x()), TDS(1, prevVec.y()), TDS(2, prevVec.z()));
       //    TVec3 teprev = TVec3(TDS(dPrevSeg.x()), TDS(dPrevSeg.y()), TDS(dPrevSeg.z()));
       TVec3 tenext = TVec3(TDS(nextVec.x()), TDS(nextVec.y()), TDS(nextVec.z()));
       //    TVec3 tenext = TVec3(TDS(0, dNextSeg.x()), TDS(1, dNextSeg.y()), TDS(2, dNextSeg.z()));
       TVec3 ttprev = teprev.normalized();
       TVec3 ttnext = tenext.normalized();
       TDS tchi = ttprev.dot(ttnext) + 1.0f;
       TVec3 tkb = TDS(2.0f)*(ttprev.cross(ttnext)) / tchi;
       TDS tk1 = TDS(0.5f)*TVec3(TDS(d2.x()), TDS(d2.y()), TDS(d2.z())).dot(tkb);
       TDS tk2 = TDS(-.5f)*TVec3(TDS(d1.x()), TDS(d1.y()), TDS(d1.z())).dot(tkb);
       //    Eigen::Matrix<float, 3, 2> t;
       //    t << tk2.getGradient(), gradK2eNext;
       //    std::cout << t << "\n\n";
       Mat3f tdiff = hessK1ePrev2 - tk1.getHessian();
       if (tdiff.cwiseAbs().maxCoeff() >= 1e-5)
         std::cout << hessK1ePrev2 << "\n" << tk1.getHessian() << "\n\n";
       // std::cout << tdiff << "\n\n";
       // std::cout << tdiff.norm() << " / " << tk1.getHessian().cwiseAbs().maxCoeff() << "\n";
      */
      ///
      
      
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
          pushBackIfNotZero(*GradFx, Triplet(3*(i-1)+j, 3*(i-1)+k, totalcoeff*block[0](j, k)));
          pushBackIfNotZero(*GradFx, Triplet(3*(i-1)+j, 3*i+k,     totalcoeff*block[1](j, k)));
          pushBackIfNotZero(*GradFx, Triplet(3*(i-1)+j, 3*(i+1)+k, totalcoeff*block[2](j, k)));
          
          pushBackIfNotZero(*GradFx, Triplet(3*i+j,     3*(i-1)+k, totalcoeff*block[3](j, k)));
          pushBackIfNotZero(*GradFx, Triplet(3*i+j,     3*i+k,     totalcoeff*block[4](j, k)));
          pushBackIfNotZero(*GradFx, Triplet(3*i+j,     3*(i+1)+k, totalcoeff*block[5](j, k)));
          
          pushBackIfNotZero(*GradFx, Triplet(3*(i+1)+j, 3*(i-1)+k, totalcoeff*block[6](j, k)));
          pushBackIfNotZero(*GradFx, Triplet(3*(i+1)+j, 3*i+k,     totalcoeff*block[7](j, k)));
          pushBackIfNotZero(*GradFx, Triplet(3*(i+1)+j, 3*(i+1)+k, totalcoeff*block[8](j, k)));
        }
      }
#endif // ifdef NOT_DEFINED
    }
    
#endif //ifdef ENABLE_BEND_AUTODIFF
  }
  
#ifdef DRAW_BENDING
  for (int i=0; i<y.numCPs(); i++) {
    ci::Vec3f f = toCi(forces.block<3,1>(3*i, 0));
    ci::Vec3f p = toCi(y.cur().points[i].pos);
    frames.push_back([f, p] (float scale) {
      ci::gl::color(0.4f, 1.0f, 1.0f);
      ci::gl::drawVector(p, p+f*scale);
    });
  }
#endif //ifdef DRAW_BENDING
  
  Profiler::stop("Bend Eval");
  return true;
}


// STRETCHING

Stretching::Stretching(const Yarn& y, EvalType et) : YarnEnergy(y, et) { }

bool Stretching::eval(VecXf* Fx, std::vector<Triplet>* GradFx, const VecXf* offset) {
  Profiler::start("Stretch Eval");
  
  for (int i=0; i<y.numSegs(); i++) {
    const Segment& seg = y.cur().segments[i];
    float restSegLength = y.rest().segments[i].length();
    Vec3f prevPoint = seg.getFirst().pos;
    Vec3f nextPoint = seg.getSecond().pos;
    
    if (offset) {
      prevPoint += offset->block<3,1>(3*i, 0);
      nextPoint += offset->block<3,1>(3*(i+1), 0);
    }
    
    Vec3f ej = nextPoint - prevPoint;
    float l = ej.norm();
    
    if (Fx) {
      Vec3f grad = ((1.0f/restSegLength) - (1.0f/l)) * ej; // Verified
      Fx->block<3,1>(3*i, 0) += y.stretchCoeff() * grad;
      Fx->block<3,1>(3*(i+1), 0) -= y.stretchCoeff() * grad;
    }
    
    if (GradFx) {
      Mat3f myHess = (Mat3f::Identity()/restSegLength -
                      (Mat3f::Identity()/l - ej * ej.transpose() / (l * l * l))); // Verified
      
      for (int j=0; j<3; j++) {
        for (int k=0; k<3; k++) {
          pushBackIfNotZero(*GradFx, Triplet(3*i+j,     3*i+k,     -y.stretchCoeff()*myHess(j,k)));
          pushBackIfNotZero(*GradFx, Triplet(3*(i+1)+j, 3*i+k,     y.stretchCoeff()*myHess(j,k)));
          pushBackIfNotZero(*GradFx, Triplet(3*i+j,     3*(i+1)+k, y.stretchCoeff()*myHess(j,k)));
          pushBackIfNotZero(*GradFx, Triplet(3*(i+1)+j, 3*(i+1)+k, -y.stretchCoeff()*myHess(j,k)));
        }
      }
    }
    
  }
  Profiler::stop("Stretch Eval");
  return true;
}



// TWISTING

Twisting::Twisting(const Yarn& y, EvalType et) : YarnEnergy(y, et) { }

bool Twisting::eval(VecXf* Fx, std::vector<Triplet>* GradFx, const VecXf* offset) {
  Profiler::start("Twist Eval");
#ifdef DRAW_TWIST
  frames.clear();
  VecXf twist = VecXf::Zero(y.numCPs()*3);
#endif // ifdef DRAW_TWIST
  for (int i=1; i<y.numCPs()-1; i++) {
    const Segment& segPrev = y.cur().segments[i-1];
    const Segment& segNext = y.cur().segments[i];
    
    Vec3f prevPoint = y.cur().points[i-1].pos;
    Vec3f curPoint  = y.cur().points[i].pos;
    Vec3f nextPoint = y.cur().points[i+1].pos;
    
    if (offset) {
      prevPoint += offset->block<3,1>(3*(i-1), 0);
      curPoint  += offset->block<3,1>(3*i,     0);
      nextPoint += offset->block<3,1>(3*(i+1), 0);
    }
    
    Vec3f ePrev = curPoint - prevPoint;
    Vec3f eNext = nextPoint - curPoint;
    float lPrev = ePrev.norm();
    float lNext = eNext.norm();
    Vec3f tPrev = ePrev.normalized();
    Vec3f tNext = eNext.normalized();
    float chi = 1.0f + tPrev.dot(tNext);
    assert(chi > 0.0f && "Segments are pointing in exactly opposite directions!");
    Vec3f curveBinorm = (2.0f * tPrev.cross(tNext)) / chi;
    float mi = segNext.getRefTwist() + (segNext.getRot() - segPrev.getRot());
    float dThetaHat = y.twistCoeff() * mi / y.restVoronoiLength(i);
    
    if (Fx) {
      Vec3f dxPrev = curveBinorm / lPrev;
      Vec3f dxNext = -curveBinorm / lNext;
      
      Fx->block<3,1>(3*(i-1), 0) += dThetaHat * dxPrev;
      Fx->block<3,1>(3*i, 0) += dThetaHat * -(dxPrev + dxNext);
      Fx->block<3,1>(3*(i+1), 0) += dThetaHat * dxNext;
      
#ifdef DRAW_TWIST
      twist.block<3,1>(3*(i-1), 0) += dThetaHat * dxPrev;
      twist.block<3,1>(3*i, 0) += dThetaHat * -(dxPrev + dxNext);
      twist.block<3,1>(3*(i+1), 0) += dThetaHat * dxNext;
#endif // ifdef DRAW_TWIST
    }
    
    if (GradFx) {
      Vec3f tTilde = (tPrev + tNext)/chi;
      Mat3f curveBinorm2 = curveBinorm * curveBinorm.transpose();
      Mat3f tPrevCross;
      tPrevCross << 0.0f, -tPrev.z(), tPrev.y(),
      tPrev.z(), 0.0f, -tPrev.x(),
      -tPrev.y(), tPrev.x(), 0.0f;
      
      Mat3f hessePrev2 = (curveBinorm2 - mi / 2.0f * (curveBinorm * (tPrev + tTilde).transpose() +
                                                      (tPrev + tTilde) * curveBinorm.transpose()))
      / lPrev / lPrev;
      Mat3f hesseNext2 = (curveBinorm2 - mi / 2.0f * (curveBinorm * (tNext + tTilde).transpose() +
                                                      (tNext + tTilde) * curveBinorm.transpose()))
      / lNext / lNext;
      Mat3f hessePreveNext = (curveBinorm2 + mi * (2.0f/chi * tPrevCross -
                                                   (curveBinorm * tTilde.transpose())));
      
      /// Populate the 3x3 block matrix:
      /// 0 1 2
      /// 3 4 5
      /// 6 7 8
      
      float coeff = -y.twistCoeff() / 2.0f / y.restVoronoiLength(i);
      Mat3f block[9];
      block[0] = coeff * hessePrev2;
      block[1] = coeff * (hessePreveNext - hessePrev2);
      block[2] = -coeff * hessePreveNext;
      block[3] = block[1].transpose();
      block[4] = coeff * (hessePrev2 + hesseNext2 - hessePreveNext - hessePreveNext.transpose());
      block[5] = coeff * (hessePreveNext - hesseNext2);
      block[6] = block[2].transpose();
      block[7] = block[5].transpose();
      block[8] = coeff * hesseNext2;
      
      for(int j=0; j<3; j++) {
        for(int k=0; k<3; k++) {
          pushBackIfNotZero(*GradFx, Triplet(3*(i-1)+j, 3*(i-1)+k, block[0](j, k)));
          pushBackIfNotZero(*GradFx, Triplet(3*(i-1)+j, 3*i+k,     block[1](j, k)));
          pushBackIfNotZero(*GradFx, Triplet(3*(i-1)+j, 3*(i+1)+k, block[2](j, k)));
          
          pushBackIfNotZero(*GradFx, Triplet(3*i+j,     3*(i-1)+k, block[3](j, k)));
          pushBackIfNotZero(*GradFx, Triplet(3*i+j,     3*i+k,     block[4](j, k)));
          pushBackIfNotZero(*GradFx, Triplet(3*i+j,     3*(i+1)+k, block[5](j, k)));
          
          pushBackIfNotZero(*GradFx, Triplet(3*(i+1)+j, 3*(i-1)+k, block[6](j, k)));
          pushBackIfNotZero(*GradFx, Triplet(3*(i+1)+j, 3*i+k,     block[7](j, k)));
          pushBackIfNotZero(*GradFx, Triplet(3*(i+1)+j, 3*(i+1)+k, block[8](j, k)));
        }
      }
    }

  }
#ifdef DRAW_TWIST
  for (int i=0; i<y.numCPs(); i++) {
    ci::Vec3f delta = toCi(twist.block<3,1>(3*i, 0));
    const Yarn* yp = &y;
    if (delta.length() > 0.0f) frames.push_back([i, delta, yp] (float scale) {
      ci::Vec3f s = toCi(yp->cur().points[i].pos);
      ci::Vec3f e = s + delta * scale;
      ci::gl::color(ci::Color(1.0f, 0.5f, 0.5f));
      ci::gl::drawVector(s, e);
    });
  }
#endif // ifdef DRAW_TWIST
  Profiler::stop("Twist Eval");
  return true;
}

// INTERNAL CONTACT

IntContact::IntContact(const Yarn& y, EvalType et) : YarnEnergy(y, et) {}

bool IntContact::eval(VecXf* Fx, std::vector<Triplet>* GradFx, const VecXf* offset) {
  Profiler::start("Int Contact Eval");
#ifdef DRAW_INT_CONTACT
  frames.clear();
#endif // ifdef DRAW_INT_CONTACT
  const float r = y.radius();
  
  for (int i=0; i<y.numSegs(); i++) {
    for (int j=i+2; j<y.numSegs(); j++) {
      if (i == 0 || j == 0 || i == y.numSegs()-1 || j == y.numSegs()-1) continue;// FIXME: eval ends
      if (i-j <= 3 && j-i <= 3) continue; // Don't evaluate close edges. FIXME: should be 1 ideally
      
      const Segment& e1 = y.cur().segments[i];
      const Segment& e2 = y.cur().segments[j];
      
      Vec3f e1p1 = e1.getFirst().pos;
      Vec3f e1p2 = e1.getSecond().pos;
      Vec3f e2p1 = e2.getFirst().pos;
      Vec3f e2p2 = e2.getSecond().pos;
      
      if (offset) {
        e1p1 += offset->block<3,1>(3*i,     0);
        e1p2 += offset->block<3,1>(3*(i+1), 0);
        e2p1 += offset->block<3,1>(3*j,     0);
        e2p2 += offset->block<3,1>(3*(j+1), 0);
      }
      
      Vec3f e1mid = (e1p1 + e1p2) / 2.0f;
      Vec3f e2mid = (e2p1 + e2p2) / 2.0f;
      
      if ((e1mid - e2mid).norm() > fmaxf(e1.length(), e2.length())) {
        continue;
      }
      
      // Splines are close, evaluate them.
      Vec3f e1p0 = y.cur().points[i-1].pos;
      Vec3f e1p3 = y.cur().points[i+2].pos;
      Vec3f e2p0 = y.cur().points[j-1].pos;
      Vec3f e2p3 = y.cur().points[j+2].pos;
      
      if (offset) {
        e1p0 += offset->block<3,1>(3*(i-1), 0);
        e1p3 += offset->block<3,1>(3*(i+2), 0);
        e2p0 += offset->block<3,1>(3*(j-1), 0);
        e2p3 += offset->block<3,1>(3*(j+2), 0);
      }
      
      Spline s1(e1p0, e1p1, e1p2, e1p3);
      Spline s2(e2p0, e2p1, e2p2, e2p3);
      
      // FIXME: hack to get an estimate of spline length.
      float l1 = e1.length();
      float l2 = e2.length();
      float coeff = contactMod * l1 * l2;
      
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
          Vec3f v = p2 - p1;
          float norm = v.norm();
          if (norm > 2*r) continue; // Quadratures are not touching, will eval to 0
          Vec4f u1(t1*t1*t1, t1*t1, t1, 1.0f);
          Vec4f u2(t2*t2*t2, t2*t2, t2, 1.0f);
          
          float dfnorm = df(norm / 2.0f / r) / 2.0f / r;
          
          if (Fx) {
            Vec3f gradBase = v/norm * dfnorm;
            
            for (int k=0; k<4; k++) {
              float c = u1.dot(constants::basis[k]);
              float d = u2.dot(constants::basis[k]);
              
              Fx->block<3,1>(3*(i-1+k), 0) += coeff * c * gradBase; // Verified
              Fx->block<3,1>(3*(j-1+k), 0) -= coeff * d * gradBase; // Verified
#ifdef DRAW_INT_CONTACT
              s1draw[k] += coeff * c * gradBase;
              s2draw[k] -= coeff * d * gradBase;
#endif // ifdef DRAW_INT_CONTACT
            }
          }
          
          if (GradFx) {
            float d2fnorm = d2f(norm / 2.0f / r) / 4.0f / r / r;
            Mat3f grad2Base = v*v.transpose()/norm/norm;
            Mat3f hessBase = (Mat3f::Identity() - grad2Base) * dfnorm / norm + grad2Base * d2fnorm;
            
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
      
      if (GradFx) {
        for (int k=0; k<4; k++) {
          for (int l=0; l<4; l++) {
            for (int p=0; p<3; p++) {
              for (int q=0; q<3; q++) {
                pushBackIfNotZero(*GradFx, Triplet(3*(i-1+k)+p, 3*(i-1+l)+q,
                                                   coeff*hess[k][l](p, q)));
                pushBackIfNotZero(*GradFx, Triplet(3*(i-1+k)+p, 3*(j-1+l)+q,
                                                   coeff*hess[k+4][l](p, q)));
                pushBackIfNotZero(*GradFx, Triplet(3*(j-1+k)+p, 3*(i-1+l)+q,
                                                   coeff*hess[k][l+4](p, q)));
                pushBackIfNotZero(*GradFx, Triplet(3*(j-1+k)+p, 3*(j-1+l)+q,
                                                   coeff*hess[k+4][l+4](p, q)));
              }
            }
          }
        }
      }
      
#ifdef DRAW_INT_CONTACT
      for (int k=0; k<4; k++) {
        ci::Vec3f s1start = toCi(y.cur().points[i-1+k].pos);
        ci::Vec3f s1vec = toCi(s1draw[k]);
        ci::Vec3f s2start = toCi(y.cur().points[j-1+k].pos);
        ci::Vec3f s2vec = toCi(s2draw[k]);
        frames.push_back([s1start, s1vec, s2start, s2vec] (float scale) {
          ci::gl::color(ci::Color(0.0f, 1.0f, 0.0f));
          ci::gl::drawVector(s1start, s1start+s1vec*scale);
          ci::gl::drawVector(s2start, s2start+s2vec*scale);
        });
      }
#endif // ifdef DRAW_INT_CONTACT
    }
  }
  
  Profiler::stop("Int Contact Eval");
  return true;
}

// PLANAR CONTACT

PlaneContact::PlaneContact(const Yarn& y, EvalType et, Vec3f normal, Vec3f origin, float stiffness)
: YarnEnergy(y, et), normal(normal.normalized()), origin(origin), stiffness(stiffness) { }

bool PlaneContact::eval(VecXf* Fx, std::vector<Triplet>* GradFx, const VecXf* offset) {
  frames.clear();
  
  for (int i=0; i<y.numCPs(); i++) {
    Vec3f pos = y.cur().points[i].pos;
    if (offset) {
      pos += offset->block<3,1>(3*i, 0);
    }
    if ((pos-origin).dot(normal) >= 0) continue;
    // else project onto the normal
    Vec3f force = (origin-pos).dot(normal)*normal;
    if (Fx) {
      Fx->block<3,1>(3*i, 0) += force * stiffness;
    }
    if (GradFx) {
      Mat3f hess = -stiffness * normal * normal.transpose();
      for (int j=0; j<3; j++) {
        for (int k=0; k<3; k++) {
          pushBackIfNotZero(*GradFx, Triplet(3*i+j, 3*i+k, hess(j, k)));
        }
      }
    }
  }
  
  return true;
}

// IMPULSE

Impulse::Impulse(const Yarn& y, EvalType et, const Clock& c, float start, float end, Vec3f force,
                 size_t index) : YarnEnergy(y, et), c(c), start(start), end(end), force(force),
index(index) { }

bool Impulse::eval(VecXf* Fx, std::vector<Triplet>* GradFx, const VecXf* offset) {
  frames.clear();
  
  if (c.time() < start || c.time() > end) return true;
  float t = sinf((c.time() - start) * constants::pi / (end - start));
  if (Fx) {
    Fx->block<3,1>(3*index, 0) += force * t;
  }
  return true;
}
