//
//  Energy.cpp
//  Visualizer
//
//  Created by eschweickart on 4/14/14.
//
//

#include "Energy.h"

void pushBackIfNotZero(std::vector<Triplet>& GradFx, Triplet& value) {
  if (value.value() != 0.0) {
    GradFx.push_back(value);
  }
}

void const RodEnergy::draw(real scale) {
  if (drawFuncs.empty()) return;
  for (std::function<void(real)> f : drawFuncs) {
    f(scale);
  }
  drawFuncs.clear();
}

// GRAVITY

Gravity::Gravity(const Rod& r, EvalType et, Vec3e& dir) : RodEnergy(r, et), dir(dir) { }

bool Gravity::eval(VecXe* Fx, std::vector<Triplet>* GradFx, const VecXe* offset) {
  if (Fx) {
    for (int i=0; i<r.numCPs(); i++) {
      Fx->segment<3>(3*i) += dir * r.getMass().diag(i);
    }
  }
  return true;
}

// SPRING

Spring::Spring(const Rod& r, EvalType et, std::size_t index, real stiffness) :
  RodEnergy(r, et), index(index), stiffness(stiffness) {}

bool Spring::eval(VecXe* Fx, std::vector<Triplet>* GradFx, const VecXe* offset) {
  // Drawing ops
#ifdef DRAW_SPRING
  Vec3c ciclamp = EtoC(clamp);
  Vec3c ciindex = EtoC(r.cur().POS(index));
  drawFuncs.push_back([ciclamp, ciindex] (real scale) {
    ci::gl::color(1.0, 0.5, 0.0);
    ci::gl::drawLine(ciclamp, ciindex);
  });
#endif // ifdef DRAW_SPRING
  
  Vec3e yPoint = r.cur().POS(index);
  if (offset) {
    yPoint += offset->segment<3>(3*index);
  }
  if (Fx) {
    Fx->segment<3>(3*index) += stiffness * (clamp - yPoint);
  }
  if (GradFx) {
    GradFx->push_back(Triplet(3*index,   3*index,   -stiffness));
    GradFx->push_back(Triplet(3*index+1, 3*index+1, -stiffness));
    GradFx->push_back(Triplet(3*index+2, 3*index+2, -stiffness));
  }
  
  return true;
}

void Spring::setClamp(Vec3e& newClamp) { clamp = newClamp; }


// MOUSE SPRING

MouseSpring::MouseSpring(const Rod& r, EvalType et, std::size_t index, real stiffness) :
  RodEnergy(r, et), index(index), stiffness(stiffness) {}

bool MouseSpring::eval(VecXe* Fx, std::vector<Triplet>* GradFx, const VecXe* offset) {
  if (!mouseDown) return true;
  assert(mouseSet && "Set the mouse position each time you call eval()!");

  Vec3e yPoint = r.cur().POS(index);
  if (offset) {
    yPoint += offset->segment<3>(3*index);
  }
  if (Fx) {
    Fx->segment<3>(3*index) += stiffness * (mouse - yPoint);
  }
  if (GradFx) {
    GradFx->push_back(Triplet(3*index,   3*index,   -stiffness));
    GradFx->push_back(Triplet(3*index+1, 3*index+1, -stiffness));
    GradFx->push_back(Triplet(3*index+2, 3*index+2, -stiffness));
  }
#ifdef DRAW_MOUSE_SPRING
  // Drawing Stuff
  Vec3c cimouse = EtoC(mouse);
  Vec3c ciindex = EtoC(yPoint);
  drawFuncs.push_back([cimouse, ciindex] (real scale) {
    ci::gl::color(1.0, 0.0, 0.0);
    ci::gl::drawLine(ciindex, cimouse);
  });
#endif // ifdef DRAW_MOUSE_SPRING
  
  return true;
}

void MouseSpring::setMouse(Vec3e& newMouse, bool newDown) {
  mouse = newMouse;
  mouseDown = newDown;
  mouseSet = true;
}


// BENDING

Bending::Bending(const Rod& r, EvalType et) : RodEnergy(r, et) { }

bool Bending::eval(VecXe* Fx, std::vector<Triplet>* GradFx, const VecXe* offset) {
  PROFILER_START("Bend Eval");
  
#ifdef DRAW_BENDING
  VecXe forces = VecXe::Zero(r.numDOF());
#endif //ifdef DRAW_BENDING
  
  for (int i=1; i<r.numCPs()-1; i++) {
    Vec3e prevPoint = r.cur().POS(i-1);
    Vec3e curPoint  = r.cur().POS(i);
    Vec3e nextPoint = r.cur().POS(i+1);
    
    if (offset) {
      prevPoint += offset->segment<3>(3*(i-1));
      curPoint  += offset->segment<3>(3*i);
      nextPoint += offset->segment<3>(3*(i+1));
    }
    
#ifdef ENABLE_BEND_AUTODIFF
    typedef Eigen::Matrix<real, 9, 1> Gradient;
    typedef Eigen::Matrix<real, 9, 9> Hessian;
    typedef DScalar2<real, 9, Gradient, Hessian> DScalar;
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
    assert(dvPrevSeg.norm() != 0.0 && dvNextSeg.norm() != 0.0 && "Edge length is 0");
    
    DVector3 dvPrevSegN = dvPrevSeg.normalized();
    DVector3 dvNextSegN = dvNextSeg.normalized();
    DScalar dotProd = dvPrevSegN.dot(dvNextSegN);
    assert(dotProd > -1.0 && "Edges are pointing in exactly opposite directions");
    
    DVector3 dvcurveBinorm = (DScalar(2.0)*dvPrevSegN.cross(dvNextSegN))/(1.0+dotProd);
    
    Vec3e prevm1 = r.cur().m1(i-1);
    Vec3e prevm2 = r.cur().m2(i-1);
    Vec3e nextm1 = r.cur().m1(i);
    Vec3e nextm2 = r.cur().m2(i);
    
    DVector3 d1prev(DScalar(prevm1.x()), DScalar(prevm1.y()), DScalar(prevm1.z()));
    DVector3 d2prev(DScalar(prevm2.x()), DScalar(prevm2.y()), DScalar(prevm2.z()));
    DVector3 d1next(DScalar(nextm1.x()), DScalar(nextm1.y()), DScalar(nextm1.z()));
    DVector3 d2next(DScalar(nextm2.x()), DScalar(nextm2.y()), DScalar(nextm2.z()));
    
    DVector2 matCurvePrev(dvcurveBinorm.dot(d2prev), -dvcurveBinorm.dot(d1prev));
    DVector2 matCurveNext(dvcurveBinorm.dot(d2next), -dvcurveBinorm.dot(d1next));
    DVector2 dvmatCurve = DScalar(0.5)*(matCurvePrev + matCurveNext);
    
    DVector2 restMatCurve(DScalar(r.restCurve(i).x()), DScalar(r.restCurve(i).y()));
    
    DVector2 curveDiff = dvmatCurve - restMatCurve;
    DVector2 fcurveDiff = r.getCS()[i].autodBend(curveDiff);
    DScalar bendEnergy = (-0.5/r.restVoronoiLength(i))*curveDiff.dot(fcurveDiff);
    
    if (Fx) {
      Fx->segment<9>(3*(i-1)) += bendEnergy.getGradient();
    }
    
    if (GradFx) {
      Hessian hess = bendEnergy.getHessian();
      for (int j=0; j<9; j++) {
        for (int k=0; k<9; k++) {
          real val = hess(j,k);
          CHECK_NAN(val);
          Triplet t(3*(i-1)+j, 3*(i-1)+k, val);
          pushBackIfNotZero(*GradFx, t);
        }
      }
    }
    
#else // ifdef ENABLE_BEND_AUTODIFF
    
    Vec3e prevVec = curPoint - prevPoint;
    Vec3e nextVec = nextPoint - curPoint;
    assert(prevVec.norm() > 0.0 && nextVec.norm() > 0.0 && "Edge length is zero!");
    
    real prevNorm = prevVec.norm();
    real nextNorm = nextVec.norm();
    Vec3e tPrev = prevVec / prevNorm;
    Vec3e tNext = nextVec / nextNorm;
    real chi = 1.0 + (tPrev.dot(tNext));
    assert(chi > 0.0 && "Edges are pointing in exactly opposite directions!");
    
    Vec3e tTilde = (tPrev + tNext)/chi;
    Vec3e d1 = r.cur().m1(i-1) + r.cur().m1(i);
    Vec3e d1tilde = d1/chi;
    Vec3e d2 = r.cur().m2(i-1) + r.cur().m2(i);
    Vec3e d2tilde = d2/chi;
    Vec3e curveBinorm = (tPrev.cross(tNext))/(0.5*chi); // Verified
    Vec2e matCurve(0.5*d2.dot(curveBinorm), -0.5*d1.dot(curveBinorm)); // Verified
    
    Vec3e gradK1ePrev = (-matCurve.x()*tTilde + tNext.cross(d2tilde)) / prevNorm; // Verified
    Vec3e gradK1eNext = (-matCurve.x()*tTilde - tPrev.cross(d2tilde)) / nextNorm; // Verified
    Vec3e gradK2ePrev = (-matCurve.y()*tTilde - tNext.cross(d1tilde)) / prevNorm; // Verified
    Vec3e gradK2eNext = (-matCurve.y()*tTilde + tPrev.cross(d1tilde)) / nextNorm; // Verified
    
    Vec2e restCurveVec = r.restCurve(i);
    // b11*2*(k1-restk1) + (b21+b12)(k2-restk2)
    // real k1coeff = 2.0 * (matCurve.x()-restCurveVec.x()); // Verified
    // b22*2*(k2-restk2) + (b21+b12)(k1-restk1)
    // real k2coeff = 2.0 * (matCurve.y()-restCurveVec.y()); // Verified
    Vec2e kcoeff = r.getCS()[i].dBend(matCurve - restCurveVec);
    real totalcoeff = 0.5 / r.restVoronoiLength(i);
    
    if (Fx) {
      Vec3e gradePrev = totalcoeff * (gradK1ePrev * kcoeff.x() + gradK2ePrev * kcoeff.y()); // Verified
      Vec3e gradeNext = totalcoeff * (gradK1eNext * kcoeff.x() + gradK2eNext * kcoeff.y()); // Verified
      
      Fx->segment<3>(3*(i-1)) += gradePrev;
      Fx->segment<3>(3*i)     += (gradeNext - gradePrev);
      Fx->segment<3>(3*(i+1)) += -gradeNext;
      
#ifdef DRAW_BENDING
      forces.segment<3>(3*(i-1)) += gradePrev;
      forces.segment<3>(3*i)     += (gradeNext - gradePrev);
      forces.segment<3>(3*(i+1)) += -gradeNext;
#endif //ifdef DRAW_BENDING
    }
    
    if (GradFx) {
      typedef Eigen::Matrix<real, 9, 1> Gradient;
      typedef Eigen::Matrix<real, 9, 9> Hessian;
      typedef DScalar2<real, 9, Gradient, Hessian> DScalar;
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
      assert(dvPrevSeg.norm() > 0.0 && dvNextSeg.norm() > 0.0 && "Edge length is 0");
      
      DVector3 dvPrevSegN = dvPrevSeg.normalized();
      DVector3 dvNextSegN = dvNextSeg.normalized();
      DScalar dotProd = dvPrevSegN.dot(dvNextSegN);
      assert(dotProd > -1.0 && "Edges are pointing in exactly opposite directions");
      
      DVector3 dvcurveBinorm = (DScalar(2.0)*dvPrevSegN.cross(dvNextSegN))/(1.0+dotProd);
      
      Vec3e prevm1 = r.cur().m1(i-1);
      Vec3e prevm2 = r.cur().m2(i-1);
      Vec3e nextm1 = r.cur().m1(i);
      Vec3e nextm2 = r.cur().m2(i);
      
      DVector3 d1prev(DScalar(prevm1.x()), DScalar(prevm1.y()), DScalar(prevm1.z()));
      DVector3 d2prev(DScalar(prevm2.x()), DScalar(prevm2.y()), DScalar(prevm2.z()));
      DVector3 d1next(DScalar(nextm1.x()), DScalar(nextm1.y()), DScalar(nextm1.z()));
      DVector3 d2next(DScalar(nextm2.x()), DScalar(nextm2.y()), DScalar(nextm2.z()));
      
      DVector2 matCurvePrev(dvcurveBinorm.dot(d2prev), -dvcurveBinorm.dot(d1prev));
      DVector2 matCurveNext(dvcurveBinorm.dot(d2next), -dvcurveBinorm.dot(d1next));
      DVector2 dvmatCurve = DScalar(0.5)*(matCurvePrev + matCurveNext);
      
      DVector2 restMatCurve(DScalar(r.restCurve(i).x()), DScalar(r.restCurve(i).y()));
      
      DVector2 curveDiff = dvmatCurve - restMatCurve;
      DVector2 fcurveDiff = r.getCS()[i].autodBend(curveDiff);
      DScalar bendEnergy = (-0.5/r.restVoronoiLength(i))*curveDiff.dot(fcurveDiff);
      
      Hessian hess = bendEnergy.getHessian();
      for (int j=0; j<9; j++) {
        for (int k=0; k<9; k++) {
          real val = hess(j,k);
          CHECK_NAN(val);
          Triplet t(3*(i-1)+j, 3*(i-1)+k, val);
          pushBackIfNotZero(*GradFx, t);
        }
      }
      
#ifdef NOT_DEFINED
      Mat3e tTilde2 = tTilde*tTilde.transpose();
      
      Mat3e tNextxd2TildextTilde = (tNext.cross(d2tilde))*tTilde.transpose();
      Mat3e tPrevxd2TildextTilde = (tPrev.cross(d2tilde))*tTilde.transpose();
      Mat3e tNextxd1TildextTilde = (tNext.cross(d1tilde))*tTilde.transpose();
      Mat3e tPrevxd1TildextTilde = (tPrev.cross(d1tilde))*tTilde.transpose();
      
      Mat3e d2TildeCross, d1TildeCross;
      d2TildeCross << 0.0, -d2tilde.z(), d2tilde.y(),
      d2tilde.z(), 0.0, -d2tilde.x(),
      -d2tilde.y(), d2tilde.x(), 0.0;
      d1TildeCross << 0.0, -d1tilde.z(), d1tilde.y(),
      d1tilde.z(), 0.0, -d1tilde.x(),
      -d1tilde.y(), d1tilde.x(), 0.0;
      
      Mat3e hessK1ePrev2 = 2.0*matCurve.x()*tTilde2-tNextxd2TildextTilde-tNextxd2TildextTilde.transpose();
      hessK1ePrev2 -= (matCurve.x()/chi)*(Mat3e::Identity() - (tPrev*tPrev.transpose()));
      hessK1ePrev2 += 0.25*(curveBinorm*prevSeg.m2().transpose()+prevSeg.m2()*curveBinorm.transpose());
      hessK1ePrev2 /= prevVec.dot(prevVec);
      
      Mat3e hessK2ePrev2 = 2.0*matCurve.y()*tTilde2+tNextxd1TildextTilde+tNextxd1TildextTilde.transpose();
      hessK2ePrev2 -= (matCurve.y()/chi)*(Mat3e::Identity() - (tPrev*tPrev.transpose()));
      hessK2ePrev2 += 0.25*(curveBinorm*prevSeg.m1().transpose()+prevSeg.m1()*curveBinorm.transpose());
      hessK2ePrev2 /= prevVec.dot(prevVec);
      
      Mat3e hessK1eNext2 = 2.0*matCurve.x()*tTilde2+tPrevxd2TildextTilde+tPrevxd2TildextTilde.transpose();
      hessK1eNext2 -= (matCurve.x()/chi)*(Mat3e::Identity() - (tNext*tNext.transpose()));
      hessK1eNext2 += 0.25*(curveBinorm*nextSeg.m2().transpose()+nextSeg.m2()*curveBinorm.transpose());
      hessK1eNext2 /= nextVec.dot(nextVec);
      
      Mat3e hessK2eNext2 = 2.0*matCurve.y()*tTilde2-tPrevxd1TildextTilde-tPrevxd1TildextTilde.transpose();
      hessK2eNext2 -= (matCurve.y()/chi)*(Mat3e::Identity() - (tNext*tNext.transpose()));
      hessK2eNext2 += 0.25*(curveBinorm*nextSeg.m1().transpose()+nextSeg.m1()*curveBinorm.transpose());
      hessK2eNext2 /= nextVec.dot(nextVec);
      
      Mat3e hessK1ePreveNext = (-matCurve.x()/chi)*(Mat3e::Identity() + (tPrev * tNext.transpose()));
      hessK1ePreveNext += (2.0*matCurve.x()*tTilde2) - tNextxd2TildextTilde + tPrevxd2TildextTilde.transpose();
      hessK1ePreveNext -= d2TildeCross;
      hessK1ePreveNext /= (nextVec.norm() * prevVec.norm());
      
      Mat3e hessK2ePreveNext = (-matCurve.y()/chi)*(Mat3e::Identity() + (tPrev * tNext.transpose()));
      hessK2ePreveNext += (2.0*matCurve.y()*tTilde2) - tNextxd1TildextTilde + tPrevxd1TildextTilde.transpose();
      hessK2ePreveNext -= d1TildeCross;
      hessK2ePreveNext /= (nextVec.norm() * prevVec.norm());
      
      /// TESTING
      /*
       typedef DScalar2<real, 3, Vec3e, Mat3e> TDS;
       typedef TDS::DVector3 TVec3;
       typedef TDS::DVector2 TVec2;
       
       TVec3 teprev = TVec3(TDS(0, prevVec.x()), TDS(1, prevVec.y()), TDS(2, prevVec.z()));
       //    TVec3 teprev = TVec3(TDS(dPrevSeg.x()), TDS(dPrevSeg.y()), TDS(dPrevSeg.z()));
       TVec3 tenext = TVec3(TDS(nextVec.x()), TDS(nextVec.y()), TDS(nextVec.z()));
       //    TVec3 tenext = TVec3(TDS(0, dNextSeg.x()), TDS(1, dNextSeg.y()), TDS(2, dNextSeg.z()));
       TVec3 ttprev = teprev.normalized();
       TVec3 ttnext = tenext.normalized();
       TDS tchi = ttprev.dot(ttnext) + 1.0;
       TVec3 tkb = TDS(2.0)*(ttprev.cross(ttnext)) / tchi;
       TDS tk1 = TDS(0.5)*TVec3(TDS(d2.x()), TDS(d2.y()), TDS(d2.z())).dot(tkb);
       TDS tk2 = TDS(-.5)*TVec3(TDS(d1.x()), TDS(d1.y()), TDS(d1.z())).dot(tkb);
       //    Eigen::Matrix<real, 3, 2> t;
       //    t << tk2.getGradient(), gradK2eNext;
       //    std::cout << t << "\n\n";
       Mat3e tdiff = hessK1ePrev2 - tk1.getHessian();
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
      
      Mat3e block[9];
      
      block[0] = k1coeff*hessK1ePrev2 + k2coeff*hessK2ePrev2 + 2.0*(gradK1ePrev*gradK1ePrev.transpose() + gradK2ePrev*gradK2ePrev.transpose());
      
      block[1] = k1coeff*(hessK1ePreveNext-hessK1ePrev2) + k2coeff*(hessK2ePreveNext-hessK2ePrev2)
      + 2.0*(-gradK1ePrev*(gradK1ePrev-gradK1eNext).transpose() + -gradK2ePrev*(gradK2ePrev-gradK2eNext).transpose());
      
      block[2] = -k1coeff*hessK1ePreveNext - k2coeff*hessK2ePreveNext + 2.0*(-gradK1ePrev*gradK1eNext.transpose() + -gradK2ePrev*gradK2eNext.transpose());
      
      block[3] = k1coeff*(hessK1ePreveNext.transpose() - hessK1ePrev2) + k2coeff*(hessK2ePreveNext.transpose() - hessK2ePrev2) + 2.0*((gradK1ePrev-gradK1eNext)*-gradK1ePrev.transpose() + (gradK2ePrev-gradK2eNext)*-gradK2ePrev.transpose());
      
      block[4] = k1coeff*(hessK1ePrev2 + hessK1eNext2 - hessK1ePreveNext - hessK1ePreveNext.transpose()) +
      k2coeff*(hessK2ePrev2 + hessK2eNext2 - hessK2ePreveNext - hessK2ePreveNext.transpose()) +
      2.0*((gradK1ePrev-gradK1eNext)*(gradK1ePrev-gradK1eNext).transpose() + (gradK2ePrev-gradK2eNext)*(gradK2ePrev-gradK2eNext).transpose());
      
      block[5] = k1coeff*(hessK1ePreveNext - hessK1eNext2) + k2coeff*(hessK2ePreveNext - hessK2eNext2)
      + 2.0*((gradK1ePrev-gradK1eNext)*gradK1eNext.transpose() + (gradK2ePrev-gradK2eNext)*gradK2eNext.transpose());
      
      block[6] = -k1coeff*hessK1ePreveNext.transpose() - k2coeff*hessK2ePreveNext.transpose()
      + 2.0*(gradK1eNext*-gradK1ePrev.transpose() + gradK2eNext*-gradK2ePrev.transpose());
      
      block[7] = k1coeff*(hessK1ePreveNext.transpose()-hessK1eNext2) + k2coeff*(hessK2ePreveNext.transpose()-hessK2eNext2) + 2.0*(gradK1eNext*(gradK1ePrev-gradK1eNext).transpose() + gradK2eNext*(gradK2ePrev-gradK2eNext).transpose());
      
      block[8] = k1coeff*hessK1eNext2 + k2coeff*hessK2eNext2 + 2.0*(gradK1eNext*gradK1eNext.transpose() + gradK2eNext*gradK2eNext.transpose());
      
      
      for(int j=0; j<3; j++) {
        for(int k=0; k<3; k++) {
          Triplet t(3*(i-1)+j, 3*(i-1)+k, totalcoeff*block[0](j, k));
          pushBackIfNotZero(*GradFx, t);
          t = Triplet(3*(i-1)+j, 3*i+k, totalcoeff*block[1](j, k));
          pushBackIfNotZero(*GradFx, t);
          t = Triplet(3*(i-1)+j, 3*(i+1)+k, totalcoeff*block[2](j, k));
          pushBackIfNotZero(*GradFx, t);
          
          t = Triplet(3*i+j, 3*(i-1)+k, totalcoeff*block[3](j, k));
          pushBackIfNotZero(*GradFx, t);
          t = Triplet(3*i+j, 3*i+k, totalcoeff*block[4](j, k));
          pushBackIfNotZero(*GradFx, t);
          t = Triplet(3*i+j, 3*(i+1)+k, totalcoeff*block[5](j, k));
          pushBackIfNotZero(*GradFx, t);
          
          t = Triplet(3*(i+1)+j, 3*(i-1)+k, totalcoeff*block[6](j, k));
          pushBackIfNotZero(*GradFx, t);
          t = Triplet(3*(i+1)+j, 3*i+k, totalcoeff*block[7](j, k));
          pushBackIfNotZero(*GradFx, t);
          t = Triplet(3*(i+1)+j, 3*(i+1)+k, totalcoeff*block[8](j, k));
          pushBackIfNotZero(*GradFx, t);
        }
      }
#endif // ifdef NOT_DEFINED
    }
    
#endif //ifdef ENABLE_BEND_AUTODIFF
  }
  
#ifdef DRAW_BENDING
  for (int i=0; i<r.numCPs(); i++) {
    Vec3c f = EtoC(forces.segment<3>(3*i));
    Vec3c p = EtoC(r.cur().POS(i));
    drawFuncs.push_back([f, p] (real scale) {
      ci::gl::color(0.4, 1.0, 1.0);
      ci::gl::drawLine(p, p+f*scale);
    });
  }
#endif //ifdef DRAW_BENDING
  
  PROFILER_STOP("Bend Eval");
  return true;
}


// STRETCHING

Stretching::Stretching(const Rod& r, EvalType et) : RodEnergy(r, et) { }

bool Stretching::eval(VecXe* Fx, std::vector<Triplet>* GradFx, const VecXe* offset) {
  PROFILER_START("Stretch Eval");
  
  for (int i=0; i<r.numEdges(); i++) {
    real restEdgeLength = r.rest().edgeLength(i);
    Vec3e prevPoint = r.cur().POS(i);
    Vec3e nextPoint = r.cur().POS(i+1);
    
    if (offset) {
      prevPoint += offset->segment<3>(3*i);
      nextPoint += offset->segment<3>(3*(i+1));
    }
    
    Vec3e ej = nextPoint - prevPoint;
    real l = ej.norm();
    
    real stretchCoeff = r.getCS()[i].stretchCoeff() + r.getCS()[i+1].stretchCoeff() / 2.0;
    
    if (Fx) {
      Vec3e grad = ((1.0/restEdgeLength) - (1.0/l)) * ej; // Verified
      Fx->segment<3>(3*i) += stretchCoeff * grad;
      Fx->segment<3>(3*(i+1)) -= stretchCoeff * grad;
    }
    
    if (GradFx) {
      Mat3e myHess = (Mat3e::Identity()/restEdgeLength -
                      (Mat3e::Identity()/l - ej * ej.transpose() / (l * l * l))); // Verified
      
      for (int j=0; j<3; j++) {
        for (int k=0; k<3; k++) {
          Triplet t(3*i+j, 3*i+k, -stretchCoeff * myHess(j,k));
          pushBackIfNotZero(*GradFx, t);
          t = Triplet(3*(i+1)+j, 3*i+k, stretchCoeff * myHess(j,k));
          pushBackIfNotZero(*GradFx, t);
          t = Triplet(3*i+j, 3*(i+1)+k, stretchCoeff * myHess(j,k));
          pushBackIfNotZero(*GradFx, t);
          t = Triplet(3*(i+1)+j, 3*(i+1)+k, -stretchCoeff * myHess(j,k));
          pushBackIfNotZero(*GradFx, t);
        }
      }
    }
    
  }
  PROFILER_STOP("Stretch Eval");
  return true;
}



// TWISTING

Twisting::Twisting(const Rod& r, EvalType et) : RodEnergy(r, et) { }

bool Twisting::eval(VecXe* Fx, std::vector<Triplet>* GradFx, const VecXe* offset) {
  PROFILER_START("Twist Eval");
#ifdef DRAW_TWIST
  VecXe twist = VecXe::Zero(r.numCPs()*3);
#endif // ifdef DRAW_TWIST
  for (int i=1; i<r.numCPs()-1; i++) {
    Vec3e prevPoint = r.cur().POS(i-1);
    Vec3e curPoint  = r.cur().POS(i);
    Vec3e nextPoint = r.cur().POS(i+1);
    
    if (offset) {
      prevPoint += offset->segment<3>(3*(i-1));
      curPoint  += offset->segment<3>(3*i);
      nextPoint += offset->segment<3>(3*(i+1));
    }
    
    Vec3e ePrev = curPoint - prevPoint;
    Vec3e eNext = nextPoint - curPoint;
    real lPrev = ePrev.norm();
    real lNext = eNext.norm();
    Vec3e tPrev = ePrev.normalized();
    Vec3e tNext = eNext.normalized();
    real chi = 1.0 + tPrev.dot(tNext);
    assert(chi > 0.0 && "Edges are pointing in exactly opposite directions!");
    Vec3e curveBinorm = (2.0 * tPrev.cross(tNext)) / chi;
    real mi = r.cur().refTwist(i) + (r.cur().rot(i) - r.cur().rot(i-1));
    real dThetaHat = r.getCS()[i].twistCoeff() * mi / r.restVoronoiLength(i);
    
    if (Fx) {
      Vec3e dxPrev = curveBinorm / lPrev;
      Vec3e dxNext = -curveBinorm / lNext;
      
      Fx->segment<3>(3*(i-1)) += dThetaHat * dxPrev;
      Fx->segment<3>(3*i) += dThetaHat * -(dxPrev + dxNext);
      Fx->segment<3>(3*(i+1)) += dThetaHat * dxNext;
      
#ifdef DRAW_TWIST
      twist.segment<3>(3*(i-1)) += dThetaHat * dxPrev;
      twist.segment<3>(3*i) += dThetaHat * -(dxPrev + dxNext);
      twist.segment<3>(3*(i+1)) += dThetaHat * dxNext;
#endif // ifdef DRAW_TWIST
    }
    
    if (GradFx) {
      Vec3e tTilde = (tPrev + tNext)/chi;
      Mat3e curveBinorm2 = curveBinorm * curveBinorm.transpose();
      Mat3e tPrevCross;
      tPrevCross << 0.0, -tPrev.z(), tPrev.y(),
      tPrev.z(), 0.0, -tPrev.x(),
      -tPrev.y(), tPrev.x(), 0.0;
      
      Mat3e hessePrev2 = (curveBinorm2 - mi / 2.0 * (curveBinorm * (tPrev + tTilde).transpose() +
                                                      (tPrev + tTilde) * curveBinorm.transpose()))
      / lPrev / lPrev;
      Mat3e hesseNext2 = (curveBinorm2 - mi / 2.0 * (curveBinorm * (tNext + tTilde).transpose() +
                                                      (tNext + tTilde) * curveBinorm.transpose()))
      / lNext / lNext;
      Mat3e hessePreveNext = (curveBinorm2 + mi * (2.0/chi * tPrevCross -
                                                   (curveBinorm * tTilde.transpose())));
      
      /// Populate the 3x3 block matrix:
      /// 0 1 2
      /// 3 4 5
      /// 6 7 8
      
      real coeff = -r.getCS()[i].twistCoeff() / 2.0 / r.restVoronoiLength(i);
      Mat3e block[9];
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
          Triplet t(3*(i-1)+j, 3*(i-1)+k, block[0](j, k));
          pushBackIfNotZero(*GradFx, t);
          t = Triplet(3*(i-1)+j, 3*i+k, block[1](j, k));
          pushBackIfNotZero(*GradFx, t);
          t = Triplet(3*(i-1)+j, 3*(i+1)+k, block[2](j, k));
          pushBackIfNotZero(*GradFx, t);
          
          t = Triplet(3*i+j, 3*(i-1)+k, block[3](j, k));
          pushBackIfNotZero(*GradFx, t);
          t = Triplet(3*i+j, 3*i+k, block[4](j, k));
          pushBackIfNotZero(*GradFx, t);
          t = Triplet(3*i+j, 3*(i+1)+k, block[5](j, k));
          pushBackIfNotZero(*GradFx, t);
          
          t = Triplet(3*(i+1)+j, 3*(i-1)+k, block[6](j, k));
          pushBackIfNotZero(*GradFx, t);
          t = Triplet(3*(i+1)+j, 3*i+k,     block[7](j, k));
          pushBackIfNotZero(*GradFx, t);
          t = Triplet(3*(i+1)+j, 3*(i+1)+k, block[8](j, k));
          pushBackIfNotZero(*GradFx, t);
        }
      }
    }

  }
#ifdef DRAW_TWIST
  for (int i=0; i<r.numCPs(); i++) {
    Vec3c delta = EtoC(twist.segment<3>(3*i));
    const Rod* rp = &r;
    if (delta.length() > 0.0) drawFuncs.push_back([i, delta, rp] (real scale) {
      Vec3c s = EtoC(rp->cur().POS(i));
      Vec3c e = s + delta * scale;
      ci::gl::color(ci::Color(1.0, 0.5, 0.5));
      ci::gl::drawLine(s, e);
    });
  }
#endif // ifdef DRAW_TWIST
  PROFILER_STOP("Twist Eval");
  return true;
}

// INTERNAL CONTACT

IntContact::IntContact(const Rod& r, EvalType et) : RodEnergy(r, et) {}

bool IntContact::eval(VecXe* Fx, std::vector<Triplet>* GradFx, const VecXe* offset) {
  PROFILER_START("Int Contact Eval");
  
  for (int i=0; i<r.numEdges(); i++) {
    for (int j=i+2; j<r.numEdges(); j++) {
      if (i == 0 || j == 0 || i == r.numEdges()-1 || j == r.numEdges()-1) continue;// FIXME: eval ends
      if (i-j <= 3 && j-i <= 3) continue; // Don't evaluate close edges. FIXME: should be 1 ideally
      
      const real rad = (r.getCS()[i].approxRadius() + r.getCS()[i+1].approxRadius()
                      + r.getCS()[j].approxRadius() + r.getCS()[j+1].approxRadius()) / 4.0;
      
      Vec3e e1p1 = r.cur().POS(i);
      Vec3e e1p2 = r.cur().POS(i+1);
      Vec3e e2p1 = r.cur().POS(j);
      Vec3e e2p2 = r.cur().POS(j+1);
      
      if (offset) {
        e1p1 += offset->segment<3>(3*i);
        e1p2 += offset->segment<3>(3*(i+1));
        e2p1 += offset->segment<3>(3*j);
        e2p2 += offset->segment<3>(3*(j+1));
      }
      
      Vec3e e1mid = (e1p1 + e1p2) / 2.0;
      Vec3e e2mid = (e2p1 + e2p2) / 2.0;
      
      // FIXME: hack to get an estimate of spline length.
      real l1 = (e1p2 - e1p1).norm();
      real l2 = (e2p2 - e2p1).norm();
      
      if ((e1mid - e2mid).norm() > fmaxf(l1, l2)) {
        continue;
      }
      
      // Splines are close, evaluate them.
      Vec3e e1p0 = r.cur().POS(i-1);
      Vec3e e1p3 = r.cur().POS(i+2);
      Vec3e e2p0 = r.cur().POS(j-1);
      Vec3e e2p3 = r.cur().POS(j+2);
      
      if (offset) {
        e1p0 += offset->segment<3>(3*(i-1));
        e1p3 += offset->segment<3>(3*(i+2));
        e2p0 += offset->segment<3>(3*(j-1));
        e2p3 += offset->segment<3>(3*(j+2));
      }
      
      Spline s1(e1p0, e1p1, e1p2, e1p3);
      Spline s2(e2p0, e2p1, e2p2, e2p3);
      
      real coeff = contactMod * l1 * l2;
      
      Mat3e hess[8][8];
      if (et == Implicit && GradFx) {
        for (int k=0; k<8; k++) {
          for (int l=0; l<8; l++) {
            hess[k][l] = Mat3e::Zero();
          }
        }
      }
      
#ifdef DRAW_INT_CONTACT
      Vec3e s1draw[4];
      Vec3e s2draw[4];
      for (int k=0; k<4; k++) {
        s1draw[k] = Vec3e::Zero();
        s2draw[k] = Vec3e::Zero();
      }
#endif // ifdef DRAW_INT_CONTACT
      
      for (int n=0; n<=nb; n++) {
        for (int m=0; m<=nb; m++) {
          real t1 = ((real) n) / nb;
          real t2 = ((real) m) / nb;
          Vec3e p1 = s1.eval(t1);
          Vec3e p2 = s2.eval(t2);
          Vec3e v = p2 - p1;
          real norm = v.norm();
          if (norm > 2*rad) continue; // Quadratures are not touching, will eval to 0
          Vec4e u1(t1*t1*t1, t1*t1, t1, 1.0);
          Vec4e u2(t2*t2*t2, t2*t2, t2, 1.0);
          
          real dfnorm = df(norm / 2.0 / rad) / 2.0 / rad;
          
          if (Fx) {
            Vec3e gradBase = v/norm * dfnorm;
            
            for (int k=0; k<4; k++) {
              real c = u1.dot(constants::basis[k]);
              real d = u2.dot(constants::basis[k]);
              
              Fx->segment<3>(3*(i-1+k)) += coeff * c * gradBase; // Verified
              Fx->segment<3>(3*(j-1+k)) -= coeff * d * gradBase; // Verified
#ifdef DRAW_INT_CONTACT
              s1draw[k] += coeff * c * gradBase;
              s2draw[k] -= coeff * d * gradBase;
#endif // ifdef DRAW_INT_CONTACT
            }
          }
          
          if (GradFx) {
            real d2fnorm = d2f(norm / 2.0 / rad) / 4.0 / rad / rad;
            Mat3e grad2Base = v*v.transpose()/norm/norm;
            Mat3e hessBase = (Mat3e::Identity() - grad2Base) * dfnorm / norm + grad2Base * d2fnorm;
            
            for (int k=0; k<4; k++) {
              for (int l=0; l<4; l++) {
                real ck = u1.dot(constants::basis[k]);
                real cl = u1.dot(constants::basis[l]);
                real dk = u2.dot(constants::basis[k]);
                real dl = u2.dot(constants::basis[l]);
                
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
                Triplet t(3*(i-1+k)+p, 3*(i-1+l)+q, coeff*hess[k][l](p, q));
                pushBackIfNotZero(*GradFx, t);
                t = Triplet(3*(i-1+k)+p, 3*(j-1+l)+q, coeff*hess[k+4][l](p, q));
                pushBackIfNotZero(*GradFx, t);
                t = Triplet(3*(j-1+k)+p, 3*(i-1+l)+q, coeff*hess[k][l+4](p, q));
                pushBackIfNotZero(*GradFx, t);
                t = Triplet(3*(j-1+k)+p, 3*(j-1+l)+q, coeff*hess[k+4][l+4](p, q));
                pushBackIfNotZero(*GradFx, t);
              }
            }
          }
        }
      }
      
#ifdef DRAW_INT_CONTACT
      for (int k=0; k<4; k++) {
        Vec3c s1start = EtoC(r.cur().POS(i-1+k));
        Vec3c s1vec = EtoC(s1draw[k]);
        Vec3c s2start = EtoC(r.cur().POS(j-1+k));
        Vec3c s2vec = EtoC(s2draw[k]);
        drawFuncs.push_back([s1start, s1vec, s2start, s2vec] (real scale) {
          ci::gl::color(ci::Color(0.0, 1.0, 0.0));
          ci::gl::drawLine(s1start, s1start+s1vec*scale);
          ci::gl::drawLine(s2start, s2start+s2vec*scale);
        });
      }
#endif // ifdef DRAW_INT_CONTACT
    }
  }
  
  PROFILER_STOP("Int Contact Eval");
  return true;
}

// PLANAR CONTACT

PlaneContact::PlaneContact(const Rod& r, EvalType et, Vec3e& normal, Vec3e& origin, real stiffness)
: RodEnergy(r, et), normal(normal.normalized()), origin(origin), stiffness(stiffness) { }

bool PlaneContact::eval(VecXe* Fx, std::vector<Triplet>* GradFx, const VecXe* offset) {
  for (int i=0; i<r.numCPs(); i++) {
    Vec3e pos = r.cur().POS(i);
    if (offset) {
      pos += offset->segment<3>(3*i);
    }
    if ((pos-origin).dot(normal) >= 0) continue;
    // else project onto the normal
    Vec3e force = (origin-pos).dot(normal)*normal;
    if (Fx) {
      Fx->segment<3>(3*i) += force * stiffness;
    }
    if (GradFx) {
      Mat3e hess = -stiffness * normal * normal.transpose();
      for (int j=0; j<3; j++) {
        for (int k=0; k<3; k++) {
          Triplet t(3*i+j, 3*i+k, hess(j, k));
          pushBackIfNotZero(*GradFx, t);
        }
      }
    }
  }
  
  return true;
}

// IMPULSE

Impulse::Impulse(const Rod& r, EvalType et, const Clock& c, real start, real end, Vec3e& force,
                 std::size_t index) : RodEnergy(r, et), c(c), start(start), end(end), force(force),
index(index) { }

bool Impulse::eval(VecXe* Fx, std::vector<Triplet>* GradFx, const VecXe* offset) {
  if (c.time() < start || c.time() > end) return true;
  real t = sinf((c.time() - start) * constants::pi / (end - start));
  if (Fx) {
    Fx->segment<3>(3*index) += force * t;
  }
  return true;
}

FEMBending::FEMBending(const Rod& r, EvalType et) : RodEnergy(r, et) {
  std::size_t n = r.numCPs();
  real l = (r.rest().POS(0) - r.rest().POS(n-1)).norm();
  real h = (r.rest().POS(0) - r.rest().POS(1)).norm();
  real modmu2 = r.youngsModulus * r.getCS()[0].areaMoment()(0, 0) / (n * l * l * l * h * h * h * h);
  
  modDxxxxTriplets.push_back(Triplet(0, 0, -1.0 * modmu2));
  modDxxxxTriplets.push_back(Triplet(0, 1, 2.0 * modmu2));
  modDxxxxTriplets.push_back(Triplet(0, 2, -1.0 * modmu2));
  modDxxxxTriplets.push_back(Triplet(1, 0, 2.0 * modmu2));
  modDxxxxTriplets.push_back(Triplet(1, 1, -5.0 * modmu2));
  modDxxxxTriplets.push_back(Triplet(1, 2, 4.0 * modmu2));
  modDxxxxTriplets.push_back(Triplet(1, 3, -1.0 * modmu2));
  
  for (int i=2; i<n-2; i++) {
    modDxxxxTriplets.push_back(Triplet(i, i-2, -1.0 * modmu2));
    modDxxxxTriplets.push_back(Triplet(i, i-1, 4.0 * modmu2));
    modDxxxxTriplets.push_back(Triplet(i, i, -6.0 * modmu2));
    modDxxxxTriplets.push_back(Triplet(i, i+1, 4.0 * modmu2));
    modDxxxxTriplets.push_back(Triplet(i, i+2, -1.0 * modmu2));
  }
  
  modDxxxxTriplets.push_back(Triplet(n-2, n-4, -1.0 * modmu2));
  modDxxxxTriplets.push_back(Triplet(n-2, n-3, 4.0 * modmu2));
  modDxxxxTriplets.push_back(Triplet(n-2, n-2, -5.0 * modmu2));
  modDxxxxTriplets.push_back(Triplet(n-2, n-1, 2.0 * modmu2));
  modDxxxxTriplets.push_back(Triplet(n-1, n-3, -1.0 * modmu2));
  modDxxxxTriplets.push_back(Triplet(n-1, n-2, 2.0 * modmu2));
  modDxxxxTriplets.push_back(Triplet(n-1, n-1, -1.0 * modmu2));
  
  modDxxxx.resize(n, n);
  modDxxxx.setFromTriplets(modDxxxxTriplets.begin(), modDxxxxTriplets.end());

  std::vector<Triplet> triplets;
  X.resize(n, n*3);
  for (int i=0; i<r.numCPs(); i++) {
    triplets.push_back(Triplet(i, 3*i, 1.0));
  }
  X.setFromTriplets(triplets.begin(), triplets.end());
}

bool FEMBending::eval(VecXe* Fx, std::vector<Triplet>* GradFx, const VecXe* offset) {
  if (Fx) {
    VecXe x = X * r.cur().pos;
    if (offset) {
      x += X * (*offset);
    }
    (*Fx) += X.transpose() * (modDxxxx * x);
  }
  
  if (GradFx) {
    for (Triplet t : modDxxxxTriplets) {
      GradFx->push_back(Triplet(t.row()*3, t.col()*3, t.value()));
    }
  }
  return true;
}
