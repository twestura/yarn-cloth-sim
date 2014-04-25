//
//  Energy.cpp
//  Visualizer
//
//  Created by eschweickart on 4/14/14.
//
//

#include "Energy.h"


void pushBackIfNotZero(std::vector<Triplet>& GradFx, Triplet value) {
  if (value.value() != 0) {
    GradFx.push_back(value);
  }
}

// GRAVITY

Gravity::Gravity(const Yarn& y, EvalType et, Vec3f dir) : YarnEnergy(y, et), dir(dir) {}

bool Gravity::eval(VecXf& Fx, std::vector<Triplet>& GradFx, const VecXf& dqdot, Clock& c) {
  assert(et == Explicit && "Unsupported EvalType");
  for (int i=0; i<y.numCPs(); i++) {
    Fx.block<3,1>(3*i, 0) -= dir * c.timestep();
  }
  return true;
}


// SPRING

Spring::Spring(const Yarn& y, EvalType et, size_t index, float stiffness) :
  YarnEnergy(y, et), index(index), stiffness(stiffness) {}

bool Spring::eval(VecXf& Fx, std::vector<Triplet>& GradFx, const VecXf& dqdot, Clock& c) {
  const float h = c.timestep();
  if (et == Explicit) {
    Fx.block<3,1>(3*index, 0) -= (clamp - y.cur().points[index].pos) * stiffness * h;
  } else {
    Fx.block<3,1>(3*index, 0) -= (clamp - (y.cur().points[index].pos +
                                           h*(y.cur().points[index].vel
                                              + dqdot.block<3,1>(3*index, 0)))) * stiffness * h;
    
    GradFx.push_back(Triplet(3*index,   3*index,   h*h*stiffness));
    GradFx.push_back(Triplet(3*index+1, 3*index+1, h*h*stiffness));
    GradFx.push_back(Triplet(3*index+2, 3*index+2, h*h*stiffness));
  }
  
  return true;
}

void Spring::setClamp(Vec3f newClamp) { clamp = newClamp; }


// MOUSE SPRING

MouseSpring::MouseSpring(const Yarn& y, EvalType et, size_t index, float stiffness) :
  YarnEnergy(y, et), index(index), stiffness(stiffness) {}

bool MouseSpring::eval(VecXf& Fx, std::vector<Triplet>& GradFx, const VecXf& dqdot, Clock& c) {
  if (!mouseDown) return true;
  assert(et == Explicit && "Unsupported EvalType");
  assert(mouseSet && "Set the mouse position each time you call eval()!");
  Fx.block<3,1>(3*index, 0) -= (mouse - y.cur().points[index].pos) * stiffness * c.timestep();
  return true;
}

void MouseSpring::setMouse(Vec3f newMouse, bool newDown) {
  mouse = newMouse;
  mouseDown = newDown;
  mouseSet = true;
}


// BENDING

Bending::Bending(const Yarn& y, EvalType et) : YarnEnergy(y, et) {
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

//#define ENABLE_AUTODIFF
#define NUM_VARS 9
bool Bending::eval(VecXf& Fx, std::vector<Triplet>& GradFx, const VecXf& dqdot, Clock& c) {
  Profiler::start("Bend Eval");
  DiffScalarBase::setVariableCount(NUM_VARS);
  float h = c.timestep();
  
  //VecXf test1 = VecXf::Zero(Fx.rows());
  //VecXf test2 = VecXf::Zero(Fx.rows());
  
  std::vector<Triplet> test1;
  std::vector<Triplet> test2;
  
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
    
    //test1.block<9,1>(3*(i-1), 0) += grad;
    
    assert(et == Implicit && "Unsupported EvalType");
    for (int j=0; j<NUM_VARS; j++) {
      // TODO: mass matrix may not be I
      Fx(3*(i-1)+j) += h*grad(j);
      
      for (int k=0; k<NUM_VARS; k++) {
        float val = h*h*hess(j,k);
        CHECK_NAN(val);
        if (val != 0) {
          GradFx.push_back(Triplet(3*(i-1)+j, 3*(i-1)+k, val));
//          test1.push_back(Triplet(3*(i-1)+j, 3*(i-1)+k, hess(j,k)));
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
    assert(chi > 0 && "Segments are pointing in exactly opposite directions!");
    Vec3f tTilde = (tPrev + tNext)/chi;
    Vec3f d1 = prevSeg.m1() + nextSeg.m1();
    Vec3f d1tilde = d1/chi;
    Vec3f d2 = prevSeg.m2() + nextSeg.m2();
    Vec3f d2tilde = d2/chi;
    Vec3f curveBinorm = (2*tPrev.cross(tNext))/chi; // Verified
    Vec2f matCurve = 0.5*Vec2f(d2.dot(curveBinorm), -d1.dot(curveBinorm)); // Verified
    
    Vec3f gradK1ePrev = (-matCurve.x()*tTilde + tNext.cross(d2tilde)) / dPrevSeg.norm(); // Verified
    Vec3f gradK1eNext = (-matCurve.x()*tTilde - tPrev.cross(d2tilde)) / dNextSeg.norm(); // Verified
    Vec3f gradK2ePrev = (-matCurve.y()*tTilde - tNext.cross(d1tilde)) / dPrevSeg.norm(); // Verified
    Vec3f gradK2eNext = (-matCurve.y()*tTilde + tPrev.cross(d1tilde)) / dNextSeg.norm(); // Verified
    
    
    // WARNING: assumes that the bending matrix is the identity.
    
    Vec2f& dRestCurve = restCurve[i-1];
    // b11*2*(k1-restk1) + (b21+b12)(k2-restk2)
    float k1coeff = 2*(matCurve.x()-dRestCurve.x()); // Verified
    // b22*2*(k2-restk2) + (b21+b12)(k1-restk1)
    float k2coeff = 2*(matCurve.y()-dRestCurve.y()); // Verified
    float totalcoeff = 0.5/voronoiCell[i-1];
    
    Vec3f gradePrev = totalcoeff * (gradK1ePrev * k1coeff + gradK2ePrev * k2coeff); // Verified
    Vec3f gradeNext = totalcoeff * (gradK1eNext * k1coeff + gradK2eNext * k2coeff); // Verified
    
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
    
    Mat3f hessK2ePrev2 = 2*matCurve.y()*tTilde2+tNextxd1TildextTilde+tNextxd1TildextTilde.transpose();
    hessK2ePrev2 -= (matCurve.y()/chi)*(Mat3f::Identity() - (tPrev*tPrev.transpose()));
    hessK2ePrev2 += 0.25*(curveBinorm*prevSeg.m1().transpose()+prevSeg.m1()*curveBinorm.transpose());
    hessK2ePrev2 /= dPrevSeg.dot(dPrevSeg);
    
    Mat3f hessK1eNext2 = 2*matCurve.x()*tTilde2+tPrevxd2TildextTilde+tPrevxd2TildextTilde.transpose();
    hessK1eNext2 -= (matCurve.x()/chi)*(Mat3f::Identity() - (tNext*tNext.transpose()));
    hessK1eNext2 += 0.25*(curveBinorm*nextSeg.m2().transpose()+nextSeg.m2()*curveBinorm.transpose());
    hessK1eNext2 /= dNextSeg.dot(dNextSeg);
    
    Mat3f hessK2eNext2 = 2*matCurve.y()*tTilde2-tPrevxd1TildextTilde-tPrevxd1TildextTilde.transpose();
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
    TDS tchi = ttprev.dot(ttnext) + 1;
    TVec3 tkb = TDS(2)*(ttprev.cross(ttnext)) / tchi;
    TDS tk1 = TDS(0.5)*TVec3(TDS(d2.x()), TDS(d2.y()), TDS(d2.z())).dot(tkb);
    TDS tk2 = TDS(-.5)*TVec3(TDS(d1.x()), TDS(d1.y()), TDS(d1.z())).dot(tkb);
    //    Eigen::Matrix<float, 3, 2> t;
    //    t << tk2.getGradient(), gradK2eNext;
    //    std::cout << t << "\n\n";
    Mat3f tdiff = tk1.getHessian() - hessK1ePrev2;
    if (tdiff.cwiseAbs().maxCoeff() >= 1)
      std::cout << tdiff << "\n\n";
    // std::cout << tdiff.norm() << " / " << tk1.getHessian().cwiseAbs().maxCoeff() << "\n";
    */
    
    
    /// Populate the 3x3 block matrix:
    /// 0 1 2
    /// 3 4 5
    /// 6 7 8
    
    Mat3f block[9];
    
    block[0] = k1coeff*hessK1ePrev2 + k2coeff*hessK2ePrev2 + 2*(gradK1ePrev*gradK1ePrev.transpose() + gradK2ePrev*gradK2ePrev.transpose());
    
    block[1] = k1coeff*(hessK1ePreveNext-hessK1ePrev2) + k2coeff*(hessK2ePreveNext-hessK2ePrev2)
    + 2*(-gradK1ePrev*(gradK1ePrev-gradK1eNext).transpose() + -gradK2ePrev*(gradK2ePrev-gradK2eNext).transpose());
    
    block[2] = -k1coeff*hessK1ePreveNext - k2coeff*hessK2ePreveNext + 2*(-gradK1ePrev*gradK1eNext.transpose() + -gradK2ePrev*gradK2eNext.transpose());
    
    block[3] = k1coeff*(hessK1ePreveNext.transpose() - hessK1ePrev2) + k2coeff*(hessK2ePreveNext.transpose() - hessK2ePrev2) + 2*((gradK1ePrev-gradK1eNext)*-gradK1ePrev.transpose() + (gradK2ePrev-gradK2eNext)*-gradK2ePrev.transpose());
    
    block[4] = k1coeff*(hessK1ePrev2 + hessK1eNext2 - hessK1ePreveNext - hessK1ePreveNext.transpose()) +
    k2coeff*(hessK2ePrev2 + hessK2eNext2 - hessK2ePreveNext - hessK2ePreveNext.transpose()) +
    2*((gradK1ePrev-gradK1eNext)*(gradK1ePrev-gradK1eNext).transpose() + (gradK2ePrev-gradK2eNext)*(gradK2ePrev-gradK2eNext).transpose());
    
    block[5] = k1coeff*(hessK1ePreveNext - hessK1eNext2) + k2coeff*(hessK2ePreveNext - hessK2eNext2)
    + 2*((gradK1ePrev-gradK1eNext)*gradK1eNext.transpose() + (gradK2ePrev-gradK2eNext)*gradK2eNext.transpose());
    
    block[6] = -k1coeff*hessK1ePreveNext.transpose() - k2coeff*hessK2ePreveNext.transpose()
    + 2*(gradK1eNext*-gradK1ePrev.transpose() + gradK2eNext*-gradK2ePrev.transpose());
    
    block[7] = k1coeff*(hessK1ePreveNext.transpose()-hessK1eNext2) + k2coeff*(hessK2ePreveNext.transpose()-hessK2eNext2) + 2*(gradK1eNext*(gradK1ePrev-gradK1eNext).transpose() + gradK2eNext*(gradK2ePrev-gradK2eNext).transpose());
    
    block[8] = k1coeff*hessK1eNext2 + k2coeff*hessK2eNext2 + 2*(gradK1eNext*gradK1eNext.transpose() + gradK2eNext*gradK2eNext.transpose());
    
    
    for(int j=0; j<3; j++) {
      for(int k=0; k<3; k++) {
        
        for(int l=0; l<9; l++) {
          if (block[l](j, k)*totalcoeff*h*h >= 1 && c.canDecreaseTimestep()) {
            Profiler::stop("Bend Eval");
            c.suggestTimestep(h/2);
            std::cerr << "Warning: Bending energy unstable. New timestep: " << c.timestep() << "\n";
            return false;
          }
        }
        
        pushBackIfNotZero(GradFx, Triplet(3*(i-1)+j, 3*(i-1)+k, h*h*totalcoeff*block[0](j, k)));
        pushBackIfNotZero(GradFx, Triplet(3*(i-1)+j, 3*i+k,     h*h*totalcoeff*block[1](j, k)));
        pushBackIfNotZero(GradFx, Triplet(3*(i-1)+j, 3*(i+1)+k, h*h*totalcoeff*block[2](j, k)));
        
        pushBackIfNotZero(GradFx, Triplet(3*i+j,     3*(i-1)+k, h*h*totalcoeff*block[3](j, k)));
        pushBackIfNotZero(GradFx, Triplet(3*i+j,     3*i+k,     h*h*totalcoeff*block[4](j, k)));
        pushBackIfNotZero(GradFx, Triplet(3*i+j,     3*(i+1)+k, h*h*totalcoeff*block[5](j, k)));
        
        pushBackIfNotZero(GradFx, Triplet(3*(i+1)+j, 3*(i-1)+k, h*h*totalcoeff*block[6](j, k)));
        pushBackIfNotZero(GradFx, Triplet(3*(i+1)+j, 3*i+k,     h*h*totalcoeff*block[7](j, k)));
        pushBackIfNotZero(GradFx, Triplet(3*(i+1)+j, 3*(i+1)+k, h*h*totalcoeff*block[8](j, k)));
        
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
    
    
    
#endif //ifdef ENABLE_AUTODIFF
    
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
#undef NUM_VARS


// STRETCHING
Stretching::Stretching(const Yarn& y, EvalType et) : YarnEnergy(y, et) { }

#define NUM_VARS 6
bool Stretching::eval(VecXf& Fx, std::vector<Triplet>& GradFx, const VecXf& dqdot, Clock& c) {
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
    if (axialStrain.getValue() >= 0.5 && c.canDecreaseTimestep()) { // Too much stretching happening here...
      Profiler::stop("Stretch Eval");
      c.suggestTimestep(h/2);
      std::cerr << "Warning: Stretching energy unstable. New timestep: " << c.timestep() << "\n";
      return false;
    }
    
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
  
  return true;
}
#undef NUM_VARS


// TWISTING

Twisting::Twisting(const Yarn& y, EvalType et) : YarnEnergy(y, et) {
  for (int i=1; i<y.numCPs()-1; i++) {
    const Segment& ePrev = y.rest().segments[i-1];
    const Segment& eNext = y.rest().segments[i];
    voronoiCell.push_back(0.5*(ePrev.length()+eNext.length()));
  }
}

bool Twisting::eval(VecXf& Fx, std::vector<Triplet>& GradFx, const VecXf& dqdot, Clock& c) {
  assert(et == Explicit && "EvalType unsupported");
  for (int i=1; i<y.numCPs()-1; i++) {
    const Segment& segPrev = y.cur().segments[i-1];
    const Segment& segNext = y.cur().segments[i];
    
    Vec3f tPrev = segPrev.vec().normalized();
    Vec3f tNext = segNext.vec().normalized();
    
    Vec3f curveBinorm = (2*tPrev.cross(tNext))/(1+tPrev.dot(tNext));
    float dThetaHat = twistMod * (segNext.getRefTwist() - (segNext.getRot() - segPrev.getRot())) / voronoiCell[i-1];
    Vec3f dxi = curveBinorm / y.rest().segments[i].length() - curveBinorm / y.rest().segments[i-1].length();
    Fx.block<3,1>(3*i, 0) += c.timestep() * dThetaHat * dxi;
  }
  return true;
}

// INTERNAL CONTACT

IntContact::IntContact(const Yarn& y, EvalType et) : YarnEnergy(y, et) {}

bool IntContact::eval(VecXf& Fx, std::vector<Triplet>& GradFx, const VecXf& dqdot, Clock& c) {
  for (int i=0; i<y.numSegs(); i++) {
    for (int j=i+2; j<y.numSegs(); j++) {
      if (i == 0 || j == 0 || i == y.numSegs()-1 || j == y.numSegs()-1) continue; // FIXME: evaluate ends!!!
      if (i-j <= 3 && j-i <= 3) continue; // Don't evaluate close edges. FIXME: should be 1 ideally

      float h = c.timestep();
      float r = constants::radius;
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
      
      Vec3f e1mid = (e1p1 + e1p2) / 2;
      Vec3f e2mid = (e2p1 + e2p2) / 2;
      
      if ((e1mid - e2mid).norm() > fmaxf(e1.length(), e2.length())) {
        if (ptd) {
          std::pair<int, int> id = ptd->id();
          if (id.first == i && id.second == j) {
            delete ptd;
            ptd = 0;
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
      int nb = 24;
      
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
      if (et == Implicit) {
        for (int k=0; k<8; k++) {
          for (int l=0; l<8; l++) {
            hess[k][l] = Mat3f::Zero();
          }
        }
      }
      
      for (int n=0; n<nb; n++) {
        for (int m=0; m<nb; m++) {
          float t1 = ((float) n) / nb;
          float t2 = ((float) m) / nb;
          Vec3f p1 = s1.eval(t1, false);
          Vec3f p2 = s2.eval(t2, false);
          Vec3f v = p2 - p1; // FIXME: shouldn't this be p1 - p2??
          float norm = v.norm();
          Vec4f u1(t1*t1*t1, t1*t1, t1, 1);
          Vec4f u2(t2*t2*t2, t2*t2, t2, 1);
          
          float dfnorm = df(norm / 2 / r) / 2 / r;
          Vec3f gradBase = v/norm;
          
          float d2fnorm = d2f(norm / 2 / r) / 4 / r / r;
          Mat3f grad2Base = v*v.transpose()/norm/norm;
          Mat3f hessBase = (Mat3f::Identity() - grad2Base) * dfnorm / norm + grad2Base * d2fnorm;
          
          for (int k=0; k<4; k++) {
            float c = u1.dot(constants::basis[k]);
            float d = -u2.dot(constants::basis[k]);
            
            Fx.block<3,1>(3*(i-1+k), 0) += h*coeff*c*gradBase;
            Fx.block<3,1>(3*(j-1+k), 0) += h*coeff*d*gradBase;
          }
          
          if (et == Implicit) {
            for (int k=0; k<4; k++) {
              for (int l=0; l<4; l++) {
                float ck = u1.dot(constants::basis[k]);
                float cl = u1.dot(constants::basis[l]);
                float dk = -u2.dot(constants::basis[k]);
                float dl = -u2.dot(constants::basis[l]);
                
                hess[k][l]     += ck * cl * hessBase;
                hess[k+4][l]   += ck * dl * hessBase;
                hess[k][l+4]   += dk * cl * hessBase;
                hess[k+4][l+4] += dk * dl * hessBase;
              }
            }
          }
        }
      }
      
      if (et == Implicit) {
        for (int k=0; k<8; k++) {
          for (int l=0; l<8; l++) {
            hess[k][l] *= h*h*coeff;
          }
        }

        for (int k=0; k<4; k++) {
          for (int l=0; l<4; l++) {
            for (int p=0; p<3; p++) {
              for (int q=0; q<3; q++) {
                pushBackIfNotZero(GradFx, Triplet(3*(i-1+k)+p, 3*(i-1+l)+q, hess[k][l](p, q)));
                pushBackIfNotZero(GradFx, Triplet(3*(i-1+k)+p, 3*(j-1+l)+q, hess[k+4][l](p, q)));
                pushBackIfNotZero(GradFx, Triplet(3*(j-1+k)+p, 3*(i-1+l)+q, hess[k][l+4](p, q)));
                pushBackIfNotZero(GradFx, Triplet(3*(j-1+k)+p, 3*(j-1+l)+q, hess[k+4][l+4](p, q)));
              }
            }
          }
        }
      }

      
    }
  }
  return true;
}
