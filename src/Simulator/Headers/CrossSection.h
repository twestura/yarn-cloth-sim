//
//  CrossSection.h
//  Visualizer
//
//  Created by eschweickart on 7/10/14.
//
//

#ifndef Visualizer_CrossSection_h
#define Visualizer_CrossSection_h

#include "autodiff.h"

typedef DScalar2<real, 9, Eigen::Matrix<real, 9, 1>, Eigen::Matrix<real, 9, 9>> DS92;
typedef DS92::DVector2 DV92;

class Shape {
protected:
  real  sArea;
  Mat2e sAreaMoment;
  
  Mat2e sBendMat;
  Mat2e sBendMat2;
  real  sTwist;
  real  sStretch;
  
  real  sRadius;
  
  void setStiffnesses(real youngs, real shear) {
// FIXME: This is pretty hacky, but it gets a more realistic yarn simulation.
#ifdef SIMULATOR
    sBendMat = Mat2e::Identity();
#else
    sBendMat = sAreaMoment * youngs;
#endif // ifdef SIMULATOR
    sBendMat2 = sBendMat + sBendMat.transpose();
    // FIXME: This may be wrong, but works for elliptic x-sections
    sTwist = sAreaMoment.sum() * shear;
    sStretch = sArea * youngs;
  }
public:
  inline const real area() const { return sArea; }
  inline const Mat2e& areaMoment() const { return sAreaMoment; }
  
  inline const Mat2e& bendMat() const { return sBendMat; }
  inline const Mat2e& bendMat2() const { return sBendMat2; }
  inline Vec2e dBend(const Vec2e& v) const { return sBendMat2 * v; }
  DV92 autodBend(const DV92& v) const {
    return DV92(v.x() * sBendMat(0, 0) + v.y() * sBendMat(0, 1),
                v.x() * sBendMat(1, 0) + v.y() * sBendMat(1, 1));
  }
  inline const real twistCoeff() const { return sTwist; }
  inline const real stretchCoeff() const { return sStretch; }
  
  inline const real approxRadius() const { return sRadius; }
  
  virtual const real calcSample(const Vec3e& ear, const Vec3e& jerk) const =0;
  virtual ~Shape() { }
};

class Ellipse : public Shape {
protected:
  real r1;
  real r2;
public:
  Ellipse(real r1, real r2, real youngs, real shear) : r1(r1), r2(r2) {
    sArea = constants::pi * r1 * r2;
    sAreaMoment << constants::pi / 4.0 * r1 * r1 * r1 * r2, 0.0,
                   0.0, constants::pi / 4.0 * r1 * r2 * r2 * r2;
    
    sRadius = r1 + r2 / 2.0;
    
    setStiffnesses(youngs, shear);
  }
  
  // FIXME: May need a better approximation for long, thin ellipses.
  const real calcSample(const Vec3e& ear, const Vec3e& jerk) const {
    // Approximate as circle
    return (constants::rhoAir * sRadius * sRadius * sRadius /
            (2.0 * constants::cAir * ear.norm() * ear.norm())) * ear.dot(jerk);
  }
};

class Rectangle : public Shape {
protected:
  real s1;
  real s2;
public:
  Rectangle(real s1, real s2, real youngs, real shear) : s1(s1), s2(s2) {
    sArea = s1 * s2;
    sAreaMoment << s1 * s1 * s1 * s2 / 12.0, 0.0,
                   0.0, s1 * s2 * s2 * s2 / 12.0;
    
    sRadius = s1 + s2 / 4.0;
    
    setStiffnesses(youngs, shear);
  }
  
  // FIXME: there should be a better approximation for this
  const real calcSample(const Vec3e& ear, const Vec3e& jerk) const {
    // Approximate as circle
    return (constants::rhoAir * sRadius * sRadius * sRadius /
            (2.0 * constants::cAir * ear.norm() * ear.norm())) * ear.dot(jerk);
  }
};

class CrossSection {
protected:
  Shape* constShape = nullptr;
  std::vector<Shape*> varShape;
public:
  CrossSection(Shape* constShape) : constShape(constShape) { }
  
  const Shape& operator[](size_t i) const {
    if (constShape) { return *constShape; }
    return *(varShape[i]);
  }
  
  ~CrossSection() {
    if (constShape) {
      delete constShape;
    } else {
      for (Shape* s : varShape) {
        delete s;
      }
    }
  }
};

#endif
