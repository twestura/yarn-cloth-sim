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

/// An abstract class representing the shape of a cross section.
class Shape {
protected:
  /// The area of the cross section in m^2.
  real  sArea;
  /// The 2D area moment of the cross section in m^4, aligned with the reference frame.
  Mat2e sAreaMoment;
  
  /// The 2D bending matrix; used for computing bending forces.
  Mat2e sBendMat;
  /// The sum of the bending matrix and its transpose, cached here for convenience.
  Mat2e sBendMat2;
  /// The twisting stiffness implied by this cross section.
  real  sTwist;
  /// The stretching stiffness implied by this cross section.
  real  sStretch;
  /// The average radius of this cross section.
  real  sRadius;
  
  /// Given material parameters, set the bending matrix, twisting stiffness, and stretch stiffness.
  void setStiffnesses(real youngs, real shear) {
    sBendMat = sAreaMoment * youngs;
    sBendMat2 = sBendMat + sBendMat.transpose();
    // FIXME: This may be wrong, but works for elliptic x-sections
    sTwist = sAreaMoment.sum() * shear;
    sStretch = sArea * youngs;
  }
public:
  /// Returns the area of the cross section in m^2.
  inline const real area() const { return sArea; }
  /// Returns the 2D area moment of the cross-section in m^4, aligned with the reference frame.
  inline const Mat2e& areaMoment() const { return sAreaMoment; }
  
  /// Returns the 2D bending matrix for the cross section.
  inline const Mat2e& bendMat() const { return sBendMat; }
  /// Returns the sum of the bending matrix and its transpose.
  inline const Mat2e& bendMat2() const { return sBendMat2; }
  /// Returns the first derivative of the bending energy with respect to 2D curvature.
  inline Vec2e dBend(const Vec2e& v) const { return sBendMat2 * v; }
  /// Returns the first and second derivatives of the bending energy when using the autodiff.
  DV92 autodBend(const DV92& v) const {
    return DV92(v.x() * sBendMat(0, 0) + v.y() * sBendMat(0, 1),
                v.x() * sBendMat(1, 0) + v.y() * sBendMat(1, 1));
  }
  /// Returns the cross section's twisting stiffness.
  inline const real twistCoeff() const { return sTwist; }
  /// Returns the cross section's stretching stiffness.
  inline const real stretchCoeff() const { return sStretch; }
  
  /// Returns the approximate (average) radius of the cross section.
  inline const real approxRadius() const { return sRadius; }
  
  /// Given the position of a listener and the jerk of the rod at a particular point, generate the
  /// pressure contribution at the listener's position.
  virtual const real calcSample(const Vec3e& ear, const Vec3e& jerk) const =0;
  virtual ~Shape() { }
};

/// An elliptic cross section.
class Ellipse : public Shape {
protected:
  /// The radius of the ellipse in the u direction.
  real r1;
  /// The radius of the ellipse in the v direction.
  real r2;
public:
  Ellipse(real r1, real r2, real youngs, real shear) : r1(r1), r2(r2) {
    sArea = constants::pi * r1 * r2;
    sAreaMoment << constants::pi / 4.0 * r1 * r2 * r2 * r2, 0.0,
                   0.0, constants::pi / 4.0 * r1 * r1 * r1 * r2;
    
    sRadius = r1 + r2 / 2.0;
    
    setStiffnesses(youngs, shear);
  }
  
  // FIXME: May need a better approximation for long, thin ellipses.
  const real calcSample(const Vec3e& ear, const Vec3e& jerk) const {
    // Approximate as circle
    return (constants::rhoAir * sRadius * sRadius * sRadius /
            (2.0 * constants::cAir * ear.dot(ear))) * ear.dot(jerk);
  }
};

/// A rectangular cross section.
class Rectangle : public Shape {
protected:
  // FIXME
  /// The length of the side orthogonal to the u direction
  real s1;
  /// The length of the side orthogonal to the v direction
  real s2;
public:
  Rectangle(real s1, real s2, real youngs, real shear) : s1(s1), s2(s2) {
    sArea = s1 * s2;
    sAreaMoment << s1 * s2 * s2 * s2 / 12.0, 0.0,
                   0.0, s1 * s1 * s1 * s2 / 12.0;
    
    sRadius = s1 + s2 / 4.0;
    
    setStiffnesses(youngs, shear);
  }
  
  // FIXME: there should be a better approximation for this
  const real calcSample(const Vec3e& ear, const Vec3e& jerk) const {
    // Approximate as circle
    return (constants::rhoAir * sRadius * sRadius * sRadius /
            (2.0 * constants::cAir * ear.dot(ear))) * ear.dot(jerk);
  }
};

/// A class representing the cross section of either the entire rod or at a single point on the rod.
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
