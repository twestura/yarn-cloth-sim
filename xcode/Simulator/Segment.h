//
//  Segment.h
//  Visualizer
//
//  Created by eschweickart on 2/12/14.
//
//

#ifndef __Visualizer__Segment__
#define __Visualizer__Segment__

#include "Eigen/Geometry"
#include "CtrlPoint.h"

class Segment
{
  /// Reference to the first control point that defines this segment.
  const CtrlPoint& first;
  /// Reference to the second control point that defines this segment.
  const CtrlPoint& second;
  /// The u vector of the twist-free reference (Bishop) frame.
  Vec3f u;
  /// Rotation in radians of the material frame.
  float rot = 0;
  /// Twist in radians from reference frame parallel transport through time.
  float refTwist = 0;
  
public:
  /// Default Segment constructor
  Segment(const CtrlPoint& a, const CtrlPoint& b) : first(a), second(b) {}
  /// Calculate the vector that represents this segment. Not necessarily unit length.
  const Vec3f inline vec() const { return second.pos - first.pos; }
  /// Calculate the length of the segment.
  const float inline length() const { return vec().norm(); }
  /// Calculate the reference (Bishop) frame vector v.
  const Vec3f inline v() const { return vec().cross(u).normalized(); }
  /// Calculate the material frave vector m1.
  const Vec3f m1() const {
    Eigen::Quaternionf q = Eigen::Quaternionf(Eigen::AngleAxisf(rot, vec().normalized()));
    return q * u;
  }
  /// Calculate the material frave vector m2.
  const Vec3f m2() const {
    Eigen::Quaternionf q = Eigen::Quaternionf(Eigen::AngleAxisf(rot, vec().normalized()));
    return q * v();
  }
  /// Get twist in reference frame from previous frame.
  const float getRefTwist() {
    return refTwist;
  }
  /// Set the reference frame by parallel transporting through time.
  void parallelTransport(Segment& prevSeg) {
    Vec3f vprev = prevSeg.vec();
    Vec3f vcur  = vec();
    Vec3f cross = vprev.cross(vcur);
    cross.normalize();
    refTwist = acos(vcur.dot(vprev)/(vcur.norm()*vprev.norm()));
    u = prevSeg.u;
    // If the angle is too small, cross becomes inaccurate. In order to prevent error propagation,
    // it's better to pretend the angle is 0.
    if (refTwist > .00001) {
      Eigen::Quaternionf q = Eigen::Quaternionf(Eigen::AngleAxisf(refTwist, cross));
      u = q * u;
    }
  }
};


#endif /* defined(__Visualizer__Segment__) */
