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
#include "Constants.h"
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
  /// The number of full (2*pi) twists in the reference frame.
  int numTwists = 0;
  
public:
  /// Default Segment constructor
  Segment(const CtrlPoint& a, const CtrlPoint& b, Vec3f f) : first(a), second(b), u(f) {}
  /// Get the first Control Point of this Segment.
  const CtrlPoint& getFirst() const { return first; }
  /// Get the second Control Point of this Segment.
  const CtrlPoint& getSecond() const { return second; }
  
  /// Calculate the vector that represents this segment. Not necessarily unit length.
  const Vec3f inline vec() const { return second.pos - first.pos; }
  /// Calculate the length of the segment.
  const float inline length() const { return vec().norm(); }
  
  // TODO: remove these?
  const Vec3f inline getU() const { return u; }
  void inline setU(const Vec3f newU) { u = newU; }
  
  /// Calculate the reference (Bishop) frame vector v.
  const Vec3f inline v() const { return vec().cross(u).normalized(); }
  /// Set the material frame rotation.
  void inline setRot(const float f) { rot = f; }
  /// Get the material frame rotation.
  const float inline getRot() const { return rot; }
  /// Calculate the material frave vector m1.
  const Vec3f inline m1() const {
    Eigen::Quaternionf q(Eigen::AngleAxisf(getRot()-getRefTwist(), vec().normalized()));
    return q * u;
  }
  /// Calculate the material frave vector m2.
  const Vec3f inline m2() const {
    Eigen::Quaternionf q(Eigen::AngleAxisf(getRot()-getRefTwist(), vec().normalized()));
    return q * v();
  }
  /// Get twist in reference frame from previous frame.
  const float getRefTwist() const { return refTwist + 2*constants::pi*numTwists; }
  
  /// Set the reference frame by parallel transporting via the given vector. No reference twist is
  /// assumed to occur.
  void parallelTransport(Segment& prevSeg) {
    Vec3f vprev = prevSeg.vec();
    Vec3f vcur  = vec();
    Vec3f cross = vprev.cross(vcur);
    cross.normalize();
    float twist = acos(vcur.dot(vprev)/(vcur.norm()*vprev.norm()));
    u = prevSeg.u;
    // If the angle is too small, cross becomes inaccurate. In order to prevent error propagation,
    // it's better to pretend the angle is 0. Also guards against NaNs from acos.
    if (twist > .00001) {
      Eigen::Quaternionf q(Eigen::AngleAxisf(twist, cross));
      u = q * u;
    }
  }
  
  /// Set the reference frame by parallel transporting via prevSeg, then update refSeg with the
  /// amount of twist accumulated by the transport.
  /// WARNING: does not account for twists greater than pi correctly.
  void parallelTransport(Segment& prevSeg, Segment& refSeg) {
    parallelTransport(prevSeg);
    Vec3f ucur = getU().normalized();
    Vec3f vcur = v().normalized();
    Vec3f uref = refSeg.getU().normalized();
    float val = ucur.dot(uref);
    float sign = vcur.dot(uref);
    float oldTwist = refSeg.refTwist;
    if (val >= 1) { // Avoid values like 1.0000000012 that introduce NaNs
      refSeg.refTwist = 0;
    } else if (val <= -1) {
      refSeg.refTwist = constants::pi;
    } else {
      refSeg.refTwist = acos(val);
    }
    if (sign < 0) { // check if we need to switch sign
      refSeg.refTwist *= -1;
    }
    
    // Account for twists >|pi|. Assumes that twists are not greater than pi between each transport.
    float diff = refSeg.refTwist - oldTwist;
    if (diff < -constants::pi) {
      refSeg.numTwists -= 1;
    } else if (diff > constants::pi) {
      refSeg.numTwists += 1;
    }
  }
};


#endif /* defined(__Visualizer__Segment__) */
