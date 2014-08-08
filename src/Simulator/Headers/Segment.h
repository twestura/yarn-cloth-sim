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
#include "Util.h"
#include "CtrlPoint.h"

class Segment
{
  /// Reference to the first control point that defines this segment.
  const CtrlPoint& first;
  /// Reference to the second control point that defines this segment.
  const CtrlPoint& second;
  /// The u vector of the twist-free reference (Bishop) frame.
  Vec3e u;
  /// Rotation in radians of the material frame.
  real rot = 0.0;
  /// Twist in radians from reference frame parallel transport through time.
  real refTwist = 0.0;
  
  bool frozen = false;
  Vec3e frozVec, frozV;
  
  const Vec3e inline compVec() const {
    return second.pos - first.pos;
  }
  
  const Vec3e inline compV() const {
    return vec().cross(u).normalized();
  }

public:
  /// Default Segment constructor
  Segment(const CtrlPoint& a, const CtrlPoint& b, Vec3e f) : first(a), second(b), u(f) {}
  /// Get the first Control Point of this Segment.
  const CtrlPoint& getFirst() const { return first; }
  /// Get the second Control Point of this Segment.
  const CtrlPoint& getSecond() const { return second; }
  
  /// Calculate the vector that represents this segment. Not necessarily unit length.
  const Vec3e inline vec() const {
    if (!frozen) { return compVec(); }
    return frozVec;
  }
  /// Calculate the length of the segment.
  const real inline length() const { return vec().norm(); }
  
  // TODO: rename these?
  const Vec3e inline getU() const { return u; }
  void inline setU(const Vec3e newU) { u = newU; }
  
  /// Calculate the reference (Bishop) frame vector v.
  const Vec3e inline v() const {
    if (!frozen) { return compV(); }
    return frozV;
  }
  /// Set the material frame rotation.
  void inline setRot(const real f) { rot = f; }
  /// Get the material frame rotation.
  const real inline getRot() const { return rot; }
  /// Calculate the material frame vector m1.
  const Vec3e inline m1() const {
    return cos(rot) * u + sin(rot) * v();
  }
  /// Calculate the material frame vector m2.
  const Vec3e inline m2() const {
    return -sin(rot) * u + cos(rot) * v();
  }
  /// Get twist in reference frame from previous frame.
  const real inline getRefTwist() const { return refTwist; }
  
  /// Returns a parallel transported vector given a previous vector and its orthogonal u component.
  Vec3e static parallelTransport(const Vec3e vecPrev, const Vec3e vecCur, const Vec3e uPrev) {
    Vec3e cross = vecPrev.cross(vecCur).normalized();
    real cosT = vecCur.dot(vecPrev)/(vecCur.norm() * vecPrev.norm());
    if (!cross.hasNaN() && cosT < 1.0 && cosT >= -1.0) {
      // Form rotation matrix
      real oneMinusCosT = 1.0 - cosT;
      real sinT = sqrt(1.0 - (cosT * cosT));
      real xyc = cross.x() * cross.y() * oneMinusCosT;
      real xzc = cross.x() * cross.z() * oneMinusCosT;
      real yzc = cross.y() * cross.z() * oneMinusCosT;
      real xs = cross.x() * sinT;
      real ys = cross.y() * sinT;
      real zs = cross.z() * sinT;
      Mat3e rotMat;
      rotMat << cosT + cross.x() * cross.x() * oneMinusCosT, xyc - zs, xzc + ys,
                xyc + zs, cosT + cross.y() * cross.y() * oneMinusCosT, yzc - xs,
                xzc - ys, yzc + xs, cosT + cross.z() * cross.z() * oneMinusCosT;
      return rotMat * uPrev;
    }
    return uPrev;
  }
  
  /// Set the reference frame by parallel transporting via the given vector. No reference twist is
  /// assumed to occur.
  void parallelTransport(const Segment& prevSeg) {
    u = parallelTransport(prevSeg.vec(), vec(), prevSeg.u);
    CHECK_NAN_VEC(u);
  }
  
  /// Parallel transports a segment through time, recording the amount of reference twist
  /// accumulated through space.
  void parallelTransport(const Segment& prevTimeSeg, const Segment& prevSpaceSeg) {
    // Parallel transport through time to update this.u
    parallelTransport(prevTimeSeg);
    // Find the amount of twist in the reference frame. First, parallel transport
    // prevSpaceSeg through space so the u vectors lie in the plane with the normal
    // parallel to this.vec.
    Vec3e uRef = parallelTransport(prevSpaceSeg.vec(), vec(), prevSpaceSeg.u);
    CHECK_NAN_VEC(uRef);
    real cosTwist = u.normalized().dot(uRef.normalized());
    // Now find the angle between the reference (space-parallel transported u) and this.u
    if (cosTwist >= 1.0) { // Avoid values like 1.0000000012 that introduce NaNs
      refTwist = 0.0;
    } else if (cosTwist <= -1.0) {
      refTwist = constants::pi;
    } else {
      refTwist = acos(cosTwist);
    }
    // Flip the sign if necessary
    if (v().normalized().dot(uRef) > 0.0) {
      refTwist = -refTwist;
    }
    CHECK_NAN(refTwist);
  }
  
  void inline setFrozen(bool f) {
    frozen = f;
    if (frozen) {
      frozVec = compVec();
      frozV = compV();
    }
  }
  
};


#endif /* defined(__Visualizer__Segment__) */
