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
  /// The number of full (2*pi) twists in the reference frame.
  int numTwists = 0;
  /// The change to numTwists this frame. Use this to keep numTwists in sync between yarns.
  int deltaTwists = 0;
  
public:
  /// Default Segment constructor
  Segment(const CtrlPoint& a, const CtrlPoint& b, Vec3e f) : first(a), second(b), u(f) {}
  /// Get the first Control Point of this Segment.
  const CtrlPoint& getFirst() const { return first; }
  /// Get the second Control Point of this Segment.
  const CtrlPoint& getSecond() const { return second; }
  
  /// Calculate the vector that represents this segment. Not necessarily unit length.
  const Vec3e inline vec() const { return second.pos - first.pos; }
  /// Calculate the length of the segment.
  const real inline length() const { return vec().norm(); }
  
  // TODO: rename these?
  const Vec3e inline getU() const { return u; }
  void inline setU(const Vec3e newU) { u = newU; }
  
  /// Calculate the reference (Bishop) frame vector v.
  const Vec3e inline v() const { return vec().cross(u).normalized(); }
  /// Set the material frame rotation.
  void inline setRot(const real f) { rot = f; }
  /// Get the material frame rotation.
  const real inline getRot() const { return rot; }
  /// Calculate the material frame vector m1.
  const Vec3e inline m1() const {
    // return cos(getRot()) * u + sin(getRot()) * v();
    // .. or, equivalently:
    Eigen::Quaternion<real> q(Eigen::AngleAxis<real>(getRot(), vec().normalized()));
    return q * u;
  }
  /// Calculate the material frame vector m2.
  const Vec3e inline m2() const {
    Eigen::Quaternion<real> q(Eigen::AngleAxis<real>(getRot(), vec().normalized()));
    return q * v();
  }
  /// Get twist in reference frame from previous frame.
  const real inline getRefTwist() const { return refTwist; } // + 2.0*constants::pi*numTwists; }
  
  /// Update numTwists this tick.
  void inline updateTwists(const Segment& oldSeg) {
    numTwists += oldSeg.deltaTwists;
  }
  
  /// Returns a parallel transported vector given a previous vector and its orthogonal u component.
  Vec3e static parallelTransport(const Vec3e vecPrev, const Vec3e vecCur, const Vec3e uPrev) {
    Vec3e cross = vecPrev.cross(vecCur).normalized();
    real twist = acos(vecCur.dot(vecPrev)/(vecCur.norm() * vecPrev.norm()));
    if (cross.allFinite() && twist > 1.0e-7) {
      Eigen::Quaternion<real> q(Eigen::AngleAxis<real>(twist, cross));
      return q * uPrev;
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
    real oldTwist = refTwist;
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
    
    // Account for twists >|pi|. Assumes that twists are not greater than pi between each transport.
    real diff = refTwist - oldTwist;
    if (diff < -constants::pi) {
      deltaTwists = 1;
      numTwists += 1;
    } else if (diff > constants::pi) {
      deltaTwists = -1;
      numTwists -= 1;
    } else {
      deltaTwists = 0;
    }
    
  }
};


#endif /* defined(__Visualizer__Segment__) */
