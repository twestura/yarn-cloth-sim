//
//  Yarn.h
//  Visualizer
//
//  Created by eschweickart on 2/12/14.
//
//

#ifndef __Visualizer__Yarn__
#define __Visualizer__Yarn__

#include "Util.h"
#include "Constants.h"
#include "CrossSection.h"

/// A Yarn at a specific point in time.
class YarnStr {
private:
  bool frozen = false;
  std::vector<Vec3e> frozV;
  std::vector<Vec3e> frozVec;
  
  const Vec3e inline compVec(size_t index) const {
    return pos.segment<3>(3*(index+1)) - pos.segment<3>(3*index);
  }
  
  const Vec3e inline compV(size_t index) const {
    return vec(index).cross(u[index]).normalized();
  }
  
public:
  VecXe pos;
  VecXe vel;
  VecXe acc;
  
  VecXe rot;
  VecXe refTwist;

  std::vector<Vec3e> u;
  
  /// Calculate a vector representing a given segment.
  const Vec3e inline vec(size_t index) const {
    if (!frozen) { return compVec(index); }
    return frozVec[index];
  }
  
  /// Calculate the reference (Bishop) frame vector v at a given segment.
  const Vec3e inline v(size_t index) const {
    if (!frozen) { return compV(index); }
    return frozV[index];
  }
  
  /// Calculate the length of a given segment.
  const real inline length(size_t index) const { return vec(index).norm(); }
  
  /// Calculate the material frame vector m1 at a given segment.
  const Vec3e inline m1(size_t index) const {
    return cos(rot(index)) * u[index] + sin(rot(index)) * v(index);
  }
  /// Calculate the material frame vector m2 at a given segment.
  const Vec3e inline m2(size_t index) const {
    return -sin(rot(index)) * u[index] + cos(rot(index)) * v(index);
  }
  
#define POS(i) pos.segment<3>(3*(i))
#define VEL(i) vel.segment<3>(3*(i))
#define ACC(i) acc.segment<3>(3*(i))
  
  /// Assures that the rod will not be modified any time soon. This allows certain commonly-used
  /// values to be pre-calculated and cached. Unfreeze the rod before modifying it!
  void inline setFrozen(bool f) {
    frozen = f;
    if (frozen) {
      if (frozVec.empty() || frozV.empty()) {
        for (int i=0; i<u.size(); i++) {
          frozVec.push_back(compVec(i));
          frozV.push_back(compV(i));
        }
      } else {
        for (int i=0; i<u.size(); i++) {
          frozVec[i] = compVec(i);
          frozV[i] = compV(i);
        }
      }
    }
  }
  
  /// Parallel transports all segments through time, recording the amount of reference twist
  /// accumulated through space. Pass the YarnStr in time as a parameter.
  void updateReferenceFrames(const YarnStr& prevTime) {
    for (int i=0; i<u.size(); i++) {
      // Parallel transport through time
      u[i] = parallelTransport(prevTime.vec(i), vec(i), prevTime.u[i]).normalized();
      
      if (i == 0) continue;
      
      // Parallel transport through space
      Vec3e uRef = parallelTransport(vec(i-1), vec(i), u[i-1]).normalized();
      
      // Record the twist between the two
      real cosTwist = u[i].dot(uRef);
      if (cosTwist >= 1.0) {
        refTwist(i) = 0.0;
      } else if (cosTwist <= -1.0) {
        refTwist(i) = constants::pi;
      } else {
        refTwist(i) = acos(cosTwist);
      }
      // Flip the sign if necessary
      if (v(i).dot(uRef) > 0.0) {
        refTwist(i) = -refTwist(i);
      }
    }
  }
};

struct Mass {
public:
  VecXe diag;
  Eigen::SparseMatrix<real> sparse;
  real total;
};

/// A Yarn stepped though points in time.
class Yarn {
private:
  size_t nCPs;
  size_t ndof;
  size_t nSegs;
  
  /// The mass matrix.
  Mass mass;
  /// The inverse of the mass matrix, cached here to save on computation.
  Mass invMass;
  
  /// Information related to the cross-section of the rod at a particular point.
  CrossSection* cs;
  
  /// The yarn at rest.
  YarnStr  restYS;
  /// The yarn at the current point in time.
  YarnStr* curYS;
  /// The yarn at the next point in time.
  YarnStr* nextYS;
  /// The rest Voronoi lengths defined at each internal control point of the yarn.
  std::vector<real> rvl;
  /// The rest curvature with respect to the previous edge (﻿﻿omega bar ^i _i)
  std::vector<Vec2e> rcp;
  /// The rest curvature with respect to the next edge (omega bar ^i _i+1)
  std::vector<Vec2e> rcn;
  /// The rest curvature defined at each internal control point of the yarn.
  std::vector<Vec2e> rc;
  
  
public:
  real shearModulus;
  real youngsModulus;
  
  /// Constructs a yarn that is initially at rest given a vector of initial positions.
  /// The U vector for the first segment is then propagated down the yarn by parallel transport
  /// through space. A vector of lumped masses may be optionally passed in; if it is null, all
  /// masses are set uniformly to 1kg. Young's and shear moduli may also be specified.
  Yarn(VecXe pos, Vec3e u0, VecXe* masses = nullptr,
       real youngs = constants::youngsModulus, real shear = constants::shearModulus) :
  youngsModulus(youngs), shearModulus(shear) {
    ndof = pos.rows();
    nCPs = ndof/3;
    assert(nCPs*3 == ndof && "Malformed position vector");
    nSegs = nCPs-1;
    
    restYS.pos = pos;
    restYS.vel = VecXe::Zero(ndof);
    restYS.acc = VecXe::Zero(ndof);
    
    restYS.rot = VecXe::Zero(nSegs);
    restYS.refTwist = VecXe::Zero(nSegs);
    
    restYS.u.push_back(u0);
    
    for(int i=1; i<nSegs; i++) {
      restYS.u.push_back(parallelTransport(restYS.vec(i-1), restYS.vec(i), restYS.u[i-1]));
    }
    
    curYS  = new YarnStr(restYS);
    nextYS = new YarnStr(restYS);
    
    restYS.setFrozen(true);
    
    for (int i=1; i<numCPs()-1; i++) {
      Vec3e ePrev = restYS.vec(i-1);
      Vec3e eNext = restYS.vec(i);
      
      Vec3e curveBinorm = 2.0*ePrev.cross(eNext) / (ePrev.norm()*eNext.norm() + ePrev.dot(eNext));
      
      CHECK_NAN_VEC(curveBinorm);
      
      Vec2e restMatCurvePrev(curveBinorm.dot(restYS.m2(i-1)), -(curveBinorm.dot(restYS.m1(i-1))));
      Vec2e restMatCurveNext(curveBinorm.dot(restYS.m2(i)), -(curveBinorm.dot(restYS.m1(i))));
      Vec2e restMatCurve = 0.5*(restMatCurvePrev + restMatCurveNext);
      
      rcp.push_back(restMatCurvePrev);
      rcn.push_back(restMatCurveNext);
      rc.push_back(restMatCurve);
      rvl.push_back(0.5*(ePrev.norm()+eNext.norm()));
    }
    
    // Set mass matrix
    if (masses && masses->rows() == nCPs) {
      mass.diag = *masses;
      mass.total = mass.diag.sum();
      std::vector<Triplet> triplets;
      triplets.reserve(ndof);
      for (int i=0; i<ndof; i++) {
        triplets.push_back(Triplet(i, i, (*masses)(i/3)));
      }
      mass.sparse.resize(ndof, ndof);
      mass.sparse.setFromTriplets(triplets.begin(), triplets.end());
      
      invMass.diag = mass.diag.cwiseInverse();
      invMass.sparse = mass.sparse.cwiseInverse();
    } else {
      std::cout << "Warning: Mass matrix set to identity.\n";
      mass.diag = VecXe::Ones(nCPs);
      mass.total = nCPs;
      mass.sparse.resize(ndof, ndof);
      mass.sparse.setIdentity();
      
      invMass.diag = VecXe::Ones(nCPs);
      invMass.sparse.resize(ndof, ndof);
      invMass.sparse.setIdentity();
    }
    
    cs = new CrossSection(new Ellipse(constants::radius, constants::radius,
                                       youngsModulus, shearModulus));
  }
  
  // WARNING: The following methods are safe as long as curYS and nextYS are always allocated upon
  // construction and are never deleted until destruction.
  /// Get a const reference to the yarn at the current point in time.
  const YarnStr& cur() const  { return *curYS;  }
  /// Get a reference to the yarn at the next point in time.
  YarnStr& next() { return *nextYS; }
  /// Get a const reference to the yarn at the next point in time.
  const YarnStr& next() const { return *nextYS; }
  /// Get a const reference to the yarn at rest.
  const YarnStr& rest() const { return restYS; }
  
  /// Get the number of control points on the yarn.
  const inline size_t numCPs() const { return nCPs; }
  /// Get the number of control points associated with 2 edges.
  const inline size_t numIntCPs() const { return nCPs > 2 ? nCPs-2 : 0; }
  /// Get the number of segments in the yarn.
  const inline size_t numSegs() const { return nSegs; }
  /// Get the number of degrees of freedom of the yarn.
  const inline size_t numDOF() const { return ndof; }
  
  /// Swaps the current and next yarns, e.g. at the end of a timestep.
  void inline swapYarns() {
    YarnStr* temp = curYS;
    curYS = nextYS;
    nextYS = temp;
    
    curYS->setFrozen(true);
    nextYS->setFrozen(false);
  }
  
  // Get the rest Voronoi length for an internal control point.
  const inline real restVoronoiLength(size_t index) const {
    assert(index > 0 && "Voronoi length undifined at this control point.");
    return rvl[index-1];
  }
  
  // Get the rest curvature for an internal control point.
  const inline Vec2e& restCurvePrev(size_t index) const {
    assert(index > 0 && "Curvature undefined at this control point.");
    return rcp[index-1];
  }
  
  // Get the rest curvature for an internal control point.
  const inline Vec2e& restCurveNext(size_t index) const {
    assert(index > 0 && "Curvature undefined at this control point.");
    return rcn[index-1];
  }
  
  // Get the rest curvature for an internal control point.
  const inline Vec2e& restCurve(size_t index) const {
    assert(index > 0 && "Curvature undefined at this control point.");
    return rc[index-1];
  }
  
  /// Get the cross-section of the rod.
  const inline CrossSection& getCS() const { return *cs; }
  
  /// Get the yarn's mass matrix.
  const inline Mass& getMass() const { return mass; }
  /// Get the yarn's inverse mass matrix.
  const inline Mass& getInvMass() const { return invMass; }
  
  /// Get the center of mass of the current yarn.
  const Vec3e getCurCoM() const {
    Vec3e com = Vec3e::Zero();
    VecXe posTimesMass = getMass().sparse * cur().pos;
    for (int i=0; i<numCPs(); i++) {
      com += posTimesMass.segment<3>(3*i);
    }
    return com / getMass().total;
  }
  
  ~Yarn() {
    delete curYS;
    delete nextYS;
  }
};

class Spline {
  Vec3e p[4];
  
public:
  
  Spline(const Vec3e p0, const Vec3e p1, const Vec3e p2, const Vec3e p3) {
    p[0] = p0; p[1] = p1; p[2] = p2; p[3] = p3;
  }

  Vec3e eval(real t) const {
    Vec4e u(t*t*t, t*t, t, 1);
    Vec4e uBasis(u.dot(constants::basis[0]),
                 u.dot(constants::basis[1]),
                 u.dot(constants::basis[2]),
                 u.dot(constants::basis[3]));
    Vec4e x(p[0].x(), p[1].x(), p[2].x(), p[3].x());
    Vec4e y(p[0].y(), p[1].y(), p[2].y(), p[3].y());
    Vec4e z(p[0].z(), p[1].z(), p[2].z(), p[3].z());
    return Vec3e(uBasis.dot(x), uBasis.dot(y), uBasis.dot(z));
  }
  
};

#endif /* defined(__Visualizer__Yarn__) */
