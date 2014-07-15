//
//  Yarn.h
//  Visualizer
//
//  Created by eschweickart on 2/12/14.
//
//

#ifndef __Visualizer__Yarn__
#define __Visualizer__Yarn__

#include "Segment.h"
#include "CrossSection.h"
#include <vector>

/// A Yarn at a specific point in time.
struct YarnStr {
public:
  /// Control points (n+1 total) that define the yarn's position.
  std::vector<CtrlPoint> points;
  /// Segments (n total) that define the yarn.
  std::vector<Segment> segments;
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
  real shearModulus;
  real youngsModulus;
  
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
  /// Constructs a yarn that is initially at rest given a vector of initial positions.
  /// The U vector for the first segment is then propagated down the yarn by parallel transport
  /// through space. A vector of lumped masses may be optionally passed in; if it is null, all
  /// masses are set uniformly to 1kg. Young's and shear moduli may also be specified.
  Yarn(std::vector<Vec3e>& points, Vec3e u0, VecXe* masses = nullptr,
       real youngs = constants::youngsModulus, real shear = constants::shearModulus) :
  youngsModulus(youngs), shearModulus(shear) {
    
    for (Vec3e p : points) {
      CHECK_NAN_VEC(p);
      CtrlPoint cp;
      cp.pos = p;
      cp.vel.setZero();
      cp.accel.setZero();
      restYS.points.push_back(cp);
    }
    curYS  = new YarnStr(restYS);
    nextYS = new YarnStr(restYS);
    for (int i=0; i<numCPs()-1; i++) {
      Segment restSeg(restYS.points[i], restYS.points[i+1], u0);
      restYS.segments.push_back(restSeg);
      
      Segment curSeg(curYS->points[i], curYS->points[i+1], u0);
      curYS->segments.push_back(curSeg);
      
      Segment nextSeg(nextYS->points[i], nextYS->points[i+1], u0);
      nextYS->segments.push_back(nextSeg);
    }
    for (int i=1; i<numSegs(); i++) {
      restYS.segments[i].parallelTransport(restYS.segments[i-1]);
      curYS->segments[i].setU(restYS.segments[i].getU());
      nextYS->segments[i].setU(restYS.segments[i].getU());
    }
    
    for (Segment& s : restYS.segments) {
      s.setFrozen(true);
    }
    
    for (int i=1; i<numCPs()-1; i++) {
      const Segment& ePrev = restYS.segments[i-1];
      const Segment& eNext = restYS.segments[i];
      
      Vec3e curveBinorm = 2.0*ePrev.vec().cross(eNext.vec()) /
      (ePrev.length()*eNext.length() + ePrev.vec().dot(eNext.vec()));
      
      CHECK_NAN_VEC(curveBinorm);
      
      Vec2e restMatCurvePrev(curveBinorm.dot(ePrev.m2()), -(curveBinorm.dot(ePrev.m1())));
      Vec2e restMatCurveNext(curveBinorm.dot(eNext.m2()), -(curveBinorm.dot(eNext.m1())));
      Vec2e restMatCurve = 0.5*(restMatCurvePrev + restMatCurveNext);
      
      rcp.push_back(restMatCurvePrev);
      rcn.push_back(restMatCurveNext);
      rc.push_back(restMatCurve);
      rvl.push_back(0.5*(ePrev.length()+eNext.length()));
    }
    
    // Set mass matrix
    if (masses && masses->rows() == numCPs()) {
      mass.diag = *masses;
      mass.total = mass.diag.sum();
      std::vector<Triplet> triplets;
      triplets.reserve(3*numCPs());
      for (int i=0; i<3*numCPs(); i++) {
        triplets.push_back(Triplet(i, i, (*masses)(i/3)));
      }
      mass.sparse.resize(3*numCPs(), 3*numCPs());
      mass.sparse.setFromTriplets(triplets.begin(), triplets.end());
      
      invMass.diag = mass.diag.cwiseInverse();
      invMass.sparse = mass.sparse.cwiseInverse();
    } else {
      std::cout << "Warning: Mass matrix set to identity.\n";
      mass.diag = VecXe::Ones(numCPs());
      mass.total = numCPs();
      mass.sparse.resize(3*points.size(), 3*points.size());
      mass.sparse.setIdentity();
      
      invMass.diag = VecXe::Ones(numCPs());
      invMass.sparse.resize(3*points.size(), 3*points.size());
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
  const inline size_t numCPs() const { return restYS.points.size(); }
  /// Get the number of control points associated with 2 edges.
  const inline size_t numIntCPs() const { return std::max((int)restYS.points.size()-2, 0); }
  /// Get the number of segments in the yarn.
  const inline size_t numSegs() const { return restYS.segments.size(); }
  
  /// Swaps the current and next yarns, e.g. at the end of a timestep.
  void inline swapYarns() {
    YarnStr* temp = curYS;
    curYS = nextYS;
    nextYS = temp;
    
    for (Segment& s : curYS->segments) {
      s.setFrozen(true);
    }
    for (Segment& s : nextYS->segments) {
      s.setFrozen(false);
    }
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
    for (int i=0; i<numCPs(); i++) {
      com += getMass().diag(i) * cur().points[i].pos;
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
  
  Spline(const CtrlPoint& cp0, const CtrlPoint& cp1, const CtrlPoint& cp2, const  CtrlPoint& cp3) :
  Spline(cp0.pos, cp1.pos, cp2.pos, cp3.pos) { }
  
  Spline(const CtrlPoint& cp0, const CtrlPoint& cp1, const CtrlPoint&cp2, bool start) :
  Spline(start ? cp0.pos : (2*cp0.pos) - cp1.pos, start ? cp1.pos : cp0.pos,
         start ? cp2.pos : cp1.pos, start ? (2*cp2.pos) - cp1.pos : cp2.pos) { }

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
