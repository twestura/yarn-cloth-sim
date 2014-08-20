//
//  Rod.h
//  Visualizer
//
//  Created by eschweickart on 2/12/14.
//
//

#ifndef __Visualizer__Rod__
#define __Visualizer__Rod__

#include "Util.h"
#include "Constants.h"
#include "CrossSection.h"

/// A rod at a specific point in time. Data that can change between updates is stored here.
class RodSnapshot {
private:
  /// True if the rod data will not change for some time; allow for caching of frequently computed
  /// valus.
  bool frozen = false;
  /// Stores v reference frame components when the rod is frozen (nEdges entries).
  std::vector<Vec3e> frozV;
  /// Stores edges of the rod when the rod is frozen (nEdges entries).
  std::vector<Vec3e> frozEdge;
  
  /// Calculate a 3-vector representing a rod edge for a given edge index.
  const Vec3e inline compEdge(size_t index) const {
    return pos.segment<3>(3*(index+1)) - pos.segment<3>(3*index);
  }
  
  /// Calculate the v component of the reference frame for a given edge index.
  const Vec3e inline compV(size_t index) const {
    return edge(index).cross(u[index]).normalized();
  }
  
public:
  /// Stores the positions of each of the control points of the rod in a numDOF-length vector.
  VecXe pos;
  /// Stores the velocities of each of the control points of the rod in a numDOF-length vector.
  VecXe vel;
  /// Stores the delta-velocities of each of the control points of the rod in a numDOF-length
  /// vector. Notice this is NOT the accelerations; it is the accelerations times the previous
  /// timestep.
  VecXe dVel;
  
  /// Convenience macro to get the 3-vector of the position of a given control point.
#define POS(i) pos.segment<3>(3*(i))
  /// Convenience macro to get the 3-vector of the velocity of a given control point.
#define VEL(i) vel.segment<3>(3*(i))
  /// Convenience macro to get the 3-vector of the delta-velocity of a given control point.
#define DVEL(i) dVel.segment<3>(3*(i))
  
  /// Stores the material frame rotation for each edge of the rod in a numEdges-length vector.
  VecXe rot;
  /// Stores the reference frame twist for each edge of the rod in a numEdges-length vector.
  VecXe refTwist;
  /// Stores the u component of the reference frame for each edge (numEdges entries).
  std::vector<Vec3e> u;
  
  /// Calculate a 3-vector representing an edge of the rod (or return a cached calculation).
  const Vec3e inline edge(size_t index) const {
    if (!frozen) { return compEdge(index); }
    return frozEdge[index];
  }
  
  /// Calculate the v component of the reference frame for a given edge index
  /// (or return a cached calculation).
  const Vec3e inline v(size_t index) const {
    if (!frozen) { return compV(index); }
    return frozV[index];
  }
  
  /// Calculate the length of a given edge.
  const real inline edgeLength(size_t index) const { return edge(index).norm(); }
  
  /// Calculate the material frame vector m1 for a given edge.
  const Vec3e inline m1(size_t index) const {
    return cos(rot(index)) * u[index] + sin(rot(index)) * v(index);
  }
  /// Calculate the material frame vector m2 at a given edge.
  const Vec3e inline m2(size_t index) const {
    return -sin(rot(index)) * u[index] + cos(rot(index)) * v(index);
  }
  
  /// Assures that the rod will not be modified any time soon. This allows certain commonly-used
  /// values to be pre-calculated and cached. Unfreeze the rod before modifying it!
  void inline setFrozen(bool f) {
    frozen = f;
    if (frozen) {
      if (frozEdge.empty() || frozV.empty()) {
        for (int i=0; i<u.size(); i++) {
          frozEdge.push_back(compEdge(i));
          frozV.push_back(compV(i));
        }
      } else {
        for (int i=0; i<u.size(); i++) {
          frozEdge[i] = compEdge(i);
          frozV[i] = compV(i);
        }
      }
    }
  }
  
  /// Parallel transports all reference frames through time, recording the amount of reference twist
  /// accumulated through space. Pass the previous RodSnapshot in time as a parameter.
  void updateReferenceFrames(const RodSnapshot& prevTime) {
    for (int i=0; i<u.size(); i++) {
      // Parallel transport through time
      u[i] = parallelTransport(prevTime.edge(i), edge(i), prevTime.u[i]).normalized();
      // Record the twist in the reference frame (not done for the first edge)
      if (i == 0) continue;
      // First, parallel transport through space
      Vec3e uRef = parallelTransport(edge(i-1), edge(i), u[i-1]).normalized();
      
      // Then record the twist between the two
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

/// Stores information about the mass distribution of a rod.
struct Mass {
public:
  /// Specifies the mass (in kg) of each control point in a numCPs-length vector.
  VecXe diag;
  /// A diagonal mass matrix of size numDOFxnumDOF.
  Eigen::SparseMatrix<real> sparse;
  /// The total mass of the rod.
  real total;
};

/// A class representing a 1D elastic rod. Stores rod data that is constant throughout the
/// simulation.
class Rod {
private:
  /// The number of control points of the rod.
  size_t nCPs;
  /// The number of degrees of freedom of the rod (3*nCPs).
  size_t ndof;
  /// The number of edges of the rod (nCPs-1).
  size_t nEdges;
  
  /// The mass configuration.
  Mass mass;
  /// The inverse mass configuration, cached here to save on computation.
  Mass invMass;
  
  /// Information related to the cross-section of the rod.
  CrossSection* cs;
  
  /// The rod configuration at rest.
  RodSnapshot restRS;
  /// The rod configuration at the current point in time.
  RodSnapshot* curRS;
  /// The rod configuration at the next point in time.
  RodSnapshot* nextRS;
  /// The rest Voronoi lengths defined at each internal control point of the rod.
  std::vector<real> rvl;
  /// The rest curvature with respect to the previous edge (﻿﻿omega bar ^i _i)
  std::vector<Vec2e> rcp;
  /// The rest curvature with respect to the next edge (omega bar ^i _i+1)
  std::vector<Vec2e> rcn;
  /// The rest curvature defined at each internal control point of the rod.
  std::vector<Vec2e> rc;
  
  
public:
  /// The shear modulus of the rod.
  real shearModulus;
  /// The Young's modulus of the rod.
  real youngsModulus;
  
  /// Constructs a rod that is initially at rest given a vector of initial positions.
  /// The u component of the reference frame for the first edge is then propagated down the rod
  /// by parallel transport through space. A vector of lumped masses may be optionally passed in;
  /// if it is null, all masses are set uniformly to 1kg.
  /// Young's and shear moduli may also be specified.
  Rod(VecXe pos, Vec3e u0, VecXe* masses = nullptr,
       real youngs = constants::youngsModulus, real shear = constants::shearModulus) :
  youngsModulus(youngs), shearModulus(shear), ndof(pos.rows()) {
    nCPs = ndof / 3;
    assert(nCPs*3 == ndof && "Malformed position vector");
    nEdges = nCPs - 1;
    
    restRS.pos = pos;
    restRS.vel = VecXe::Zero(ndof);
    restRS.dVel = VecXe::Zero(ndof);
    
    restRS.rot = VecXe::Zero(nEdges);
    restRS.refTwist = VecXe::Zero(nEdges);
    
    // Set (twist-free) reference frames
    restRS.u.push_back(u0);
    for(int i=1; i<nEdges; i++) {
      restRS.u.push_back(parallelTransport(restRS.edge(i-1), restRS.edge(i), restRS.u[i-1]));
    }
    
    // Allocate current and next rod configurations
    curRS  = new RodSnapshot(restRS);
    nextRS = new RodSnapshot(restRS);
    
    // The rest configuration doesn't change, so freeze it
    restRS.setFrozen(true);
    
    // Calculate some properies of the rest configuration
    for (int i=1; i<numCPs()-1; i++) {
      Vec3e ePrev = restRS.edge(i-1);
      Vec3e eNext = restRS.edge(i);
      
      Vec3e curveBinorm = 2.0*ePrev.cross(eNext) / (ePrev.norm()*eNext.norm() + ePrev.dot(eNext));
      
      CHECK_NAN_VEC(curveBinorm);
      
      Vec2e restMatCurvePrev(curveBinorm.dot(restRS.m2(i-1)), -(curveBinorm.dot(restRS.m1(i-1))));
      Vec2e restMatCurveNext(curveBinorm.dot(restRS.m2(i)), -(curveBinorm.dot(restRS.m1(i))));
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
  
  // WARNING: The following methods are safe as long as curRS and nextRS are always allocated upon
  // construction and are never deleted until destruction.
  /// Get a const reference to the rod configuration at the current point in time.
  const RodSnapshot& cur() const  { return *curRS;  }
  /// Get a reference to the rod configuration at the next point in time.
  RodSnapshot& next() { return *nextRS; }
  /// Get a const reference to the rod configuration at the next point in time.
  const RodSnapshot& next() const { return *nextRS; }
  /// Get a const reference to the rod configuration at rest.
  const RodSnapshot& rest() const { return restRS; }
  
  /// Get the number of control points on the rod.
  const inline size_t numCPs() const { return nCPs; }
  /// Get the number of control points associated with 2 edges.
  const inline size_t numIntCPs() const { return nCPs > 2 ? nCPs-2 : 0; }
  /// Get the number of edges in the rod.
  const inline size_t numEdges() const { return nEdges; }
  /// Get the number of degrees of freedom of the rod.
  const inline size_t numDOF() const { return ndof; }
  
  /// Swaps the current and next rod configurations, e.g. at the end of a timestep.
  void inline swapRods() {
    RodSnapshot* temp = curRS;
    curRS = nextRS;
    nextRS = temp;
    
    curRS->setFrozen(true);
    nextRS->setFrozen(false);
  }
  
  /// Get the rest Voronoi length for an internal control point.
  const inline real restVoronoiLength(size_t index) const {
    assert(index > 0 && "Voronoi length undifined at this control point.");
    return rvl[index-1];
  }
  
  /// Get the rest curvature for an internal control point.
  const inline Vec2e& restCurvePrev(size_t index) const {
    assert(index > 0 && "Curvature undefined at this control point.");
    return rcp[index-1];
  }
  
  /// Get the rest curvature for an internal control point.
  const inline Vec2e& restCurveNext(size_t index) const {
    assert(index > 0 && "Curvature undefined at this control point.");
    return rcn[index-1];
  }
  
  /// Get the rest curvature for an internal control point.
  const inline Vec2e& restCurve(size_t index) const {
    assert(index > 0 && "Curvature undefined at this control point.");
    return rc[index-1];
  }
  
  /// Get the cross-section of the rod.
  const inline CrossSection& getCS() const { return *cs; }
  
  /// Get the rod's mass configuration.
  const inline Mass& getMass() const { return mass; }
  /// Get the rod's inverse mass configuration.
  const inline Mass& getInvMass() const { return invMass; }
  
  /// Get the center of mass of the current rod.
  const Vec3e getCurCoM() const {
    Vec3e com = Vec3e::Zero();
    VecXe posTimesMass = getMass().sparse * cur().pos;
    for (int i=0; i<numCPs(); i++) {
      com += posTimesMass.segment<3>(3*i);
    }
    return com / getMass().total;
  }
  
  ~Rod() {
    delete curRS;
    delete nextRS;
  }
};

/// A class representing a Catmull-Rom spline segment parameterized by 4 3D points.
class Spline {
protected:
  /// The points that define the spline.
  Vec3e p[4];
  
public:
  /// Construct the spline from 4 control points.
  Spline(const Vec3e p0, const Vec3e p1, const Vec3e p2, const Vec3e p3) {
    p[0] = p0; p[1] = p1; p[2] = p2; p[3] = p3;
  }

  /// Returns a point on the spline given a time value between 0.0 and 1.0.
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

#endif /* defined(__Visualizer__Rod__) */
