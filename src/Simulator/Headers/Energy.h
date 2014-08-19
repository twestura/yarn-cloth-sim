//
//  Energy.h
//  Visualizer
//
//  Created by eschweickart on 3/7/14.
//
//

#ifndef Visualizer_Energy_h
#define Visualizer_Energy_h

#include "autodiff.h"
#include "Rod.h"
#include "Constants.h"
#include "Clock.h"
#include "PTDetector.h"

/// How the force should be evaluated.
enum EvalType {
  Implicit,
  Explicit,
};

/// The origin of the force relative to the rod.
enum EnergySource {
  Internal,
  External,
};

/// An abstract class representing and energy that acts on the rod.
class RodEnergy {
protected:
  /// A reference to the particular rod this enenergy affects.
  const Rod& r;
  /// How the force should be evaluated.
  EvalType et;
  /// Lambda functions that define how the energy will be drawn, using the timestep as a parameter.
  std::vector<std::function<void(real)>> drawFuncs;
public:
  RodEnergy(const Rod& r, EvalType et) : r(r), et(et) { }
  virtual ~RodEnergy() { }
  
  /// Evaluate the forces generated by the energy based on the current rod configuration.
  virtual bool eval(VecXe*, std::vector<Triplet>* = nullptr, const VecXe* = nullptr) =0;
  /// Returns the EvalType of the energy.
  const EvalType inline evalType() const { return et; }
  /// Returns the energy's source relative to the rod.
  virtual EnergySource inline const energySource() =0;
  /// Sets the EvalType of the energy.
  void inline setEvalType(EvalType newET) { et = newET; }
  /// Draw the energy, optionally passing the timestep as a parameter.
  virtual void const draw(real scale = 1.0);
};

/// A class representing gravity.
class Gravity : public RodEnergy {
protected:
  /// The direction and magnitude of the gravity force.
  Vec3e dir;
public:
  Gravity(const Rod& r, EvalType et, Vec3e dir);
  bool eval(VecXe* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXe* offset = nullptr);
  EnergySource inline const energySource() { return External; }
};

/// A class representing a general spring energy.
class Spring : public RodEnergy {
protected:
  /// The position of the spring clamp (the fixed end of the spring).
  Vec3e clamp;
  /// The index of the point on the rod to which the spring is attached.
  uint32 index;
  /// The stiffness of the spring.
  real stiffness;
public:
  Spring(const Rod& r, EvalType et, uint32 index, real stiffness);
  bool eval(VecXe* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXe* offset = nullptr);
  void setClamp(Vec3e newClamp);
  EnergySource inline const energySource() { return External; }
};

/// A class representing a spring between a point on the rod and the user's mouse.
class MouseSpring : public RodEnergy {
private:
  /// True if the mouse is being held down.
  bool mouseDown = false;
  /// True if the mouse position has been set (ie setMouse has been called)
  bool mouseSet = false;
  /// The 3D position of the mouse.
  Vec3e mouse;
  /// The index of the point on the rod which is affected by the mouse spring.
  uint32 index;
  /// The stiffness of the spring.
  real stiffness;
public:
  MouseSpring(const Rod& r, EvalType et, uint32 index, real stiffness);
  bool eval(VecXe* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXe* offset = nullptr);
  /// Set the 3D mouse position, and record whether the mouse is being held down. Call this every
  /// timestep.
  void setMouse(Vec3e newMouse, bool newDown);
  EnergySource inline const energySource() { return External; }
};

/// The internal bending energy of the rod.
class Bending : public RodEnergy {
private:
  
public:
  Bending(const Rod& r, EvalType et);
  bool eval(VecXe* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXe* offset = nullptr);
  EnergySource inline const energySource() { return Internal; }
};

/// Use the finite element method to simulate stiff-rodlike bending in a single dimension.
class FEMBending : public RodEnergy {
private:
  std::vector<Triplet> modDxxxxTriplets;
  Eigen::SparseMatrix<real> modDxxxx;
  Eigen::SparseMatrix<real> X;
public:
  FEMBending(const Rod& r, EvalType et);
  bool eval(VecXe* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXe* offset = nullptr);
  EnergySource inline const energySource() { return Internal; }
};

/// The internal stretching energy of the rod.
class Stretching : public RodEnergy {
private:
  
public:
  Stretching(const Rod& r, EvalType et);
  bool eval(VecXe* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXe* offset = nullptr);
  EnergySource inline const energySource() { return Internal; }
};

/// The internal twisting energy of the rod.
class Twisting : public RodEnergy {
private:

public:
  Twisting(const Rod& r, EvalType et);
  bool eval(VecXe* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXe* offset = nullptr);
  EnergySource inline const energySource() { return Internal; }
};

/// Energy associated with the rod contacting itself.
class IntContact : public RodEnergy {
private:
  /// The normalized contact function. Takes the relative distance between the two centerlines
  /// (1.0 meaning the rods are just touching, 0.0 meaning the rods are fully compressed into
  /// each other) and returns a scalar stiffness.
  template<typename T>
  static inline T f(T d) {
    return d >= 1.0 ? 0.0 : 1.0/d/d + d*d - 2.0;
  }
  
  /// The derivative of the normalized contact function.
  template<typename T>
  static inline T df(T d) {
    return d >= 1.0 ? 0.0 : -2.0/d/d/d + 2.0*d;
  }
  
  /// The 2nd derivative of the normalized contact function.
  template<typename T>
  static inline T d2f(T d) {
    return d >= 1.0 ? 0.0 : 6.0/d/d/d/d + 2.0;
  }
  
  /// The number of quadrature points per rod edge.
  const int nb = constants::numQuadPoints;
  /// A global tunable contact stiffness.
  const real contactMod = 0.6;
  
public:
  IntContact(const Rod& r, EvalType et);
  bool eval(VecXe* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXe* offset = nullptr);
  EnergySource inline const energySource() { return Internal; }
};

/// Energy associated with the rod contacting a plane. Uses a simple spring model.
class PlaneContact : public RodEnergy {
  /// The normal of the plane.
  Vec3e normal;
  /// A single point on the plane.
  Vec3e origin;
  /// The spring stiffness.
  real stiffness;
public:
  PlaneContact(const Rod& r, EvalType et, Vec3e normal, Vec3e origin, real stiffness);
  bool eval(VecXe* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXe* offset = nullptr);
  EnergySource inline const energySource() { return External; }
};

/// A half-sine impulse applied to a point on the rod.
class Impulse : public RodEnergy {
  /// A reference to the simulation clock.
  const Clock& c;
  /// The start and end times of the impulse.
  real start, end;
  /// The direction and magnitude of the impulse (direction is constant over the duration).
  Vec3e force;
  /// The index of the control point on the rod this impulse affects.
  uint32 index;
public:
  Impulse(const Rod& r, EvalType et, const Clock& c, real start, real end, Vec3e force,
          uint32 index);
  bool eval(VecXe* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXe* offset = nullptr);
  EnergySource inline const energySource() { return External; }
};


#endif
