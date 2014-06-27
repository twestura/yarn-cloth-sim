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
#include "Yarn.h"
#include "Constants.h"
#include "Clock.h"
#include "PTDetector.h"

enum EvalType {
  Implicit,
  Explicit,
};

enum EnergySource {
  Internal,
  External,
};

class YarnEnergy {
protected:
  const Yarn& y;
  EvalType et;
  std::vector<std::function<void(real)>> frames;
public:
  YarnEnergy(const Yarn& y, EvalType et) : y(y), et(et) { }
  virtual ~YarnEnergy() { }
  virtual bool eval(VecXe*, std::vector<Triplet>* = nullptr, const VecXe* = nullptr) =0;
  virtual void suggestTimestep(Clock&) { }
  const EvalType inline evalType() const { return et; }
  virtual EnergySource inline const energySource() =0;
  void inline setEvalType(EvalType newET) { et = newET; }
  virtual void const draw(real scale = 1.0);
};

class Gravity : public YarnEnergy {
protected:
  Vec3e dir;
public:
  Gravity(const Yarn& y, EvalType et, Vec3e dir);
  bool eval(VecXe* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXe* offset = nullptr);
  EnergySource inline const energySource() { return External; }
};

class Spring : public YarnEnergy {
protected:
  Vec3e clamp;
  size_t index;
  real stiffness;
public:
  Spring(const Yarn& y, EvalType et, size_t index, real stiffness);
  bool eval(VecXe* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXe* offset = nullptr);
  void setClamp(Vec3e newClamp);
  EnergySource inline const energySource() { return External; }
};

class MouseSpring : public YarnEnergy {
private:
  bool mouseDown = false;
  bool mouseSet = false;
  Vec3e mouse;
  size_t index;
  real stiffness;
public:
  MouseSpring(const Yarn& y, EvalType et, size_t index, real stiffness);
  bool eval(VecXe* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXe* offset = nullptr);
  void setMouse(Vec3e newMouse, bool newDown);
  EnergySource inline const energySource() { return External; }
};

class Bending : public YarnEnergy {
private:
  
public:
  Bending(const Yarn& y, EvalType et);
  bool eval(VecXe* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXe* offset = nullptr);
  EnergySource inline const energySource() { return Internal; }
};

class Stretching : public YarnEnergy {
private:
  
public:
  Stretching(const Yarn& y, EvalType et);
  bool eval(VecXe* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXe* offset = nullptr);
  EnergySource inline const energySource() { return Internal; }
};


class Twisting : public YarnEnergy {
private:

public:
  Twisting(const Yarn& y, EvalType et);
  bool eval(VecXe* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXe* offset = nullptr);
  EnergySource inline const energySource() { return Internal; }
};

class IntContact : public YarnEnergy {
private:
  template<typename T>
  static inline T f(T d) {
    return d >= 1.0 ? 0.0 : 1.0/d/d + d*d - 2.0;
  }
  
  template<typename T>
  static inline T df(T d) {
    return d >= 1.0 ? 0.0 : -2.0/d/d/d + 2.0*d;
  }
  
  template<typename T>
  static inline T d2f(T d) {
    return d >= 1.0 ? 0.0 : 6.0/d/d/d/d + 2.0;
  }
  
  const int nb = constants::numQuadPoints;
  const real contactMod = 0.6;
  
public:
  IntContact(const Yarn& y, EvalType et);
  bool eval(VecXe* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXe* offset = nullptr);
  EnergySource inline const energySource() { return Internal; }
};

class PlaneContact : public YarnEnergy {
  Vec3e normal;
  Vec3e origin;
  real stiffness;
public:
  PlaneContact(const Yarn& y, EvalType et, Vec3e normal, Vec3e origin, real stiffness);
  bool eval(VecXe* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXe* offset = nullptr);
  EnergySource inline const energySource() { return External; }
};

class Impulse : public YarnEnergy {
  const Clock& c;
  real start, end;
  Vec3e force;
  size_t index;
public:
  Impulse(const Yarn& y, EvalType et, const Clock& c, real start, real end, Vec3e force,
          size_t index);
  bool eval(VecXe* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXe* offset = nullptr);
  EnergySource inline const energySource() { return External; }
};


#endif
