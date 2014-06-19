//
//  Energy.h
//  Visualizer
//
//  Created by eschweickart on 3/7/14.
//
//

#ifndef Visualizer_Energy_h
#define Visualizer_Energy_h

#include "Eigen/Sparse"
#include "autodiff.h"
#include "Yarn.h"
#include "Constants.h"
#include "Clock.h"
#include "PTDetector.h"

typedef Eigen::Triplet<float> Triplet;
typedef Eigen::VectorXf VecXf;
typedef Eigen::Matrix3f Mat3f;

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
  std::vector<std::function<void(float)>> frames;
public:
  YarnEnergy(const Yarn& y, EvalType et) : y(y), et(et) { }
  virtual bool eval(VecXf*, std::vector<Triplet>* = nullptr, const VecXf* = nullptr) =0;
  virtual void suggestTimestep(Clock&) { }
  const EvalType inline evalType() const { return et; }
  virtual EnergySource inline const energySource() =0;
  void inline setEvalType(EvalType newET) { et = newET; }
  virtual void const draw(float scale = 1.0f);
};

class Gravity : public YarnEnergy {
protected:
  Vec3f dir;
public:
  Gravity(const Yarn& y, EvalType et, Vec3f dir);
  bool eval(VecXf* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXf* offset = nullptr);
  EnergySource inline const energySource() { return External; }
};

class Spring : public YarnEnergy {
protected:
  Vec3f clamp;
  size_t index;
  float stiffness;
public:
  Spring(const Yarn& y, EvalType et, size_t index, float stiffness);
  bool eval(VecXf* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXf* offset = nullptr);
  void setClamp(Vec3f newClamp);
  EnergySource inline const energySource() { return External; }
};

class MouseSpring : public YarnEnergy {
private:
  bool mouseDown = false;
  bool mouseSet = false;
  Vec3f mouse;
  size_t index;
  float stiffness;
public:
  MouseSpring(const Yarn& y, EvalType et, size_t index, float stiffness);
  bool eval(VecXf* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXf* offset = nullptr);
  void setMouse(Vec3f newMouse, bool newDown);
  EnergySource inline const energySource() { return External; }
};

class Bending : public YarnEnergy {
private:
  
public:
  Bending(const Yarn& y, EvalType et);
  bool eval(VecXf* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXf* offset = nullptr);
  EnergySource inline const energySource() { return Internal; }
};

class Stretching : public YarnEnergy {
private:
  
public:
  Stretching(const Yarn& y, EvalType et);
  bool eval(VecXf* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXf* offset = nullptr);
  EnergySource inline const energySource() { return Internal; }
};


class Twisting : public YarnEnergy {
private:

public:
  Twisting(const Yarn& y, EvalType et);
  bool eval(VecXf* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXf* offset = nullptr);
  EnergySource inline const energySource() { return Internal; }
};

class IntContact : public YarnEnergy {
private:
  static inline float f(float d) {
    return d >= 1.0f ? 0.0f : 1.0f/d/d + d*d - 2.0f;
  }
  
  static inline float df(float d) {
    return d >= 1.0f ? 0.0f : -2.0f/d/d/d + 2.0f*d;
  }
  
  static inline float d2f(float d) {
    return d >= 1.0f ? 0.0f : 6.0f/d/d/d/d + 2.0f;
  }
  
  const int nb = constants::numQuadPoints;
  const float contactMod = 0.6f;
  
public:
  IntContact(const Yarn& y, EvalType et);
  bool eval(VecXf* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXf* offset = nullptr);
  EnergySource inline const energySource() { return Internal; }
};

class PlaneContact : public YarnEnergy {
  Vec3f normal;
  Vec3f origin;
  float stiffness;
public:
  PlaneContact(const Yarn& y, EvalType et, Vec3f normal, Vec3f origin, float stiffness);
  bool eval(VecXf* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXf* offset = nullptr);
  EnergySource inline const energySource() { return External; }
};

class Impulse : public YarnEnergy {
  const Clock& c;
  float start, end;
  Vec3f force;
  size_t index;
public:
  Impulse(const Yarn& y, EvalType et, const Clock& c, float start, float end, Vec3f force,
          size_t index);
  bool eval(VecXf* Fx, std::vector<Triplet>* GradFx = nullptr, const VecXf* offset = nullptr);
  EnergySource inline const energySource() { return External; }
};


#endif
