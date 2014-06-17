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

enum EvalType {
  Implicit,
  Explicit,
};

class YarnEnergy {
protected:
  const Yarn& y;
  EvalType et;
  std::vector<std::function<void(void)>> frames;
public:
  YarnEnergy(const Yarn& y, EvalType et) : y(y), et(et) { }
  virtual bool eval(const VecXf&, Clock&, VecXf&, std::vector<Triplet>* = nullptr) =0;
  virtual void suggestTimestep(Clock&) { }
  const EvalType inline evalType() const { return et; }
  void inline setEvalType(EvalType newET) { et = newET; }
  virtual void const draw();
};

class Gravity : public YarnEnergy {
protected:
  Vec3f dir;
public:
  Gravity(const Yarn& y, EvalType et, Vec3f dir);
  bool eval(const VecXf& dqdot, Clock& c, VecXf& Fx, std::vector<Triplet>* GradFx = nullptr);
};

class Spring : public YarnEnergy {
protected:
  Vec3f clamp;
  size_t index;
  float stiffness;
public:
  Spring(const Yarn& y, EvalType et, size_t index, float stiffness);
  bool eval(const VecXf& dqdot, Clock& c, VecXf& Fx, std::vector<Triplet>* GradFx = nullptr);
  void setClamp(Vec3f newClamp);
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
  bool eval(const VecXf& dqdot, Clock& c, VecXf& Fx, std::vector<Triplet>* GradFx = nullptr);
  void setMouse(Vec3f newMouse, bool newDown);
};

class Bending : public YarnEnergy {
private:
  
public:
  Bending(const Yarn& y, EvalType et);
  bool eval(const VecXf& dqdot, Clock& c, VecXf& Fx, std::vector<Triplet>* GradFx = nullptr);
};

class Stretching : public YarnEnergy {
private:
  
public:
  Stretching(const Yarn& y, EvalType et);
  void suggestTimestep(Clock&);
  bool eval(const VecXf& dqdot, Clock& c, VecXf& Fx, std::vector<Triplet>* GradFx = nullptr);
};


class Twisting : public YarnEnergy {
private:

public:
  Twisting(const Yarn& y, EvalType et);
  bool eval(const VecXf& dqdot, Clock& c, VecXf& Fx, std::vector<Triplet>* GradFx = nullptr);
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

  PTDetector* ptd = nullptr;
  
public:
  IntContact(const Yarn& y, EvalType et);
  void suggestTimestep(Clock&);
  bool eval(const VecXf& dqdot, Clock& c, VecXf& Fx, std::vector<Triplet>* GradFx = nullptr);
};

class PlaneContact : public YarnEnergy {
  Vec3f normal;
  Vec3f origin;
  float stiffness;
public:
  PlaneContact(const Yarn& y, EvalType et, Vec3f normal, Vec3f origin, float stiffness);
  bool eval(const VecXf& dqdot, Clock& c, VecXf& Fx, std::vector<Triplet>* GradFx = nullptr);
};

class Impulse : public YarnEnergy {
  float start, end;
  Vec3f force;
  size_t index;
public:
  Impulse(const Yarn& y, EvalType et, float start, float end, Vec3f force, size_t index);
  bool eval(const VecXf& dqdot, Clock& c, VecXf& Fx, std::vector<Triplet>* GradFx = nullptr);
};


#endif
