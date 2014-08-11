//
//  Integrator.h
//  Visualizer
//
//  Gathers energies (either implicit or explicit) and integrates the given
//  object by a given amount using an IMEX scheme.
//
//  Created by eschweickart on 2/24/14.
//
//

#ifndef Visualizer_IMEXIntegrator_h
#define Visualizer_IMEXIntegrator_h

#include "Integrator.h"
#include "Util.h"
#include "autodiff.h"
#include "Eigen/Sparse"

// An integrator that evaluates both implicit and explicit forces.
class IMEXIntegrator : public Integrator {
private:
  Eigen::SparseMatrix<real> hessBase;
public:
  IMEXIntegrator(std::vector<RodEnergy*>& energies, Rod& r);
  
  bool integrate(Clock& c);
  /// Quasistatically set the twist of the rod. Uses Newton's method to calculate the material
  /// rotation for each edge, assuming the ends are clamped. Returns false if the method does not
  /// converge.
  bool setRotations() const;  
};


#endif
