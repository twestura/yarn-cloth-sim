//
//  ExIntegrator.cpp
//  Visualizer
//
//  Created by eschweickart on 6/10/14.
//
//

#include "ExIntegrator.h"

ExIntegrator::ExIntegrator(Yarn& y, std::vector<YarnEnergy*>& energies) : Integrator(y, energies) {
  alpha1 = 6.615e-7;
  alpha2 = 0.882;
  
  std::vector<Triplet> triplets;
  size_t NumEqs = y.numCPs() * 3;
  for (YarnEnergy* e : energies) {
    if (e->energySource() == Internal) {
      e->eval(nullptr, &triplets);
    }
  }
  stiffness.resize(NumEqs, NumEqs);
  stiffness.setFromTriplets(triplets.begin(), triplets.end());
  stiffness *= -1;
  
  // Eigen::SelfAdjointEigenSolver<Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>> saes(-stiffness.toDense());
  // std::cout << saes.eigenvalues() << "\n\n";
  
  // Rayleigh damping matrix, stiffness not rotated.
  damping = alpha1 * stiffness + alpha2 * y.getMass().sparse;
  
  /*
  affineRest.resize(3*y.numCPs());
  Vec3e com = y.getCurCoM();
  for (int i=0; i<y.numCPs(); i++) {
    affineRest.segment<3>(3*i) = y.rest().points[i].pos - com;
  }
  */
}

bool ExIntegrator::integrate(Clock& c) {
  PROFILER_START("Integrate");
  
  size_t NumEqs = y.numCPs() * 3;
  VecXe forces = VecXe::Zero(NumEqs);
  
  for (YarnEnergy* e : energies) {
    if (!e->eval(&forces)) return false;
  }
  
  // Damping calculations
  /*
  VecXe vel(NumEqs);
  Mat3e Apq = Mat3e::Zero();
  Vec3e com = y.getCurCoM();
  for (int i=0; i<y.numCPs(); i++) {
    Apq += y.getMass().diag(i) * (y.cur().points[i].pos - com) * affineRest.segment<3>(3*i).transpose();
    vel.segment<3>(3*i) = y.cur().points[i].vel;
  }
  Eigen::SelfAdjointEigenSolver<Mat3e> saes;
  saes.compute(Apq.transpose() * Apq);
  Mat3e invSqrt = saes.operatorInverseSqrt();
  if (invSqrt.hasNaN()) { // Polar decomposition failed; matrix is singular.
    forces -= damping * vel;
    if (c.getTicks() % 1000 == 0) {
      std::vector<Triplet> triplets;
      for (YarnEnergy* e : energies) {
        if (e->energySource() == Internal) {
          e->eval(nullptr, &triplets);
        }
      }
      stiffness.setFromTriplets(triplets.begin(), triplets.end());
      stiffness *= -1;
      damping = alpha1 * stiffness + alpha2 * y.getMass().sparse;
    }
  } else {
    Mat3e R;
    R.noalias() = Apq * invSqrt; // rotates affine rest to affine cur
    CHECK_NAN_VEC(R);
    Mat3e RT = R.transpose(); // rotates affine cur to affine rest
    for (int i=0; i<y.numCPs(); i++) {
      vel.segment<3>(3*i) = RT * vel.segment<3>(3*i);
    }
    vel = alpha1 * stiffness * vel;
    for (int i=0; i<y.numCPs(); i++) {
      vel.segment<3>(3*i) = R * vel.segment<3>(3*i);
    }
    vel = alpha2 * y.getMass().sparse * vel;
    forces -= vel;
    
    for (int i=0; i<y.numCPs(); i++) {
      Vec3e aCur = y.cur().points[i].pos - com;
      Vec3e aRest = affineRest.segment<3>(3*i);
      aCur = RT * aCur;
      frames.push_back([aCur, aRest] () {
        for (int i=0; i<8; i++) {
          ci::gl::color(1, 0, 0, 0.5);
          ci::gl::drawSphere(EtoC(aCur) + Vec3c(0, 10, 0), 0.1);
          ci::gl::color(0, 0, 1, 0.5);
          ci::gl::drawSphere(EtoC(aRest) + Vec3c(0, 10, 0), 0.1);
        }
      });
    }
  }
  */
  
  VecXe vel(NumEqs);
  for (int i=0; i<y.numCPs(); i++) {
    vel.segment<3>(3*i) = y.cur().points[i].vel;
  }
  if (c.getTicks() % 1000 == 0 && c.getTicks() != 0) { // May need to update stiffness matrix
    std::vector<Triplet> triplets;
    for (YarnEnergy* e : energies) {
      if (e->energySource() == Internal) {
        e->eval(nullptr, &triplets);
      }
    }
    stiffness.setFromTriplets(triplets.begin(), triplets.end());
    stiffness *= -1;
    damping = alpha1 * stiffness + alpha2 * y.getMass().sparse;
  }
  forces -= damping * vel;
  
  forces = y.getInvMass().sparse * forces;
  
  // Symplectic Euler
  for (int i=0; i<y.numCPs(); i++) {
    Vec3e dqdot = forces.segment<3>(3*i) * c.timestep();
    y.next().points[i].accel = dqdot;
    y.next().points[i].vel = y.cur().points[i].vel + dqdot;
    y.next().points[i].pos = y.cur().points[i].pos + y.next().points[i].vel * c.timestep();
  }
  
  // Explicit Newmark -- NOT STABLE
  /*
  real lambda = 0.5;
  for (int i=0; i<y.numCPs(); i++) {
    Vec3e dqdot = forces.segment<3>(3*i) * c.timestep();
    y.next().points[i].accel = dqdot;
    y.next().points[i].vel = y.cur().points[i].vel + (1.0-lambda) * y.cur().points[i].accel + lambda * dqdot;
    y.next().points[i].pos = y.cur().points[i].pos + c.timestep() * y.cur().points[i].vel + 0.5 * c.timestep() * y.cur().points[i].accel;
  }
   */
  
  
  PROFILER_STOP("Integrate");
  return true;
}