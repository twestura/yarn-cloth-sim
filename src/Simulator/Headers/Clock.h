//
//  Clock.h
//  Visualizer
//
//  Created by eschweickart on 2/24/14.
//
//

#ifndef Visualizer_Clock_h
#define Visualizer_Clock_h

#include "Constants.h"

/// A class that controls the passage of time within a simulation. Decides and keeps track of the
/// timestep, and records how much time has passed since the simulation began.
class Clock
{
private:
  /// Current time of the simulation in seconds.
  real t = 0.0;
  /// Default (maximum) timestep for this clock.
  real defaultTimestep;
  /// Timestep of the model for the next step.
  real h;
  /// Number of times the clock has been incremented.
  size_t ticks = 0;
  
public:
  Clock(real defaultTimestep = constants::INITIAL_TIMESTEP) : defaultTimestep(defaultTimestep),
  h(defaultTimestep) { }
  /// The timestep will not decrease beyond this value.
  const real minTimestep = 1.0e-6;
  
  /// Get the current time of the simulation in seconds.
  const real inline time() const { return t; }
  /// Get the size of the next timestep.
  const real inline timestep() const { return h; }
  /// Get the number of times the clock has been incremented.
  const size_t inline getTicks() const { return ticks; }

  /// Suggest the size of the next timestep. If larger than the current timestep, the request
  /// will be ignored. The size of the timestep will not decrease beyond the minimum timestep.
  void inline suggestTimestep(real s) { h = std::max(minTimestep, std::min(h, s)); }
  /// Increment the timer by its current timestep.
  void inline increment() {
    if (h <= 0.0) throw;
    t += h;
    h = defaultTimestep;
    ticks++;
  }
  /// Returns true if the timestep can be made smaller.
  bool inline canDecreaseTimestep() const { return h > minTimestep; }
};

#endif
