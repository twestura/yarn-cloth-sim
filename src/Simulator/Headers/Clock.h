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

class Clock
{
private:
  /// Current time of the simulation.
  float t = 0.0f;
  /// Default (maximum) timestep for this clock.
  float defaultTimestep;
  /// Timestep of the model for the next step.
  float h;
  /// Number of times the clock has been incremented.
  size_t ticks = 0;
  
public:
  Clock(float defaultTimestep = constants::INITIAL_TIMESTEP) : defaultTimestep(defaultTimestep),
  h(defaultTimestep) { }
  /// The timestep will not decrease beyond this value.
  const float minTimestep = 1e-6f;
  
  /// Get the current time of the simulation.
  const float inline time() const { return t; }
  /// Get the size of the next timestep.
  const float inline timestep() const { return h; }
  /// Get the number of ticks so far.
  const size_t inline getTicks() const { return ticks; }

  /// Suggest the size of the next timestep.
  void inline suggestTimestep(float s) { h = fmax(minTimestep, fmin(h, s)); }
  /// Increment the timer by its current timestep.
  void inline increment() {
    if (h <= 0.0f) throw;
    t += h;
    h = defaultTimestep;
    ticks++;
  }
  /// Returns true if the timestep can be made smaller.
  bool inline canDecreaseTimestep() const { return h > minTimestep; }
};

#endif
