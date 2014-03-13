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
  float t = 0;
  /// Timestep of the model for the next step.
  float h = constants::INITIAL_TIMESTEP;
  
public:
  /// Get the current time of the simulation.
  const float inline time() const;
  /// Get the size of the next timestep.
  const float inline timestep() const;
  /// Suggest the size of the next timestep.
  void inline suggestTimestep(float s);
  /// increment the timer by its current timestep.
  void inline increment();
};

const float inline Clock::time() const      { return t; }
const float inline Clock::timestep() const  { return h; }
void inline Clock::suggestTimestep(float s) { h = fmin(h, s); }

void inline Clock::increment() {
  if (h <= 0) throw;
  t += h;
  h = constants::INITIAL_TIMESTEP; // TODO: might increase this as simulation runs
}



#endif
