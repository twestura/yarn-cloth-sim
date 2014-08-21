//
//  FrameExporter.h
//  Visualizer
//
//  Created by eschweickart on 6/12/14.
//
//

#ifndef Visualizer_FrameExporter_h
#define Visualizer_FrameExporter_h

#include "Constants.h"
#include "Clock.h"

class FrameExporter {
  real fr;
  real lastFrame;
  uint32 frameCount = 0;
public:
  FrameExporter(real framerate = 1.0/60.0) : fr(framerate), lastFrame(-fr) { }
  
  void const inline suggestTimestep(Clock& c) const {
    // NOTE: this may be negative, but clock accounts for this.
    c.suggestTimestep(nextTimestep(c));
  }
  
  real const inline nextTimestep(const Clock& c) const {
    return fr - (c.time() - lastFrame);
  }
  
  void record(const Clock& c) {
    if (c.time() - lastFrame < fr)  return;
    if (c.time() - lastFrame > fr * 1.01) {
      std::cout << "Warning: jumpy frames in Framewriter\n";
    }
    std::stringstream path;
    path << constants::ResultPath << "png/" << frameCount << ".png";
    writeImage(path.str(), ci::app::copyWindowSurface());
    frameCount++;
    lastFrame = c.time();
  }
  
  void writeMPEG(char const* fname) {
    std::string mm = "../../src/RodSound/Source/MovieMaker.sh " + constants::ResultPath + " " + fname;
    system(mm.data());
  }
};

#endif
