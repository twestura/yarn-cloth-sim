//
//  Compressor.h
//  Visualizer
//
//  Created by eschweickart on 12/16/13.
//
//

#ifndef __Visualizer__Compressor__
#define __Visualizer__Compressor__

#include <iostream>
#include "Common.h"

// A group of control points within an AABB.
struct AAGroup {
  ci::Vec3f aabb[2];
  std::vector<uint32_t> indices;
  bool ok = false;
};

// Compress a range of frames. Create a .comp file if possible; otherwise create a .key file.
void compress(AppData& ad);

// Decompress a given frame (.comp or .key) and load it into
// the circular buffer.
void decompressFrame(AppData& ad, const int frame);

#endif /* defined(__Visualizer__Compressor__) */

