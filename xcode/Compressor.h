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

#endif /* defined(__Visualizer__Compressor__) */

struct Bvh {
  ci::Vec3f aabb[2];
  std::vector<uint32_t> indices;
  
};