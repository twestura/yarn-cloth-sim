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
  std::vector<Bvh> children;
  
  int num_nodes() {
    int ret = 1;
    for (Bvh b : children) {
      ret += b.num_nodes();
    }
    return ret;
  }
  
  int num_leaf_nodes() {
    if (children.empty()) return 1;
    int sum = 0;
    for (Bvh c : children) {
      sum += c.num_leaf_nodes();
    }
    return sum;
  }
  
  int num_leaf_indices() {
    if (children.empty()) return indices.size();
    int sum = 0;
    for (Bvh c : children) {
      sum += c.num_leaf_indices();
    }
    return sum;
  }

  
};