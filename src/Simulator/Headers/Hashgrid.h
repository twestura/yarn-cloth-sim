//
//  Hashgrid.h
//  Visualizer
//
//  Created by eschweickart on 7/21/14.
//
//

#ifndef __Visualizer__Hashgrid__
#define __Visualizer__Hashgrid__

#include "Clock.h"
#include <unordered_map>

// TODO: Write documentation.

// Avoid abstract class structure to remove allocation, since this needs to be instantiated fast.
class CollisionObject {
public:
  void (*collide)(CollisionObject& other);
  CollisionObject(void (*collide)(CollisionObject&)) : collide(collide) {  }
};

struct GridEntry {
  uint32 timestamp = 0;
  std::vector<CollisionObject> cos;
};

class HashGrid {
private:
  struct Hasher {
  private:
    const static uint32 p1 = 73856093;
    const static uint32 p2 = 19349663;
    const static uint32 p3 = 83492791;
    static real l;
    
  public:
    inline uint32 operator()(Vec3e v) const {
      uint32 x = (uint32) (v.x() / l);
      uint32 y = (uint32) (v.y() / l);
      uint32 z = (uint32) (v.z() / l);
      return (x * p1) ^ (y * p2) ^ (z * p3);
    }
    
    inline static void setSize(real newl) { l = newl; }
    
  };
  
  std::unordered_map<Vec3e, GridEntry, Hasher> grid;
  
public:
  HashGrid(real l) { Hasher::setSize(l); }
  
  void insert(CollisionObject co, Vec3e& v, Clock& c) {
    GridEntry& ge = grid[v]; // inserts a new GridEntry if it does not already exist
    if (ge.timestamp < c.getTicks()) {
      ge.cos.clear();
      ge.timestamp = c.getTicks();
    } else {
      for (CollisionObject& other : ge.cos) {
        co.collide(other);
      }
    }
    
    ge.cos.push_back(co);
  }
  
};


#endif /* defined(__Visualizer__Hashgrid__) */
