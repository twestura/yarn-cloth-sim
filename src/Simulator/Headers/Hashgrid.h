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
  std::size_t timestamp = 0;
  std::vector<CollisionObject> cos;
};

class HashGrid {
private:
  struct Hasher {
  private:
    const static std::size_t p1 = 73856093;
    const static std::size_t p2 = 19349663;
    const static std::size_t p3 = 83492791;
    static real l;
    
  public:
    inline std::size_t operator()(Vec3e v) const {
      std::size_t x = (std::size_t) (v.x() / l);
      std::size_t y = (std::size_t) (v.y() / l);
      std::size_t z = (std::size_t) (v.z() / l);
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
