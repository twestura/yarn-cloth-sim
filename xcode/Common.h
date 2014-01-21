//
//  Common.h
//  Visualizer
//
//  Created by eschweickart on 1/7/14.
//
//

#ifndef Visualizer_common_h
#define Visualizer_common_h

#include <boost/circular_buffer.hpp>

using namespace std;
using namespace ci;

//#define CREATE_POS_FILES
//#define CREATE_RESID_FILES
//#define COMPRESS
//#define DECOMPRESS
//#define RADIUS_ONE_HALF

#define MIN_RES 0
#ifdef RADIUS_ONE_HALF
#define RADIUS 0.5
#define MAX_RES 0.08 // TODO: This should not be hard-coded.
#else
#define RADIUS 1
#define MAX_RES 2.07
#endif // ifdef RADIUS_ONE_HALF

#define NUM_FRAMES 7
#define START_FRAME 0
#define END_FRAME 70

#define RESID_PATH "../../../../afghan/resid/"
#define POS_PATH "../../../../afghan/pos/"

// A struct representing one frame of the animation.
// Optionally stores residual information.
struct Frame {
  float minResid;
  float maxResid;
  vector<float> resid;
  vector<Vec3f> pos;
  
  Frame () {}
  
  Frame (vector<Vec3f> p) {
    pos = p;
    resid.clear();
    minResid = maxResid = 0;
  }
  
  Frame (vector<Vec3f> p, vector<float> r, float min, float max) {
    pos = p;
    resid = r;
    minResid = min;
    maxResid = max;
  }
  
  Vec3f& operator[]( const uint32_t index )
  {
    return pos[index];
  }
  
  const Vec3f& operator[]( const uint32_t index ) const
  {
    return pos[index];
  }
  
  const uint32_t size() const
  {
    return pos.size();
  }
};

// A struct for use in the KD-tree. Finds all control points within
// a given squared distance.
struct NeighborLookupProc {
  vector<uint32_t> neighbors;
  void process( uint32_t id, float distSqrd, float &maxDistSqrd )
  {
    neighbors.push_back(id);
  }
};

// A modified circular buffer to load in frames dynamically. Always allows
// access to the first element of the buffer.
template <class T>
struct ModCircularBuffer {
private:
  bool is_empty = true;
  int offset = 0;
  T zero;
  boost::circular_buffer<T> cb;
public:
  
  void init( const uint32_t size, int start_frame )
  {
    cb = boost::circular_buffer<T>(size);
    offset = (start_frame == 0 ? 0 : -start_frame+1);
  }
  
  T& operator[]( const uint32_t index )
  {
    if (index == 0) {
      return zero;
    } else {
      return cb[(index-1) + offset];
    }
  }
  
  const T& operator[]( const uint32_t index ) const
  {
    if (index == 0) {
      return zero;
    } else {
      return cb[(index-1) + offset];
    }
  }
  
  void push_back( T elt )
  {
    if (is_empty) {
      zero = elt;
      is_empty = false;
    } else {
      if (cb.full()) offset--;
      cb.push_back(elt);
    }
  }
  
  void push_front( T elt )
  {
    if (cb.full()) offset++;
    cb.push_front(elt);
  }
  
  const int left_buffer_size( int frame ) const
  {
    return cb.capacity() - frame - offset;
  }
  
  const int right_buffer_size( int frame ) const
  {
    return frame-1 + offset;
  }
  
};

struct AppData {
  ModCircularBuffer<Frame> frames;
  int currentFrame = START_FRAME;
};

// Given a group of points, compute the residual from the least squares
// best fit transformation from the current frame to the next. If retMax
// is true, the maximum error is returned; otherwise, the total error squared
// is returned.
float getResidual(const AppData& ad, const vector<uint32_t>& indices, const int curFrame, const int nextFrame, const bool retMax);

// Similar to the method above, in an attempt to get things to be more accurate.
// Not entirely successful.
/*
float newGetResidual(const vector<uint32_t> indices, const bool retMax)
{
  using namespace Eigen;
  int n = indices.size();
  
  // Compute centroids
  Vec3f c = Vec3f();
  Vec3f cTilde = Vec3f();
  for (int index : indices) {
    c += points[currentFrame][index];
    cTilde += points[currentFrame][index]; //here
  }
  c /= n;
  cTilde /= n;
  
  // Create least squares system--see James and Twigg, Appendix A.
  Matrix<float, 12, 12> M;
  M << cTilde[0]*cTilde[0], cTilde[1]*cTilde[0], cTilde[2]*cTilde[0], 0, 0, 0, 0, 0, 0, cTilde[0], 0, 0,
  cTilde[0]*cTilde[1], cTilde[1]*cTilde[1], cTilde[2]*cTilde[1], 0, 0, 0, 0, 0, 0, 0, cTilde[0], 0,
  cTilde[0]*cTilde[2], cTilde[1]*cTilde[2], cTilde[2]*cTilde[2], 0, 0, 0, 0, 0, 0, 0, 0, cTilde[0],
  0, 0, 0, cTilde[0]*cTilde[0], cTilde[1]*cTilde[0], cTilde[2]*cTilde[0], 0, 0, 0, cTilde[1], 0, 0,
  0, 0, 0, cTilde[0]*cTilde[1], cTilde[1]*cTilde[1], cTilde[2]*cTilde[1], 0, 0, 0, 0, cTilde[1], 0,
  0, 0, 0, cTilde[0]*cTilde[2], cTilde[1]*cTilde[2], cTilde[2]*cTilde[2], 0, 0, 0, 0, 0, cTilde[1],
  0, 0, 0, 0, 0, 0, cTilde[0]*cTilde[0], cTilde[1]*cTilde[0], cTilde[2]*cTilde[0], cTilde[2], 0, 0,
  0, 0, 0, 0, 0, 0, cTilde[0]*cTilde[1], cTilde[1]*cTilde[1], cTilde[2]*cTilde[1], 0, cTilde[2], 0,
  0, 0, 0, 0, 0, 0, cTilde[0]*cTilde[2], cTilde[1]*cTilde[2], cTilde[2]*cTilde[2],  0, 0, cTilde[2],
  cTilde[0], cTilde[1], cTilde[2], 0, 0, 0, 0, 0, 0, 1, 0, 0,
  0, 0, 0, cTilde[0], cTilde[1], cTilde[2], 0, 0, 0, 0, 1, 0,
  0, 0, 0, 0, 0, 0, cTilde[0], cTilde[1], cTilde[2], 0, 0, 1;
  
  Matrix<float, 12, 1> B;
  B << c[0]*cTilde[0], c[1]*cTilde[0], c[2]*cTilde[0], c[0]*cTilde[1], c[1]*cTilde[1], c[2]*cTilde[1],
  c[0]*cTilde[2], c[1]*cTilde[2], c[2]*cTilde[2], c[0], c[1], c[2];
  
  // Solve least squares, convert to affine transform matrix
  Matrix<float, 12, 1> X = M.fullPivLu().solve(B);
  
  Matrix<float, 3, 4> Transform;
  for (int i=0; i<9; i++) {
    Transform(i/3, i%3) = X(i);
  }
  Transform(0,3) = X(9);
  Transform(1,3) = X(10);
  Transform(2,3) = X(11);
  
  cout << Transform << "\n";
  
  // Compute residuals
  float ret = 0;
  Vec3f temp;
  Vector4f input;
  Vector3f output;
  for (int index : indices) {
    temp = points[currentFrame][index] - c;
    input << temp.x, temp.y, temp.z, 1;
    output = Transform * input;
    for (int j=0; j<3; j++) {
      float diff = abs(output(j) + cTilde[j] - points[currentFrame][index][j]); //here
      if (retMax) {
        ret = max(ret, diff);
      } else {
        ret += diff * diff;
      }
    }
  }
  
  return ret;
}
*/

// Load a .pos file into the circular buffer.
void loadFrame(AppData& ad, const int frame, const bool back);

// A fairly unsafe method that writes contiguous memory to a given file.
// Abstracted here to confine the dragons and black magic.
void writeBinary(const void* data, const uint size, ofstream& outFile);

#endif
