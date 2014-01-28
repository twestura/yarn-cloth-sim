#version 120

uniform mat3x3 worldProxyJoints[101];
uniform vec3   worldProxyJointsT[101];

attribute vec3 in_Vertex;
attribute vec4 proxyJointIndices;
attribute vec4 proxyJointWeights;

//varying vec4 vcolor;

void main()
{
 // vcolor        = gl_ColorRGBA;
  
  vec3 position = in_Vertex * (1-proxyJointWeights[0]-proxyJointWeights[1]-proxyJointWeights[2]-proxyJointWeights[3]);
  for (int i=0; i<4; i++) {
    mat3x4 pjMat = mat3x4(worldProxyJoints[i], worldProxyJointsT[i]);
    position += (pjMat * vec4(in_Vertex, 1)) * proxyJointWeights[i];
  }
	gl_Position   = gl_ModelViewProjectionMatrix * position;
	
}