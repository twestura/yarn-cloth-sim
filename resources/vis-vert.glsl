#version 120

uniform mat3x3 worldProxyJoints[150];
uniform vec3   worldProxyJointsT[150];

attribute vec4 proxyJointIndices;
attribute vec4 proxyJointWeights;

//varying vec4 vcolor;

void main()
{
 // vcolor        = gl_ColorRGBA;
  
  vec3 position = vec3(0, 0, 0);
  for (int i=0; i<4; i++) {
    mat4x3 pjMat = mat4x3(worldProxyJoints[int(proxyJointIndices[i])], worldProxyJointsT[int(proxyJointIndices[i])]);
//    pjMat = mat4x3(vec3(1, 0, 0), vec3(0, 1, 0), vec3(0, 0, 1), vec3(0, 0, 0));
    position += (pjMat * gl_Vertex) * proxyJointWeights[i];
  }
	gl_Position   = gl_ModelViewProjectionMatrix * vec4(position, 1);
  
}