#version 120

varying vec3 N;
varying vec4 eyePosition;
varying vec2 texCoord;

void main()
{
  gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
  eyePosition = gl_ModelViewMatrix * gl_Vertex;
  texCoord = gl_MultiTexCoord0.xy;
  
  N = normalize(gl_NormalMatrix* normalize(gl_Normal));
}