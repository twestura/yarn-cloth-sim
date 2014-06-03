#version 120

// Standard diffuse shader. Code originally from Edward Angel's OpenGL: A Primer.

varying vec3 N;
varying vec4 eyePosition;

void main()
{
  gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
  eyePosition = gl_ModelViewMatrix * gl_Vertex;
  
  
  N = normalize(gl_NormalMatrix* normalize(gl_Normal));
}