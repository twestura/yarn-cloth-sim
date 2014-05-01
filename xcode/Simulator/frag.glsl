#version 120

// Standard diffuse shader. Code originally from Edward Angel's OpenGL: A Primer.

varying vec3 N;
varying vec4 eyePosition;

void main()
{
  vec3 normal  = normalize(N);
  vec3 light   = normalize(gl_LightSource[0].position.xyz - eyePosition.xyz);
  vec3 ambient = gl_FrontMaterial.ambient.xyz * gl_LightSource[0].ambient.xyz;
  vec3 diffuse = max(dot(normal, light), 0.0) * gl_FrontMaterial.diffuse.xyz*gl_LightSource[0].diffuse.xyz;

  
  //vec3 temp = (normal - vec3(0.5, 0.5, 0.5)) * 2;
  gl_FragColor = vec4(ambient + diffuse, 1); //vec4(temp, 1);
}