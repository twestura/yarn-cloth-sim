#version 120

uniform sampler2D texture;

varying vec3 N;
varying vec4 eyePosition;
varying vec2 texCoord;

void main()
{
  vec3 normal   = normalize(N);
  vec3 light    = normalize(gl_LightSource[0].position.xyz - eyePosition.xyz);
  vec3 ambient  = gl_FrontMaterial.ambient.xyz * gl_LightSource[0].ambient.xyz;
  vec3 diffuse  = max(dot(normal, light), 0.0) * gl_FrontMaterial.diffuse.xyz*gl_LightSource[0].diffuse.xyz;
  vec4 texColor = texture2D(texture, texCoord);
  
  gl_FragColor = vec4((ambient + diffuse) * texColor.xyz, 1);
}