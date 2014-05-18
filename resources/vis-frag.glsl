#version 120

//varying vec4 vcolor;

/*
float LinearizeDepth(float zoverw)
{
	float n = 1.0;
	float f = 10000.0;
	
	return(2.0 * n) / (f + n - zoverw * (f - n));
}

uniform sampler2D depthImg;
*/

void main()
{
  
	//float depth = texture2D(depthImg, gl_TexCoord[0].xy).r;
	//depth = LinearizeDepth(depth)*77;
  
	// gl_FragColor = vcolor;
  gl_FragColor = vec4(0, 0.4, 1, 0.6);
}



