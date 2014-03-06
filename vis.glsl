#version 100

#ifdef GL_ES
precision mediump float;
#endif

#define EPSILON .001
#define PI 3.1415926535897932384626
#define MAX_SOURCES 2

varying vec2 vTextureCoordinates;

uniform float time[MAX_SOURCES];
uniform vec4 b0[MAX_SOURCES];
uniform vec4 b1[MAX_SOURCES];
uniform float nspeed;
uniform float hspeed;
uniform float magnitude;
uniform float scale;
uniform float harmonic;
uniform float spread;
uniform vec2 mouse[MAX_SOURCES];

//simplex noise from:
	//
	// Description : Array and textureless GLSL 2D simplex noise function.
	//      Author : Ian McEwan, Ashima Arts.
	//  Maintainer : ijm
	//     Lastmod : 20110822 (ijm)
	//     License : Copyright (C) 2011 Ashima Arts. All rights reserved.
	//               Distributed under the MIT License. See LICENSE file.
	//               https://github.com/ashima/webgl-noise
	// 

vec4 mod289(vec4 x)
{
  return x - floor(x * (1.0 / 289.0)) * 289.0;
}

vec3 mod289(vec3 x) {
  return x - floor(x * (1.0 / 289.0)) * 289.0;
}

vec2 mod289(vec2 x) {
  return x - floor(x * (1.0 / 289.0)) * 289.0;
}

vec4 permute(vec4 x){
  return mod289(((x*34.0)+1.0)*x);
}

vec3 permute(vec3 x) {
  return mod289(((x*34.0)+1.0)*x);
}

vec4 taylorInvSqrt(vec4 r)
{
  return 1.79284291400159 - 0.85373472095314 * r;
}

float snoise(vec3 v)
  { 
  const vec2  C = vec2(1.0/6.0, 1.0/3.0) ;
  const vec4  D = vec4(0.0, 0.5, 1.0, 2.0);

// First corner
  vec3 i  = floor(v + dot(v, C.yyy) );
  vec3 x0 =   v - i + dot(i, C.xxx) ;

// Other corners
  vec3 g = step(x0.yzx, x0.xyz);
  vec3 l = 1.0 - g;
  vec3 i1 = min( g.xyz, l.zxy );
  vec3 i2 = max( g.xyz, l.zxy );

  //   x0 = x0 - 0.0 + 0.0 * C.xxx;
  //   x1 = x0 - i1  + 1.0 * C.xxx;
  //   x2 = x0 - i2  + 2.0 * C.xxx;
  //   x3 = x0 - 1.0 + 3.0 * C.xxx;
  vec3 x1 = x0 - i1 + C.xxx;
  vec3 x2 = x0 - i2 + C.yyy; // 2.0*C.x = 1/3 = C.y
  vec3 x3 = x0 - D.yyy;      // -1.0+3.0*C.x = -0.5 = -D.y

// Permutations
  i = mod289(i); 
  vec4 p = permute( permute( permute( 
             i.z + vec4(0.0, i1.z, i2.z, 1.0 ))
           + i.y + vec4(0.0, i1.y, i2.y, 1.0 )) 
           + i.x + vec4(0.0, i1.x, i2.x, 1.0 ));

// Gradients: 7x7 points over a square, mapped onto an octahedron.
// The ring size 17*17 = 289 is close to a multiple of 49 (49*6 = 294)
  float n_ = 0.142857142857; // 1.0/7.0
  vec3  ns = n_ * D.wyz - D.xzx;

  vec4 j = p - 49.0 * floor(p * ns.z * ns.z);  //  mod(p,7*7)

  vec4 x_ = floor(j * ns.z);
  vec4 y_ = floor(j - 7.0 * x_ );    // mod(j,N)

  vec4 x = x_ *ns.x + ns.yyyy;
  vec4 y = y_ *ns.x + ns.yyyy;
  vec4 h = 1.0 - abs(x) - abs(y);

  vec4 b0 = vec4( x.xy, y.xy );
  vec4 b1 = vec4( x.zw, y.zw );

  //vec4 s0 = vec4(lessThan(b0,0.0))*2.0 - 1.0;
  //vec4 s1 = vec4(lessThan(b1,0.0))*2.0 - 1.0;
  vec4 s0 = floor(b0)*2.0 + 1.0;
  vec4 s1 = floor(b1)*2.0 + 1.0;
  vec4 sh = -step(h, vec4(0.0));

  vec4 a0 = b0.xzyw + s0.xzyw*sh.xxyy ;
  vec4 a1 = b1.xzyw + s1.xzyw*sh.zzww ;

  vec3 p0 = vec3(a0.xy,h.x);
  vec3 p1 = vec3(a0.zw,h.y);
  vec3 p2 = vec3(a1.xy,h.z);
  vec3 p3 = vec3(a1.zw,h.w);

//Normalise gradients
  vec4 norm = taylorInvSqrt(vec4(dot(p0,p0), dot(p1,p1), dot(p2, p2), dot(p3,p3)));
  p0 *= norm.x;
  p1 *= norm.y;
  p2 *= norm.z;
  p3 *= norm.w;

// Mix final noise value
  vec4 m = max(0.6 - vec4(dot(x0,x0), dot(x1,x1), dot(x2,x2), dot(x3,x3)), 0.0);
  m = m * m;
  return 42.0 * dot( m*m, vec4( dot(p0,x0), dot(p1,x1), 
                                dot(p2,x2), dot(p3,x3) ) );
  }
/*
float snoise(vec2 v){
  const vec4 C = vec4(0.211324865405187,  // (3.0-sqrt(3.0))/6.0
                      0.366025403784439,  // 0.5*(sqrt(3.0)-1.0)
                     -0.577350269189626,  // -1.0 + 2.0 * C.x
                      0.024390243902439); // 1.0 / 41.0
// First corner
  vec2 i  = floor(v + dot(v, C.yy) );
  vec2 x0 = v -   i + dot(i, C.xx);

// Other corners
  vec2 i1;
  //i1.x = step( x0.y, x0.x ); // x0.x > x0.y ? 1.0 : 0.0
  //i1.y = 1.0 - i1.x;
  i1 = (x0.x > x0.y) ? vec2(1.0, 0.0) : vec2(0.0, 1.0);
  // x0 = x0 - 0.0 + 0.0 * C.xx ;
  // x1 = x0 - i1 + 1.0 * C.xx ;
  // x2 = x0 - 1.0 + 2.0 * C.xx ;
  vec4 x12 = x0.xyxy + C.xxzz;
  x12.xy -= i1;

// Permutations
  i = mod289(i); // Avoid truncation effects in permutation
  vec3 p = permute( permute( i.y + vec3(0.0, i1.y, 1.0 ))
		+ i.x + vec3(0.0, i1.x, 1.0 ));

  vec3 m = max(0.5 - vec3(dot(x0,x0), dot(x12.xy,x12.xy), dot(x12.zw,x12.zw)), 0.0);
  m = m*m ;
  m = m*m ;

// Gradients: 41 points uniformly over a line, mapped onto a diamond.
// The ring size 17*17 = 289 is close to a multiple of 41 (41*7 = 287)

  vec3 x = 2.0 * fract(p * C.www) - 1.0;
  vec3 h = abs(x) - 0.5;
  vec3 ox = floor(x + 0.5);
  vec3 a0 = x - ox;

// Normalise gradients implicitly by scaling m
// Approximation of: m *= inversesqrt( a0*a0 + h*h );
  m *= 1.79284291400159 - 0.85373472095314 * ( a0*a0 + h*h );

// Compute final noise value at P
  vec3 g;
  g.x  = a0.x  * x0.x  + h.x  * x0.y;
  g.yz = a0.yz * x12.xz + h.yz * x12.yw;
  return 130.0 * dot(m, g);
}*/

float aura(vec2 pos){
	float value = 0.;
	for(int i=0; i<MAX_SOURCES; i++){
		vec2 delta = (1./512.)*(pos - mouse[i]);
		float dist = sqrt(delta.x*delta.x+delta.y*delta.y);
		float env = /*sum*/exp(-10.*dist);
		value += env;
	}
	return value;
}

float vis(float t){
// TODO: move in a circle
	vec3 pos = vec3(scale*(1./512.)*gl_FragCoord.xy, t);
	
	vec4 w0 = magnitude*vec4(
		.5*snoise( pos ),
		.66*snoise( 2.*pos ), 
		.75*snoise( 4.*pos ),
		1.*snoise( 8.*pos )
		);
	vec4 w1 = magnitude*vec4(
		1.125*snoise( 16.*pos ),
		1.25*snoise( 32.*pos ),
		1.5*snoise( 64.*pos ),
		2.*snoise( 128.*pos )
		);
	
	float value = 0.;
	
	for(int i=0; i<MAX_SOURCES; i++){
		vec4 b = b0[i]+b1[i];
		float sum = b.x+b.y+b.z+b.w;
			
		vec4 s0 = w0*b0[i];
		vec4 s1 = w1*b1[i];
		vec4 s = s0+s1;
		float warp = s.x+s.y+s.z+s.w;
		
		vec2 delta = (1./512.)*(gl_FragCoord.xy - mouse[i]);
		float dist = sqrt(delta.x*delta.x+delta.y*delta.y);
		float env = exp(-10.*dist);
		value += sum*exp(-10.*(dist+abs(env*warp)));//.5*env*(1.+sin(harmonic*dist - hspeed*t + warp));
	}

	return value;//aura(gl_FragCoord.xy+vec2(value, value));
}

void main(void){
	//compute positions for this fragment
	float t = 0.;
	//interpolate times between positions
	float norm = EPSILON;
	for(int i=0; i<MAX_SOURCES; i++){
		vec2 delta = (1./512.)*(gl_FragCoord.xy - mouse[i]);
		float dist = sqrt(delta.x*delta.x+delta.y*delta.y);
		float env = /*exp(-10.*dist);*/1./(dist+EPSILON);
		t += env * fract(nspeed*time[i]);
		norm += env;
	} t /= norm;
	//select nearest
/*	float mind = 99.;
	for(int i=0; i<MAX_SOURCES; i++){
		vec2 delta = (1./512.)*(gl_FragCoord.xy - mouse[i]);
		float dist = sqrt(delta.x*delta.x+delta.y*delta.y);
		if (dist<mind){
			t = time[i];
			mind=dist;
		}
	}*/
	//t=time[0];
	
	float r = vis(t);
	float g = r;//vis(time-spread*(1.+sin(.13*hspeed*time)));
	float b = r;//vis(time+spread*(1.+sin(.11*hspeed*time)));
	
	//monochrome
	vec3 c = vec3(r,g,b);
	//draw fragment
	gl_FragColor = vec4(c, 1.0);
}