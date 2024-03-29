#extension GL_OES_standard_derivatives : enable

#ifdef GL_ES
precision mediump float;
#endif

#define PI 3.1415926535897932384626
#define IMPULSE_CAP 100
#define IMOD 4096

varying vec2 vTextureCoordinates;
varying vec3 position;

uniform float gridSize; //side length of square grid used for evaluation
uniform float density; //number of impulses / kernel area (accuracy)
uniform vec2 origin; //offset of image space from texture space
uniform float sync; //nonrandomness of phase
uniform vec4 harmonic; //annular sector in frequency domain: min freq, max freq, min orientation, max orientation

mat2 filter, sigma_f_plus_g_inv;
float a, a_prime_square;

float lambda = density*(1./PI); //mean impulses per grid cell

ivec2 bound_grid(ivec2 gpos){
	return gpos + IMOD*(1 - gpos/IMOD);
}

//hash based on Blum, Blum & Shub 1986
//and Sharpe http://briansharpe.wordpress.com/2011/10/01/gpu-texture-free-noise/
const float bbsm = 137023.;//magic product of primes chosen to have high period without float precision issues
vec4 bbsmod( vec4 a ) {
	return a - floor( a * ( 1.0 / bbsm ) ) * bbsm;
}
vec4 bbs(vec4 a) {
	return bbsmod(a*a);
}
vec4 bbsopt( vec4 a ) {
	return fract( a * a * ( 1.0 / bbsm ) ) * bbsm;
}
float seed(in ivec2 p){
	vec4 h = vec4(p.xy, p.xy+ivec2(1));
	vec4 h0 = bbs(h);
	vec4 h1 = bbsopt(h0.xzxz);
	vec4 h2 = bbsopt(h0.yyww+h1);
	//vec4 h3 = bbsopt(h0.xzxz+h2);
	return h2.x*(1./bbsm);
}

//permutation polynomial
//based on Gustavson/McEwan https://github.com/ashima/webgl-noise/
//and Sharpe http://briansharpe.wordpress.com/2011/10/01/gpu-texture-free-noise/
const float pp_epsilon = .01;
float nextRand(inout float u){//rng
	u = fract(((u*34.0*289.)+1.0)*u+pp_epsilon);
	return fract(7.*u);
}

//approximate poisson distribution uniform u
//from Galerne, Lagae, Lefebvre, Drettakis
const float poisson_epsilon = .001;
int poisson(inout float u, in float m){
	float u1 = nextRand(u);
	float u2 = nextRand(u);
	float x = sqrt(-2.*log(u1+poisson_epsilon))*cos(2.*PI*u2);
	return int(m+sqrt(m)*x+.5);
}

//Gabor noise based on Lagae, Lefebvre, Drettakis, Dutre 2011

//evaluate the contribution to this fragment by a single impulse at displacement delta, in cell with seed u
float eval_impulse(inout float u, in vec2 delta){
	//impulse frequency, orientation - uniform distribution on input ranges
	float ifreq = mix(harmonic.x, harmonic.y, nextRand(u)); 
	float iorientation = mix(harmonic.z, harmonic.w, nextRand(u));
	//evaluate kernel, accumulate fragment value
	vec2 mu = ifreq*vec2(cos(iorientation), sin(iorientation));
	float phi = nextRand(u); //phase - uniform dist [0, 1]
	float k_prime = a_prime_square/(a*a)*exp(-.5*dot(mu, sigma_f_plus_g_inv*mu));
	return k_prime*exp(-PI*a_prime_square*dot(delta,delta))*cos(2.*PI*(dot(delta, filter*mu)+phi));
	//e ^ -(disp^2 / 2*sigma^2)
}

//evaluate the contribution to this fragment by impulses in the cell at displacement dnbr from this fragment's cell
 float eval_cell(in vec2 cpos, in ivec2 gpos, in ivec2 dnbr){
	float u = seed(bound_grid(gpos+dnbr)); //deterministic seed for nbr cell
	//int u = seed(gpos+dnbr);
	int impulses = poisson(u, lambda); //number of impulses in nbr cell
	float acc = 0.;
	//for impulses
	for(int k=0; k<IMPULSE_CAP; k++){
		if(k>=impulses){
			break;
		}else{ //mysterious bug on windows requires this else
			//position of impulse in cell space - uniform distribution
			vec2 ipos = vec2(nextRand(u), nextRand(u));
			//displacement to fragment
			vec2 delta = (cpos - ipos - vec2(dnbr))*gridSize;
			acc += eval_impulse(u, delta);
		}
	}
	return acc;
 }
 
 float det2x2(mat2 m){
	return (m[0][0]*m[1][1] - m[0][1]*m[1][0]);
 }
 mat2 inv2x2(mat2 m){
	return (1./det2x2(m))*mat2(m[1][1], -m[0][1], -m[1][0], m[0][0]);
 }
 mat2 id2x2(){
	return mat2(1.,0.,0.,1.);
 }
 
void main(void){
	//compute positions for this fragment
	vec2 pos = vTextureCoordinates.xy+origin.xy;
	vec2 temp = pos/gridSize; 
	vec2 cpos = fract(temp);
	ivec2 gpos = bound_grid(ivec2(floor(temp)));
	
	float norm = .33/log2(lambda+2.); //attempt to normalize color value
	
	a = 1./gridSize;
	float filterSigma = 1.;
	
	mat2 jacob = mat2(dFdx(vTextureCoordinates.xy),dFdy(vTextureCoordinates.xy));
	mat2 jacob_t = mat2(jacob[0][0], jacob[1][0], jacob[0][1], jacob[1][1]);
	mat2 sigma_f_inv = (4.*PI*PI*filterSigma*filterSigma)*(jacob*jacob_t);
	mat2 sigma_f = inv2x2(sigma_f_inv);
	mat2 sigma_g_inv = (2.*PI/(a*a))* id2x2();
	mat2 sigma_g = inv2x2(sigma_g_inv);
	mat2 sigma_fg_inv = sigma_f_inv + sigma_g_inv;
	mat2 sigma_fg = inv2x2(sigma_fg_inv);
	
	filter = sigma_fg * sigma_g_inv;	
	sigma_f_plus_g_inv = inv2x2(sigma_f + sigma_g);
	
	a_prime_square = 2.*PI*sqrt(det2x2(sigma_fg));
	/*
	float u = seed(gpos);
	float value = u;
	value = nextRand(u);
	*/
	float value = 
		eval_cell(cpos, gpos, ivec2(-1, -1)) +
		eval_cell(cpos, gpos, ivec2(-1,  0)) +
		eval_cell(cpos, gpos, ivec2(-1,  1)) +
		eval_cell(cpos, gpos, ivec2( 0, -1)) +
		eval_cell(cpos, gpos, ivec2( 0,  0)) +
		eval_cell(cpos, gpos, ivec2( 0,  1)) +
		eval_cell(cpos, gpos, ivec2( 1, -1)) +
		eval_cell(cpos, gpos, ivec2( 1,  0)) +
		eval_cell(cpos, gpos, ivec2( 1,  1));
	
	//normalize / clamp
	value*=norm;
	value= value*.5+.5;
	
	//monochrome
	vec3 c = vec3(value,value,value);//.5*vec3(a_prime_square/(a*a));
	//draw fragment
	gl_FragColor = vec4(c, 1.0);
}