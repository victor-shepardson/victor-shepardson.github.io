#extension GL_OES_standard_derivatives : enable

#ifdef GL_ES
precision mediump float;
#endif

#define PI 3.1415926535897932384626
#define IMPULSE_CAP 100
#define EPSILON .001
#define C1 4099.
#define C2 4011.

varying vec2 vTextureCoordinates;
varying vec2 position;

uniform float gridSize; //side length of square grid used for evaluation
uniform float density; //number of impulses / kernel area (accuracy)
uniform vec2 origin; //offset of image space from texture space
uniform float sync; //nonrandomness of phase
uniform vec4 harmonic; //annular sector in frequency domain: min freq, max freq, min orientation, max orientation
uniform float filterSigma; //sigma for anisotropic filtering

vec2 pos; //fragment position in image space, measured in pixels
vec2 cpos; //fragment position in cell space: continuous in [0, 1]
ivec2 gpos; //cell position in grid space: integer in [-inf, inf]
mat2 filter, sigma_f_plus_g_inv;
float a, a_prime_square;

float lambda = density/PI; //mean impulses per grid cell

float norm = .33/log2(lambda+2.); //attempt to normalize color value

float nextRand(inout float u){//rng
    u = fract(C2*cos(u + C1 * cos(u)));
	return u;
}
float seed(in vec2 p){
	vec2 temp = p;
	float temp2 = nextRand(temp.x) + temp.y;
	return nextRand(temp2);
}
 
int poisson(inout float u, in float m){//from Galerne, Lagae, Lefebvre, Drettakis
	float u1 = nextRand(u);
	float u2 = nextRand(u);
	float x = sqrt(-2.*log(u1+EPSILON))*cos(2.*PI*u2);
	return int(m+sqrt(m)*x+.5);
}

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
 float eval_cell(in ivec2 dnbr){
	float u = seed(vec2(gpos+dnbr)); //deterministic seed for nbr cell
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
	pos = vTextureCoordinates.xy+origin.xy;
	vec2 temp = pos/gridSize; 
	cpos = fract(temp);
	gpos = ivec2(floor(temp));
	
	a = 1./gridSize;
	
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
	
	
	float value = 
		eval_cell(ivec2(-1, -1)) +
		eval_cell(ivec2(-1,  0)) +
		eval_cell(ivec2(-1,  1)) +
		eval_cell(ivec2( 0, -1)) +
		eval_cell(ivec2( 0,  0)) +
		eval_cell(ivec2( 0,  1)) +
		eval_cell(ivec2( 1, -1)) +
		eval_cell(ivec2( 1,  0)) +
		eval_cell(ivec2( 1,  1));
	
	//normalize / clamp
	value*=norm;
	value= value*.5+.5;
	
	//monochrome
	vec3 c = vec3(value,value,value);//.5*vec3(a_prime_square/(a*a));
	//draw fragment
	gl_FragColor = vec4(c, 1.0);
}