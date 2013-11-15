#ifdef GL_ES
precision mediump float;
#endif

#define PROCESSING_COLOR_SHADER

#define PI 3.1415926535897932384626
#define IMPULSE_CAP 1000
#define EPSILON .001

uniform float gridSize; //side length of square grid used for evaluation
uniform float density; //number of impulses / kernel area (accuracy)
uniform vec2 origin; //offset of image space from texture space
uniform int offset; //offset to rng seeds
uniform float sync; //nonrandomness of phase
uniform vec4 harmonic; //annular sector in frequency domain: min freq, max freq, min orientation, max orientation

//mean impulses per grid cell
float lambda = density/PI;

float expcoeff = -PI/(gridSize*gridSize);
float norm = .33/log2(lambda);

//Borosh and Niederreiter 1983
uint nextRand(uint lastRand){//rng
    return 3039177861U*lastRand;
}
float consumeFloat(inout uint r){ //return a rand float in [0, 1), advance rng
	//float ret = float(r)/4.2949673E9;
	float ret = float(r) * 2.3283064E-10; //precise enough?
	r = nextRand(r);
	return ret;
}

uint rowmajorSeed(ivec2 pos){
	uint ret = (uint(pos.x) & 0x0000ffffU) | (uint(pos.y) << 16);
	//for(int i=0; i<2; i++)
	//	ret = nextRand(ret);
	return ret + uint(offset);
}

uint mortonSeed(ivec2 pos){
	//morton order seed,
	//interleave bits
	int mask = 1;
	int ret = 0;
	pos.y = pos.y << 1;
	for(int i=0; i<16; i++){
		ret = ret | (pos.x & mask);
		mask = mask << 1;
		ret = ret | (pos.y & mask);
		mask = mask << 1;
		pos.x = pos.x << 1;
		pos.y = pos.y << 1;
	}
	return uint(ret+offset);
 }
 
 float apprexp(float x){
	/*float t1 = .5*x;
	float t2 = .1*x*x;
	float t3 = .0833333333*x*t2;
	float t_o = t1+t3;
	float t_e = 1.0+t2;
	return max(0.0,(t_e+t_o)/(t_e-t_o));*/
	float x2 = x*x*.5;
	float x3 = x2*x*.333333333333;
	return max(0.0, 1.0+x+x2+x3);
 
 }
 
 int poisson(inout uint u, float m){//from Galerne, Lagae, Lefebvre, Drettakis
	float u1 = consumeFloat(u);
	float u2 = consumeFloat(u);
	float x = sqrt(-2*log(u1+EPSILON))*cos(2*PI*u2);
	return int(m+sqrt(m)*x+.5);
 
 }
void main(void){
	//fragment position in image space: continuous [-inf, inf]
	vec2 pos = gl_FragCoord.xy+origin.xy;

	vec2 temp; 
	vec2 cpos = modf(pos/gridSize, temp); //fragment position in cell space: continuous [0, 1]
	//correct for negative coordinates
	vec2 mc = vec2(cpos.x<0, cpos.y<0);
	cpos+= mc; temp-= mc;
	ivec2 gpos = ivec2(temp); //cell position in grid space: integer in [-inf, inf]
	
	float v =0; //fragment value
	
	//for local cells
	ivec2 dnbr; //offset to nbr cell in grid space
	for(dnbr.x=-1; dnbr.x<=1; dnbr.x++){
	for(dnbr.y=-1; dnbr.y<=1; dnbr.y++){
		uint u = mortonSeed(gpos+dnbr); //deterministic seed for nbr cell
		//uint u = nextRand(s);
		int impulses = poisson(u, lambda); //number of impulses in nbr cell
		int k=0;
		//for impulses
		for(int k=0; k<IMPULSE_CAP; k++){
			if(k>=impulses) break;
			//position of impulse in cell space - uniform distribution
			vec2 ipos = vec2(consumeFloat(u), consumeFloat(u));
			//displacement to fragment
			vec2 delta = (cpos - ipos - vec2(dnbr))*gridSize;
			//impulse frequency, orientation - uniform distribution on input ranges
			float ifreq = mix(harmonic.x, harmonic.y, consumeFloat(u)); 
			float iorientation = mix(harmonic.z, harmonic.w, consumeFloat(u));
			//evaluate kernel, accumulate fragment value
			vec2 omega = vec2(cos(iorientation), sin(iorientation));
			//phase - uniform dist [0, 1]
			float corientation = mix(harmonic.z, harmonic.w, .5);
			float null;
			vec2 disp = pos-delta; //position of impulse in image space
			float phi = modf(dot(vec2(cos(corientation), sin(corientation)),disp)*ifreq, null);//consumeFloat(u);
			phi = mix(consumeFloat(u), phi, sync);
			v+= exp(dot(delta,delta)/*expcoeff)*/*-PI*ifreq*ifreq)
				*cos(2*PI*(ifreq*dot(delta, omega)+phi));
		}
	}}
	//normalize / clamp
	v*=norm;
	v=v*.5+.5;
	//monochrome
	vec3 c = vec3(v,v,v);
	//draw fragment
	gl_FragColor = vec4(c, 1.0);
}