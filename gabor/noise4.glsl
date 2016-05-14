#extension GL_OES_standard_derivatives : enable

#ifdef GL_ES
precision mediump float;
precision mediump int;
#endif

#define PI 3.1415926535897932384626
#define IMPULSE_CAP 100
#define IMOD 4096

varying vec2 vTextureCoordinates;
varying vec3 position;

uniform float gridSize; //side length of square grid used for evaluation
uniform float density; //number of impulses / kernel area (accuracy)
uniform vec2 origin; //offset of image space from texture space
uniform vec4 sector; //annular sector in frequency domain: fundamental freq, octaves, min orientation, max orientation

uniform vec4 wsector;
uniform float warp;

struct gnoise_params{
	mat2 filter, sigma_f_plus_g_inv;
	float ainv, a, a_prime_square, density, filterSigma, octaves;
	vec4 sector;
	mat2 jacob;
};

ivec2 bound_grid(ivec2 gpos){
	return gpos + IMOD*(1 - gpos/IMOD);
}

//hash based on Blum, Blum & Shub 1986
//and Sharpe http://briansharpe.wordpress.com/2011/10/01/gpu-texture-free-noise/
const float bbsm = 137023.;//magic product of two primes chosen to have high period without float precision issues
vec4 bbsmod( vec4 a ) {
	return a - floor( a * ( 1.0 / bbsm ) ) * bbsm;
}
vec4 bbs(vec4 a) {
	return bbsmod(a*a);
}
vec4 bbsopt( vec4 a ) {
	return fract( a * a * ( 1.0 / bbsm ) ) * bbsm;
}
vec4 seed(in ivec2 p){
	vec4 h = vec4(p.xy, p.xy+ivec2(IMOD/2));
	vec4 h0 = bbs(h);
	vec4 h1 = bbsopt(h0.xzxz);
	vec4 h2 = bbsopt(h0.yyww+h1);
	//vec4 h3 = bbsopt(h0.xzxz+h2);
	return h2*(1./bbsm);
}

//permutation polynomial
//based on Gustavson/McEwan https://github.com/ashima/webgl-noise/
//and Sharpe http://briansharpe.wordpress.com/2011/10/01/gpu-texture-free-noise/
const float pp_epsilon = .01;
vec4 nextRand(inout vec4 u){//rng
	u = fract(((u*34.0*289.)+1.0)*u+pp_epsilon);
	return fract(7.*u);
}

//approximate poisson distribution uniform u
//from Galerne, Lagae, Lefebvre, Drettakis
const float poisson_epsilon = .001;
int poisson(inout vec4 u, in float m){
	float u1 = nextRand(u).x;
	float u2 = nextRand(u).x;
	float x = sqrt(-2.*log(u1+poisson_epsilon))*cos(2.*PI*u2);
	return int(m+sqrt(m)*x+.5);
}

//Gabor noise based on Lagae, Lefebvre, Drettakis, Dutre 2011
 vec4 eval_cell(in vec2 cpos, in ivec2 gpos, in ivec2 dnbr, gnoise_params params){
	vec4 u = seed(bound_grid(gpos+dnbr)); //deterministic seed for nbr cell
	int impulses = int(.5+params.density);//poisson(u, params.density); //number of impulses in nbr cell
	vec4 h = params.sector; //annular sector
	float a = params.a; //bandwidth
	float aps = params.a_prime_square; //intermediate calculations for filtering
	float filt_scale = aps*params.ainv*params.ainv;
	vec2 fpos = cpos - vec2(dnbr);//fragment position in cell space
	
	vec4 acc = vec4(0.);
	//for impulses
	for(int k=0; k<IMPULSE_CAP; k++){
		if(k<impulses){
			//displacement of impulse to fragment
			vec4 delta_x = (fpos.xxxx - nextRand(u))*params.ainv;
			vec4 delta_y = (fpos.yyyy - nextRand(u))*params.ainv;
			vec4 delta_ab = vec4(delta_x.xy, delta_y.xy).xzyw;
			vec4 delta_cd = vec4(delta_x.zw, delta_y.zw).xzyw;
			//impulse frequency, orientation - uniform distribution on input ranges
			vec4 mfreq = pow(vec4(2.), nextRand(u)*params.octaves);
			vec4 ifreq = h.x*mfreq; 
			vec4 iorientation = mix(h.zzzz, h.wwww, nextRand(u));
			//evaluate kernel, accumulate fragment value
			vec4 mu_ab = ifreq*vec4(cos(iorientation.xy), sin(iorientation.xy)).xzyw;
			vec4 mu_cd = ifreq*vec4(cos(iorientation.zw), sin(iorientation.zw)).xzyw;
			vec4 phi = nextRand(u); //phase - uniform dist [0, 1]
			vec4 filt_exp = -.5*vec4(
				dot(mu_ab.xy, params.sigma_f_plus_g_inv*mu_ab.xy),
				dot(mu_ab.zw, params.sigma_f_plus_g_inv*mu_ab.zw),
				dot(mu_cd.xy, params.sigma_f_plus_g_inv*mu_cd.xy),
				dot(mu_cd.zw, params.sigma_f_plus_g_inv*mu_cd.zw)
				);
			vec4 delta_ab2 = delta_ab*delta_ab;
			vec4 delta_cd2 = delta_cd*delta_cd;
			vec4 dd = vec4(
				delta_ab2.xz+delta_ab2.yw,
				delta_cd2.xz+delta_cd2.yw
				);
			vec4 dm = vec4(
				dot(delta_ab.xy,params.filter*mu_ab.xy),
				dot(delta_ab.zw,params.filter*mu_ab.zw),
				dot(delta_cd.xy,params.filter*mu_cd.xy),
				dot(delta_cd.zw,params.filter*mu_cd.zw)
			);
			acc+= filt_scale/mfreq*exp(-PI*aps*dd+filt_exp)*cos(2.*PI*(dm+phi));
		}else {break;}
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
 
 //annular sector of pink noise
vec4 gnoise(vec2 pos, gnoise_params params) {
	//compute positions for this fragment
	vec2 temp = pos*params.a; 
	vec2 cpos = fract(temp);
	ivec2 gpos = bound_grid(ivec2(floor(temp)));
		
	mat2 jacob = params.jacob;
	mat2 jacob_t = mat2(jacob[0][0], jacob[1][0], jacob[0][1], jacob[1][1]);
	mat2 sigma_f_inv = (4.*PI*PI*params.filterSigma*params.filterSigma)*(jacob*jacob_t);
	mat2 sigma_f = inv2x2(sigma_f_inv);
	mat2 sigma_g_inv = (2.*PI*params.ainv*params.ainv)* id2x2();
	mat2 sigma_g = inv2x2(sigma_g_inv);
	mat2 sigma_fg_inv = sigma_f_inv + sigma_g_inv;
	mat2 sigma_fg = inv2x2(sigma_fg_inv);
	
	//filter params
	params.filter = sigma_fg * sigma_g_inv;	
	params.sigma_f_plus_g_inv = inv2x2(sigma_f + sigma_g);
	params.a_prime_square = 2.*PI*sqrt(det2x2(sigma_fg));
	
	params.octaves = params.sector.y;//log2(params.sector.y/params.sector.x);

	vec4 value = 
		eval_cell(cpos, gpos, ivec2(-1, -1), params) +
		eval_cell(cpos, gpos, ivec2(-1,  0), params) +
		eval_cell(cpos, gpos, ivec2(-1,  1), params) +
		eval_cell(cpos, gpos, ivec2( 0, -1), params) +
		eval_cell(cpos, gpos, ivec2( 0,  0), params) +
		eval_cell(cpos, gpos, ivec2( 0,  1), params) +
		eval_cell(cpos, gpos, ivec2( 1, -1), params) +
		eval_cell(cpos, gpos, ivec2( 1,  0), params) +
		eval_cell(cpos, gpos, ivec2( 1,  1), params);
	
	//float lambda = params.a*params.a*params.density;
	
	//ad hoc attempt to normalize
	//value/=log2(2.+2.*params.density);
	value*=(1./PI)*pow(params.density, -.5);
	//value*=(1./3.)/log2(params.density+2.);
	float octexp = pow(2., params.octaves);
	value*= (1.+params.octaves)*octexp/(2.*octexp-1.);
	
	return value;
}
 
 /*
 vec3 hsv2rgb(vec3 c)
{
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}
float rectify( const float f ){
	if( f > 0.0031308 )
		return 1.055 * ( pow( f, 1./2.4 ) ) - 0.055;
	else
		return 12.92 * f;
}
vec3 xyz2rgb(vec3 c){
	return vec3(
		rectify(c.x *  3.2406 + c.y * -1.5372 + c.z * -0.4986),
		rectify(c.x * -0.9689 + c.y *  1.8758 + c.z *  0.0415),
		rectify(c.x *  0.0557 + c.y * -0.2040 + c.z *  1.0570)
		);
}
*/
 
void main(void){
	vec2 pos = vTextureCoordinates.xy+origin.xy;	
	gnoise_params params;
	params.ainv = gridSize;
	params.a = 1./params.ainv;
	params.filterSigma = 1.;
	params.jacob = mat2(dFdx(vTextureCoordinates.xy),dFdy(vTextureCoordinates.xy));
	params.sector = sector;
	params.density = density*(1./PI); //mean impulses per grid cell
/*
	gnoise_params wparams;
	wparams.ainv = params.ainv;
	wparams.a = params.a;
	wparams.filterSigma = params.filterSigma;
	wparams.jacob = params.jacob;
	wparams.sector = wsector;
	wparams.density = params.density;
	*/

	vec4 value = gnoise(pos, params); 
	value= value*.5+.5;
	
	//monochrome
	vec3 c = value.xyz;//.5*vec3(a_prime_square/(a*a));
	//draw fragment
	gl_FragColor = vec4(c, 1.0);
}