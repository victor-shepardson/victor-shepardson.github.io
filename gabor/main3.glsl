#extension GL_OES_standard_derivatives : enable

varying vec2 vTextureCoordinates;
varying vec3 position;

uniform float bandwidth; //side length of square grid used for evaluation
uniform float density; //number of impulses / kernel area (accuracy)
uniform vec2 origin; //offset of image space from texture space
uniform vec4 sector; //annular sector in frequency domain: fundamental freq, octaves, orientation, isotropy

uniform vec4 wsector;
uniform float wbandwidth;
uniform float warp, wdirection, wspread;

void main(void){
	vec2 pos = vTextureCoordinates.xy+origin.xy;	
	gnoise_params params;
	params.a = bandwidth;
	params.filterSigma = 1.;
	params.jacob = mat2(dFdx(vTextureCoordinates.xy),dFdy(vTextureCoordinates.xy));
	params.sector = sector;
	params.density = density; //mean impulses per grid cell

	gnoise_params wparams;
	wparams.a = wbandwidth;
	wparams.filterSigma = params.filterSigma;
	wparams.jacob = params.jacob;
	wparams.sector = wsector;
	wparams.density = params.density;
	
	vec2 p0 = vec2(512.,512.);
	float osx = gnoise(pos+p0, wparams);
	//float osy = gnoise(pos-p0, params);
	pos=pos+warp*vec2(cos(wspread*osx+wdirection), sin(wspread*osx+wdirection));

	float value = gnoise(pos, params); 
	value= value*.5+.5;
	
	//monochrome
	vec3 c = vec3(value,value,value);
	//draw fragment
	gl_FragColor = vec4(c, 1.0);
}