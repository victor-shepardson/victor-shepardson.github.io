var genshader_main3 = "#extension GL_OES_standard_derivatives : enable\n\nvarying vec2 vTextureCoordinates;\nvarying vec3 position;\n\nuniform float bandwidth; //side length of square grid used for evaluation\nuniform float density; //number of impulses / kernel area (accuracy)\nuniform vec2 origin; //offset of image space from texture space\nuniform vec4 sector; //annular sector in frequency domain: fundamental freq, octaves, orientation, isotropy\n\nuniform vec4 wsector;\nuniform float wbandwidth;\nuniform float warp, wdirection, wspread;\n\nvoid main(void){\nvec2 pos = vTextureCoordinates.xy+origin.xy;\ngnoise_params params;\nparams.a = bandwidth;\nparams.filterSigma = 1.;\nparams.jacob = mat2(dFdx(vTextureCoordinates.xy),dFdy(vTextureCoordinates.xy));\nparams.sector = sector;\nparams.density = density; //mean impulses per grid cell\n\ngnoise_params wparams;\nwparams.a = wbandwidth;\nwparams.filterSigma = params.filterSigma;\nwparams.jacob = params.jacob;\nwparams.sector = wsector;\nwparams.density = params.density;\n\nvec2 p0 = vec2(1.,1.);\nfloat osx = gnoise(pos+p0, wparams);\n//float osy = gnoise(pos-p0, params);\npos=pos+warp*vec2(cos(wspread*osx+wdirection), sin(wspread*osx+wdirection));\n\nfloat value = gnoise(pos, params);\nvalue= value*.5+.5;\n\n//monochrome\nvec3 c = vec3(value,value,value);\n//draw fragment\ngl_FragColor = vec4(c, 1.0);\n}\n";
