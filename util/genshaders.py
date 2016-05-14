import sys

# script which converts source code to a javascript file defining a string:
# shadername.glsl -> shadername.glsl.js whuch defines var genshader_shadername
# may be broken for source code containing '\'

def gen(fname):
	lines = [ line.strip() for line in open(fname, 'r').readlines()]
	outname = fname + '.js'
	outfile = open(outname, 'w')
	outfile.write('var genshader_'+outname.split('.')[0]+' = "')
	for line in lines:
		outfile.write(line+'\\n')
	outfile.write('";\n')

if __name__ == "__main__":
	fnames = sys.argv[1:]
	print fnames
	for fname in fnames:
		gen(fname)