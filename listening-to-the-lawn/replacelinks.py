import sys, string

# remove links from svg hackily

def gen(fname):
	lines = [ line.strip() for line in open(fname, 'r').readlines()]
	outname = fname+'.replaced' 
	outfile = open(outname, 'w')
	cur_id='lol'
	for line in lines:
		line0 = string.replace(line, '<a', '<g')
		line1 = string.replace(line0, '</a', '</g')
		outfile.write(line1+'\n')

if __name__ == "__main__":
	fnames = sys.argv[1:]
	print fnames
	for fname in fnames:
		gen(fname)