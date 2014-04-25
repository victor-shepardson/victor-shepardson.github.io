import sys

# remove links from svg hackily

def gen(fname):
	lines = [ line.strip() for line in open(fname, 'r').readlines()]
	outname = fname+'.stripped' 
	outfile = open(outname, 'w')
	cur_id='lol'
	for line in lines:
		if line.strip()[:2] == '<a' or line.strip()[:3] == '</a':
			pass
		else:
			outfile.write(line+'\n')
			

if __name__ == "__main__":
	fnames = sys.argv[1:]
	print fnames
	for fname in fnames:
		gen(fname)