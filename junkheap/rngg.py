import itertools

def test(a,b,c,m):
	sum = 0;
	v = 0
	dict = {}
	for i in xrange(m):
		nv = (v*v*a + v*b + c ) % m
		dict[v] = nv
		if (nv==m) or (i==m-1 and nv!=0) or (i<m-1 and nv in dict):
			return (False, dict, sum)
		v = nv
		sum+= abs(v-nv)
	return (True, dict, sum)
	
#print test(34, 1, 0, 289)

for m in xrange(100):
	for (a,b,c) in itertools.product(xrange(m), xrange(m), xrange(m)):
		val = test(a,b,c,m)
		if val[0] and val[2] > .2*m*m:
			print (a, b, c, m)+val
			break;