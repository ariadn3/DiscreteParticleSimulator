_DEFAULT_N = 1000
_N = (250, 500, 2000)
_DEFAULT_L = 20000
_DEFAULT_r = 1
_DEFAULT_STEPS = 1000

def writeFile(N, L, r, steps):
	with open('{}-{}-{}-{}.in'.format(N, L, r, steps), 'w') as testOut:
		testOut.write('{}\n{}\n{}\n{}\nperf\n'.format(N, L, r, steps))

writeFile(_DEFAULT_N, _DEFAULT_L, _DEFAULT_r, _DEFAULT_STEPS)

for n in _N:
	writeFile(n, _DEFAULT_L, _DEFAULT_r, _DEFAULT_STEPS)
