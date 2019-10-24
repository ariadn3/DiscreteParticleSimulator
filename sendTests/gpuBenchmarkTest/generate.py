_DEFAULT_N = 1000
_N = (2000, 3000, 4000, 6000, 8000, 12000, 16000, 24000, 32000)
_DEFAULT_L = 20000
# _L = (5000, 7500, 10000, 15000, 30000, 40000)
_DEFAULT_r = 1
# _r = (2, 3, 4, 6, 8, 16)
_DEFAULT_STEPS = 1000
# _STEPS = (250, 375, 500, 750, 1500, 2000, 3000)

_WRITE_TO_DIR = 'test_in/'

def writeFile(N, L, r, steps):
	with open(_WRITE_TO_DIR + '{}-{}-{}-{}.in'.format(N, L, r, steps), 'w') as testOut:
		testOut.write('{}\n{}\n{}\n{}\nperf\n'.format(N, L, r, steps))

writeFile(_DEFAULT_N, _DEFAULT_L, _DEFAULT_r, _DEFAULT_STEPS)

for n in _N:
	writeFile(n, _DEFAULT_L, _DEFAULT_r, _DEFAULT_STEPS)
# for l in _L:
# 	writeFile(_DEFAULT_N, l, _DEFAULT_r, _DEFAULT_STEPS)
# for r in _r:
# 	writeFile(_DEFAULT_N, _DEFAULT_L, r, _DEFAULT_STEPS)
# for s in _STEPS:
# 	writeFile(_DEFAULT_N, _DEFAULT_L, _DEFAULT_r, s)
