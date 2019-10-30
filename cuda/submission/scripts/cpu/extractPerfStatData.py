import os

_OUTPUT_FILE = 'cpuData.csv'
_OUTPUT_PARAM_STRINGS = ('N', 'L', 'r', 'steps', 'threads', 'machine')

def getPerfStatRuntime(file):
	return float(file.readlines()[-2].strip().split()[0])

dataDict = {}
for root, dirs, files in os.walk('.', topdown = False):
	for name in files:
		rootSplit = root.split('/')
		if len(rootSplit) < 2 or rootSplit[1] != 'benchmark':
			continue

		N, L, r, steps, threads, _ = name.split('-')
		N = int(N)
		L = float(L)
		r = float(r)
		steps = int(steps)
		threads = int(threads)
		machine = root.split('/')[2]

		paramTuple = (N, L, r, steps, threads, machine)

		file = open(root + '/' + name, 'r')
		runtime = getPerfStatRuntime(file)
		file.close()

		if paramTuple not in dataDict or runtime < dataDict[paramTuple]:
			dataDict[paramTuple] = runtime

with open(_OUTPUT_FILE, 'w') as fileOut:
	fileOut.write(','.join(_OUTPUT_PARAM_STRINGS) + ',' + 'time' + '\n')
	for k, v in dataDict.items():
		fileOut.write(','.join(str(i) for i in k) + ',' + str(v) + '\n')