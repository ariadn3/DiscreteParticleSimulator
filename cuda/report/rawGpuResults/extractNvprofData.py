import os

unitsMultiplier = {'m': 1e-3}
_OUTPUT_FILE = 'gpuData.csv'
_OUTPUT_PARAM_STRINGS = ('N', 'L', 'r', 'steps', 'config', 'machine')

def getNvprofRuntime(file):
	for _ in range(4):
		file.readline()
	units = file.readline().strip().split(',')[2]
	totalTime = 0
	for line in file.readlines():
		line = line.strip()
		if line == '':
			break
		totalTime += float(line.split(',')[2])
	return unitsMultiplier.get(units[:-1], 1) * totalTime

dataDict = {}
for root, dirs, files in os.walk('.', topdown = False):
	for name in files:
		rootSplit = root.split('/')
		if len(rootSplit) < 2 or rootSplit[1] != 'gpu':
			continue

		N, L, r, steps, progType, _ = name.split('-')
		N = int(N)
		L = float(L)
		r = float(r)
		steps = int(steps)
		machine = root.split('/')[2]

		paramTuple = (N, L, r, steps, progType, machine)

		file = open(root + '/' + name, 'r')
		runtime = getNvprofRuntime(file)
		file.close()

		if paramTuple not in dataDict or runtime < dataDict[paramTuple]:
			dataDict[paramTuple] = runtime

with open(_OUTPUT_FILE, 'w') as fileOut:
	fileOut.write(','.join(_OUTPUT_PARAM_STRINGS) + ',' + 'time' + '\n')
	for k, v in dataDict.items():
		fileOut.write(','.join(str(i) for i in k) + ',' + str(v) + '\n')