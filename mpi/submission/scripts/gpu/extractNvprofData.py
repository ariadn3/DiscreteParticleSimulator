import os

unitsMultiplier = {'m': 1e-3}
_OUTPUT_FILE = 'gpuData.csv'
_OUTPUT_PARAM_STRINGS = ('N', 'L', 'r', 'steps', 'config', 'machine')
_CLOCK_TAG = 'time'
_STATS_TAG = ('checkWallCollision', 'checkCollision', 'settleCollision', 'updateParticles')
_STATS_PARAM_STRINGS = (_CLOCK_TAG,) + _STATS_TAG

def getNvprofStats(file):
	stats = {}
	while True:
		line = file.readline()
		if line[0] != '=':
			break

	units = file.readline().strip().split(',')[2]
	adjustmentRatio = unitsMultiplier.get(units[:-1], 1)

	for line in file.readlines():
		line = line.strip()
		if line == '':
			break
		tag = line.split(',')[-1].replace('"', '')[:-6]
		if tag in _STATS_TAG:
			stats[tag] = adjustmentRatio * float(line.split(',')[2])
	return stats

def getTimeStats(file):
	file.readline()
	realString = file.readline().strip().split('\t')[1].split('m')
	realString[1] = realString[1][:-1]
	return int(realString[0])*60 + float(realString[1])

dataDict = {}
for root, dirs, files in os.walk('.', topdown = False):
	for name in files:
		rootSplit = root.split('/')
		if len(rootSplit) < 2 or rootSplit[1] != 'gpu' or name[-5:] == '.time' or name[-4:] == '.out':
			continue

		N, L, r, steps, progType, _ = name.split('-')
		N = int(N)
		L = float(L)
		r = float(r)
		steps = int(steps)
		machine = root.split('/')[2]

		paramTuple = (N, L, r, steps, progType, machine)

		file = open(root + '/' + name, 'r')
		runStats = getNvprofStats(file)
		file.close()

		timeFile = open(root + '/' + name + '.time', 'r')
		runStats[_CLOCK_TAG] = getTimeStats(timeFile)
		timeFile.close()

		if paramTuple not in dataDict or runStats[_CLOCK_TAG] < dataDict[paramTuple][_CLOCK_TAG]:
			dataDict[paramTuple] = runStats

with open(_OUTPUT_FILE, 'w') as fileOut:
	fileOut.write(','.join(_OUTPUT_PARAM_STRINGS) + ',' + ','.join(_STATS_PARAM_STRINGS) + '\n')
	for k, v in dataDict.items():
		fileOut.write(','.join(str(i) for i in k) + ',' + ','.join(str(v[i]) for i in _STATS_PARAM_STRINGS) + '\n')