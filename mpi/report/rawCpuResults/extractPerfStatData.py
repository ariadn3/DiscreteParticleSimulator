import os
import re

_INTERESTED_STATS = {'task-clock': float,
						'context-switches': int,
						'cycles': int,
						'instructions': int,
						'branches': int,
						'branch-misses':int
						}
_CLOCK_DTYPE = float
_CLOCK_STRING = 'wall-clock-time'

_OUTPUT_FILE = 'parData.csv'
_PARAM_STRINGS = ('N', 'L', 'r', 'steps', 'processes', 'configType')

dataStruct = re.compile('^.+?(?=#)')
dataDict = {}
for root, dirs, files in os.walk('.', topdown = False):
	for name in files:
		configType = root.split('/')[-1]
		if configType != 'PureXeS' and configType != 'PureI7':
			continue

		N, L, r, steps, processes, _ = name.split('-')
		N = int(N)
		L = float(L)
		r = float(r)
		steps = int(steps)
		processes = int(processes)
		paramTuple = (N, L, r, steps, processes, configType)

		with open(root + '/' + name) as fileOut:
			lines = fileOut.readlines()
			wallClockTime = _CLOCK_DTYPE(lines[-2].split()[0])

			if paramTuple not in dataDict or wallClockTime < dataDict[paramTuple][_CLOCK_STRING]:
				curStats = {_CLOCK_STRING: wallClockTime}
				for line in lines:
					line = line.strip()
					regexSearch = dataStruct.findall(line)
					if regexSearch:
						regexSearch = regexSearch[0].strip().replace(',', '')
						linePart = regexSearch.split()
						if any(stat in linePart for stat in _INTERESTED_STATS):
							curStats[linePart[1]] = _INTERESTED_STATS[linePart[1]](linePart[0])
				dataDict[paramTuple] = curStats

with open(_OUTPUT_FILE, 'w') as outFile:
	outFile.write(','.join(_PARAM_STRINGS) + ',' + _CLOCK_STRING + ',' + ','.join(_INTERESTED_STATS) + '\n')
	for k, v in dataDict.items():
		outFile.write(','.join(str(i) for i in k) + ',' + ','.join(str(v[i]) for i in v) + '\n')