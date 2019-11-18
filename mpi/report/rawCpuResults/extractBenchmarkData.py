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

_OUTPUT_FILE_DICT = {'SeqBenchmark': 'seqBenchmarkData.csv', 'ParBenchmark': 'parBenchmarkData.csv'}
_PARAM_STRINGS = ('N', 'L', 'r', 'steps', 'threads', 'configType')

_MAP_DICT = {'XeS4114': 'PureXeS', '7700k': 'PureI7'}

outputDict = {i:{} for i in _OUTPUT_FILE_DICT}
dataStruct = re.compile('^.+?(?=#)')
for root, dirs, files in os.walk('.', topdown = False):
	for name in files:
		if len(root.split('/')) < 2 or (root.split('/')[1] not in _OUTPUT_FILE_DICT):
			continue
		# print(root, dirs, name)
		benchmarkType = root.split('/')[1]

		if benchmarkType == 'ParBenchmark':
			N, L, r, steps, threads, _ = name.split('-')
			N = int(N)
			L = float(L)
			r = float(r)
			steps = int(steps)
			threads = int(threads)
			machine = _MAP_DICT[root.split('/')[2]]
			paramTuple = (N, L, r, steps, threads, machine)
		else:
			N, L, r, steps, _ = name.split('-')
			N = int(N)
			L = float(L)
			r = float(r)
			steps = int(steps)
			threads = 1
			machine = _MAP_DICT[root.split('/')[2]]
			paramTuple = (N, L, r, steps, threads, machine)

		dataDict = outputDict[benchmarkType]
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

for benchmarkType, outputFile in _OUTPUT_FILE_DICT.items():
	dataDict = outputDict[benchmarkType]
	with open(outputFile, 'w') as outFile:
		outFile.write(','.join(_PARAM_STRINGS) + ',' + _CLOCK_STRING + ',' + ','.join(_INTERESTED_STATS) + '\n')
		for k, v in dataDict.items():
			outFile.write(','.join(str(i) for i in k) + ',' + ','.join(str(v[i]) for i in v) + '\n')