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

_SEQUENTIAL_OUTPUT_FILE = 'seqData.csv'
_SEQUENTIAL_PARAM_STRINGS = ('N', 'L', 'r', 'steps', 'machine')
_PARALLEL_OUTPUT_FILE = 'parData.csv'
_PARALLEL_PARAM_STRINGS = ('N', 'L', 'r', 'steps', 'threads', 'machine')
_PARALLEL_2_OUTPUT_FILE = 'parData2.csv'
_PARALLEL_2_PARAM_STRINGS = ('N', 'L', 'r', 'steps', 'threads', 'gridSize', 'machine')

fileExtensionDict = {'.seq': {}, '.par': {}, '.par2': {}}
dataStruct = re.compile('^.+?(?=#)')

for root, dirs, files in os.walk('.', topdown = False):
	for name in files:
		fileExtension = '.' + name.rpartition('.')[-1]
		if fileExtension not in fileExtensionDict:
			continue

		if fileExtension == '.par2':
			N, L, r, steps, threads, gridSize, _ = name.split('-')
			threads = int(threads)
			gridSize = int(gridSize)
		if fileExtension == '.par':
			N, L, r, steps, threads, _ = name.split('-')
			threads = int(threads)
		elif fileExtension == '.seq':
			N, L, r, steps, _ = name.split('-')

		N = int(N)
		L = float(L)
		r = float(r)
		steps = int(steps)
		machine = root.split('/')[2]

		if fileExtension == '.par2':
			paramTuple = (N, L, r, steps, threads, gridSize, machine)
		elif fileExtension == '.par':
			paramTuple = (N, L, r, steps, threads, machine)
		elif fileExtension == '.seq':
			paramTuple = (N, L, r, steps, machine) 

		dictToWrite = fileExtensionDict[fileExtension]
		with open(root + '/' + name) as fileOut:
			lines = fileOut.readlines()
			wallClockTime = _CLOCK_DTYPE(lines[-2].split()[0])

			if paramTuple not in dictToWrite or wallClockTime < dictToWrite[paramTuple][_CLOCK_STRING]:
				curStats = {_CLOCK_STRING: wallClockTime}
				for line in lines:
					line = line.strip()
					regexSearch = dataStruct.findall(line)
					if regexSearch:
						regexSearch = regexSearch[0].strip().replace(',', '')
						linePart = regexSearch.split()
						if any(stat in linePart for stat in _INTERESTED_STATS):
							curStats[linePart[1]] = _INTERESTED_STATS[linePart[1]](linePart[0])
				dictToWrite[paramTuple] = curStats

outputFiles = ((_SEQUENTIAL_OUTPUT_FILE, _SEQUENTIAL_PARAM_STRINGS, fileExtensionDict['.seq']),
	(_PARALLEL_OUTPUT_FILE, _PARALLEL_PARAM_STRINGS, fileExtensionDict['.par']),
	(_PARALLEL_2_OUTPUT_FILE, _PARALLEL_2_PARAM_STRINGS, fileExtensionDict['.par2']))

for filePath, paramString, dataDict in outputFiles:
	with open(filePath, 'w') as outFile:
		outFile.write(','.join(paramString) + ',' + _CLOCK_STRING + ',' + ','.join(_INTERESTED_STATS) + '\n')
		for k, v in dataDict.items():
			outFile.write(','.join(str(i) for i in k) + ',' + ','.join(str(v[i]) for i in v) + '\n')