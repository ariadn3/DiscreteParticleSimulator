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

seqData = {}
parData = {}
dataStruct = re.compile('^.+?(?=#)')

for root, dirs, files in os.walk('.', topdown = False):
	for name in files:
		fileExtension = name[-4:]
		if fileExtension != '.seq' and fileExtension != '.par':
			continue
		isSeq = (fileExtension == '.seq')

		if not isSeq:
			N, L, r, steps, threads, _ = name.split('-')
			threads = int(threads)
		else:
			N, L, r, steps, _ = name.split('-')

		N = int(N)
		L = float(L)
		r = float(r)
		steps = int(steps)
		machine = root.split('\\')[2]
		paramTuple = (N, L, r, steps, machine) if isSeq else (N, L, r, steps, threads, machine)

		dictToWrite = seqData if isSeq else parData
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

for filePath, paramString, dataDict in ((_SEQUENTIAL_OUTPUT_FILE, _SEQUENTIAL_PARAM_STRINGS, seqData), (_PARALLEL_OUTPUT_FILE, _PARALLEL_PARAM_STRINGS, parData)):
	with open(filePath, 'w') as outFile:
		outFile.write(','.join(paramString) + ',' + _CLOCK_STRING + ',' + ','.join(_INTERESTED_STATS) + '\n')
		for k, v in dataDict.items():
			outFile.write(','.join(str(i) for i in k) + ',' + ','.join(str(v[i]) for i in v) + '\n')