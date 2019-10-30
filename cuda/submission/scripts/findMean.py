_SIZE = (1, 2, 4, 6, 8, 12, 16, 24, 32)

dataDict = {}
for s in _SIZE:
    with open('filterCollisionsAnalysis/N{}.out'.format(s), 'r') as fileIn:
        lines = fileIn.readlines()
        totalTime = sum(float(l.strip().split()[-1]) for l in lines)
        size = len(lines)
        dataDict[s*1000] = totalTime/size

with open('filterTimes.csv', 'w') as fileOut:
    fileOut.write('N,avgTime\n')
    for k, v in dataDict.items():
        fileOut.write('{},{}\n'.format(k, v))
