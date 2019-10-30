import os

dataDict = {}

for root, dirs, files in os.walk('./collisionTypeAnalysis', topdown = False):
	for name in files:
		N, _, _, _, compileType, _ = name.strip().split('-')
		N = int(N)
		with open(root + '/' + name, 'r') as fileIn:
			for _ in range(N):
				fileIn.readline()
			totalPColl, totalWColl = 0, 0
			for _ in range(N):
				*_, pColl, wColl = fileIn.readline().strip().split()
				totalPColl += int(pColl)
				totalWColl += int(wColl)
		dataDict[(N, compileType)] = (totalPColl, totalWColl)

with open('./collisionType.csv', 'w') as fileOut:
	fileOut.write('N,type,ppCollisions,pwCollisions\n')
	fileOut.write('\n'.join(','.join(str(i) for i in k) + ',' + ','.join(str(i) for i in v) for k, v in dataDict.items()))
	fileOut.write('\n')