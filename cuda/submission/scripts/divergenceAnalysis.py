normFile = open('divergenceAnalysis/divergence-cpu.out', 'r')
fmadFile = open('divergenceAnalysis/divergence-fmad.out', 'r')
nofmadFile = open('divergenceAnalysis/divergence-nofmad.out', 'r')

_N = 1000
_s = 1000

fmadArray = [0 for _ in range(_s + 1)]
nofmadArray = [0 for _ in range(_s + 1)]

for curStep in range(_s+1):
    for _ in range(_N):
        normLine = normFile.readline()
        fmadLine = fmadFile.readline()
        nofmadLine = nofmadFile.readline()
        if normLine != fmadLine:
            fmadArray[curStep] += 1
        if normLine != nofmadLine:
            nofmadArray[curStep] += 1

normFile.close()
fmadFile.close()
nofmadFile.close()

with open('divergedParticles.csv', 'w') as fileOut:
    fileOut.write('step,fmad,nofmad\n')
    for curStep in range(_s+1):
        fileOut.write('{},{},{}\n'.format(curStep, fmadArray[curStep], nofmadArray[curStep]))
