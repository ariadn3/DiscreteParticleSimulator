import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from math import pi
from time import time

filepath = input('Enter the input/output file you would like to visualize (xxx.in & xxx.out): ')
# filepath = 'testAnim'

startTime = time()

for extension in ('.in', '.out'):
	try:
		open(filepath + extension).close()
	except FileNotFoundError:
		print('{}{} not found, exiting...'.format(filepath, extension))
		exit(0)

dataDict = {}

with open(filepath + '.in', 'r') as particleIn:
	N = int(particleIn.readline().strip())
	L = float(particleIn.readline().strip())
	r = float(particleIn.readline().strip())
	steps = int(particleIn.readline().strip())

with open(filepath + '.out', 'r') as particleOut:
	for line in particleOut.readlines():
		data = line.strip().split()
		step, index, x, y, _, _ = data
		step = int(step)
		index = int(index)
		x = float(x)
		y = float(y)

		if step not in dataDict:
			dataDict[step] = {}
		dataDict[step][index] = (x, y)

for i in range(steps+1):
	if i not in dataDict:
		print('WARNING: Timestep {} not in output'.format(i))
	if i in dataDict and len(dataDict[i]) != N:
		print('WARNING: Unequal number of particles (detected at timestep {})'.format(i))

fig, ax = plt.subplots(figsize=(8, 8))

def init():
	ax.set_xlim(0, L)
	ax.set_ylim(0, L)
	return []

circlePatch = []
def plotPoints(step):
	circles = []
	for patch in circlePatch:
		patch.remove()
	circlePatch.clear()
	for v in dataDict[step].values():
		newPatch = plt.Circle(v, r)
		circlePatch.append(newPatch)
		circles.append(ax.add_patch(newPatch))
	# scatter.set_offsets(list(dataDict[step].values()))
	return circles

ani = FuncAnimation(fig, plotPoints, frames = list(range(0, steps+1)),
                    init_func = init, blit = True, interval = 10)
ani.save('{}.gif'.format(filepath), writer = 'imagemagick', fps=200)

endTime = time()

print('Written to {}.gif!\n----------'.format(filepath))
print('Time taken: {:.3f}s'.format(endTime-startTime))