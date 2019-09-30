from math import sqrt

with open('std-print.in') as fileIn:
	N = int(fileIn.readline())
	L = float(fileIn.readline())
	s = float(fileIn.readline())
	steps = int(fileIn.readline())

maxVelocity = L/4
print('--------- THEORETICAL VELOCITY LIMIT (L/4): {:.14f} ----------\n'.format(maxVelocity))

max1DVelo = -1
max2DVelo = -1
with open('std.out') as fileOut:
	for line in fileOut.readlines():
		step, index, _, _, v_x, v_y, *_ = line.strip().split(' ')
		v_x = float(v_x)
		v_y = float(v_y)
		max1DVelo = max(max1DVelo, v_x, v_y)
		v = sqrt(v_x**2 + v_y**2)
		max2DVelo = max(max2DVelo, v)
		if v > maxVelocity:
			step = int(step)
			index = int(index)
			print('\t Violation: Index {:4d} - Particle {:4d} - Velocity = {:4.8f}'.format(index, index, v))
print()
print('Max 1D velocity: {:4.8f}'.format(max1DVelo))
print('Max 2D velocity: {:4.8f}'.format(max2DVelo))