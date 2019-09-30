from math import sqrt

_FILEPATH = 'slingshot.in'

_N = 4
_L = 12
_r = 1
_STEPS = 20

x = _r + 0.5001
y = _r + 2 + 0.0001
_EPS = 0.1

with open(_FILEPATH, 'w') as slingOut:
    slingOut.write('{}\n{}\n{}\n{}\nprint\n'.format(_N, _L, _r, _STEPS))
    slingOut.write('0 {} {} 3 0\n'.format(x-_EPS, y+2))
    slingOut.write('1 {} {} 0 3\n'.format(x, y-_EPS))
    slingOut.write('2 {} {} -3 0\n'.format(x+6+sqrt(2)+_EPS, y+2-sqrt(2)))
    slingOut.write('3 {} {} 0 3\n'.format(x+6+sqrt(2), y-sqrt(2)-_EPS))
