_WRITE_TO_DIR = 'machineFiles/'

_MACHINES = (1, 2, 4, 8)
_XeS_HOST = False

_NODE_SUFFIX = 'soctf-pdc-'
_XeS_NODES = tuple('{:03d}'.format(i) for i in range(1, 9))
_i7_NODES = tuple('{:03d}'.format(i) for i in range(9, 17))
_XeS_MAX_PROC = 20
_i7_MAX_PROC = 8

for m in _MACHINES:
	if _XeS_HOST:
		with open(_WRITE_TO_DIR + '{}.mf'.format(m*_XeS_MAX_PROC), 'w') as XeSPure:
			XeSPure.write('\n'.join((_NODE_SUFFIX + n) for n in _XeS_NODES[:m]))
			XeSPure.write('\n')
	else:
		with open(_WRITE_TO_DIR + '{}.mf'.format(m*_i7_MAX_PROC), 'w') as i7Pure:
			i7Pure.write('\n'.join((_NODE_SUFFIX + n) for n in _i7_NODES[:m]))
			i7Pure.write('\n')
	with open(_WRITE_TO_DIR + '{}.mf'.format(m*(_i7_MAX_PROC + _XeS_MAX_PROC)), 'w') as mixed:
		mixed.write('\n'.join((_NODE_SUFFIX + n) for pair in (zip(_XeS_NODES, _i7_NODES) if _XeS_HOST else zip(_i7_NODES, _XeS_NODES)) for n in pair))
		mixed.write('\n')
