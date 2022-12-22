import numpy as np
import sys

try:
	ifn = sys.argv[1]
	ofn = sys.argv[2]
except IndexError:
	exit(-1)

data = np.loadtxt(ifn)
data[:, :-1] *= (np.pi / 180)
np.savetxt(ofn, data, delimiter=" ", fmt="%f")