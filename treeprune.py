import numpy as np
import sys
import matplotlib.pyplot as plt

try:
	fn = sys.argv[1]
except IndexError:
	exit(-1)

N = 48

two_spin = np.zeros((N, N))
three_spin = np.zeros((N, N))
four_spin = np.zeros((N, N))

# read the treepruning output file from Tourbillon.
# this file is arranged 
#   spin_i, spin_j, two_spin_score, three_spin_score, four_spin_score
# we place these in the two_spin/three_spin/four_spin arrays
with open(fn, "r") as f:
	for l in f:
		k = l.strip().split()
		i = int(k[0])
		j = int(k[1])
		tws = float(k[2])
		ths = float(k[3])
		fos = float(k[4])
		two_spin[i, j] = tws
		two_spin[j, i] = tws
		three_spin[i, j] = ths
		three_spin[j, i] = ths
		four_spin[i, j] = fos
		four_spin[j, i] = fos
		

spin_interested = 0

scores = []

#with open("spinpairs.nn", "w") as outnn:
# and then loop over the spins
for spin_interested in range(0, 48):
	#print("Spin %d" % (spin_interested))
	two_spin_ss = two_spin[spin_interested, :]
	three_spin_ss = three_spin[spin_interested, :]
	four_spin_ss = four_spin[spin_interested, :]
	# average them to account for there being many more three spin than two spin, and more four than three.
	sf_three_spin = np.mean(two_spin_ss) / np.mean(three_spin_ss)
	sf_four_spin = np.mean(two_spin_ss) / np.mean(four_spin_ss)
	
	three_spin_ss *= sf_three_spin
	four_spin_ss *= sf_four_spin
	
	# calculate pair score
	total = two_spin_ss + three_spin_ss + four_spin_ss
	#fig, ax = plt.subplots()
	mins = 0.2
	
	# add to score array the spin considered, the neighbour, and the score. 
	for i in range(0, len(total)):
		if (i > spin_interested):
			scores.append([spin_interested, i, total[i]])

scores = np.array(scores)
print(np.shape(scores))
scores = scores[scores[:, 2].argsort()[::-1]]
print(scores[:, 0])

# and now output the nearest spin pairs.
nmax = 16 * N // 2 # number of operators. 
# in Perras, they considered eg 16 nearest neighbours for each site. 
# but each nearest neighbour will count to two spins, hence the //2
with open("spinpairs.nn", "w") as outnn:
	for n in range(0, nmax):
		outnn.write("%d %d\n" % (scores[n,0], scores[n,1]))
	



