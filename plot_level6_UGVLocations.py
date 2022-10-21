import numpy as np
import matplotlib.pyplot as plt
import csv
import string
import seaborn as sns
import math
import random
import os 
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
fontprops = fm.FontProperties(size=14)

if __name__ == '__main__':
	print("hello")

	#
	# plot
	#
	fig = plt.figure(figsize=(10,4))
	ax1 = plt.subplot(121)
	ax2 = plt.subplot(122)

	Dists = []
	for i in range(6):
		dist = 12. + 2.*float(i)
		Dists.append(dist)
	
	Dists = np.array(Dists)


	Thetas = []
	for i in range(9):
		theta = 10.*float(i)/180.*math.pi
		Thetas.append(theta)

		xs = Dists*math.cos(theta)
		ys = Dists*math.sin(theta)
		ax1.plot(xs, ys, label = r'$N_{h,t}=$'+str(i), lw=2, ls='', marker='s')

	ax1.set_xlim(-14,25)
	ax1.set_ylim(-14,25)
	ax1.set_xlabel('X (m)',fontdict={'family' : 'Times New Roman', 'size': 12})
	ax1.set_ylabel('Y (m)',fontdict={'family' : 'Times New Roman', 'size': 12})
	ax1.set_title('Distribution of Shields')
	ax1.legend(frameon=True)

	ax1.text(20, -10, '(a)', fontsize=18)

	#
	Thetas = np.array(Thetas)
	for i in range(6):
		xs = []
		ys = []
		for j in range(9):
			shieldNum = j
			theta = Thetas[shieldNum]
			x = Dists[i]*math.cos(theta)
			y = Dists[i]*math.sin(theta)
			xs.append(x)
			ys.append(y)
		ax2.plot(xs, ys, label = r'UGV h='+str(i), lw=2, ls='', marker='s')

	ax2.set_xlim(-14,25)
	ax2.set_ylim(-14,25)
	ax2.set_xlabel('X (m)',fontdict={'family' : 'Times New Roman', 'size': 12})
	ax2.set_ylabel('Y (m)',fontdict={'family' : 'Times New Roman', 'size': 12})
	ax2.set_title('Distribution of UGVs')
	ax2.legend(frameon=True)

	ax2.text(20, -10, '(b)', fontsize=18)

	#
	ax1.set_aspect(1)
	ax2.set_aspect(1)
	scalebar1 = AnchoredSizeBar(ax1.transData,
								1, '1 m', 'lower left',
								pad=0.2,
								color='black',
								sep=5,
								frameon=True,
								size_vertical=1,
								fontproperties=fontprops)
	scalebar2 = AnchoredSizeBar(ax2.transData,
								1, '1 m', 'lower left',
								pad=0.2,
								color='black',
								sep=5,
								frameon=True,
								size_vertical=1,
								fontproperties=fontprops)

	ax1.add_artist(scalebar1)
	ax2.add_artist(scalebar2)

	plt.tight_layout()
	plt.savefig('figure-level6-DistribuionOfUGVs.png',dpi=300)

	plt.show()
