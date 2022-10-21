import numpy as np
import matplotlib.pyplot as plt
import csv
import string
import seaborn as sns
import math
import random
import os 
from scipy.optimize import curve_fit

def ReadFile(path, fileName):
	print("Data File Path : {}".format(path))
	print("File Name : {}".format(fileName))

	# read
	f = open(path+fileName)
	lines = f.readlines()

	NumShields_t1 = []
	NumShields_t2 = []

	for line in lines:
		line = line.strip().split()
		NumShields_t1.append(float(line[1]))
		NumShields_t2.append(float(line[2]))

	return NumShields_t1, NumShields_t2

if __name__ == '__main__':
	print("hello")

	path = 'level6-intertemporalShiledAllocatin/step4-intertemporalEquilibrium/'

	#UGVID_t1 = np.array([1,2,3,4,5,6])
	#UGVID_t2 = UGVID_t1+6

	UGVID_t1 = []
	UGVID_t2 = []
	UGVID_t1.append('P1h1')
	UGVID_t1.append('P1h2')
	UGVID_t1.append('P1h3')
	UGVID_t1.append('P1h4')
	UGVID_t1.append('P1h5')
	UGVID_t1.append('P1h6')
	UGVID_t2.append('P2h1')
	UGVID_t2.append('P2h2')
	UGVID_t2.append('P2h3')
	UGVID_t2.append('P2h4')
	UGVID_t2.append('P2h5')
	UGVID_t2.append('P2h6')

	#
	# plot
	#
	fig = plt.figure(figsize=(7,5))

	for i in range(4):
		num = 12 + 12*i
		fileName = 'IntertemporalAllocationResults_'+str(num)+'.txt'
		NumShields_t1, NumShields_t2 = ReadFile(path, fileName)

		plt.plot(UGVID_t1, NumShields_t1,  label = r'$P1 N_{Total}$='+str(num), lw=2, ls='-', marker='o')
		plt.plot(UGVID_t2, NumShields_t2,  label = r'$P2 N_{Total}$='+str(num), lw=2, ls='--', marker='s')

	plt.legend(frameon=True)
	plt.xlabel('UGV Index in Two Time Periods',fontdict={'family' : 'Times New Roman', 'size': 12})
	plt.ylabel('Number of Shields',fontdict={'family' : 'Times New Roman', 'size': 12})
	plt.title('Intertemporal Allocation of Shields')
	#plt.xlim(-0.5,11)
	#plt.ylim(0,2500)
	plt.tight_layout()
	plt.savefig('figure-level6_allocationShields.png',dpi=300)
	plt.show()
