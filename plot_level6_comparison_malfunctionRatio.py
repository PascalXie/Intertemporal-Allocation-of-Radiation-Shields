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

	MeanOfMR_All_t1 = []
	MeanOfMR_All_t2 = []
	StdDevOfMR_All_t1 = []
	StdDevOfMR_All_t2 = []

	for line in lines:
		line = line.strip().split()
		MeanOfMR_All_t1.	append(float(line[1]))
		StdDevOfMR_All_t1.	append(float(line[2]))
		MeanOfMR_All_t2.	append(float(line[3]))
		StdDevOfMR_All_t2.	append(float(line[4]))

	return MeanOfMR_All_t1, MeanOfMR_All_t2, StdDevOfMR_All_t1, StdDevOfMR_All_t2


if __name__ == '__main__':
	print("hello")

	TotalNumShields = []
	for i in range(4):
		num = 12 + 12*i
		TotalNumShields.append(num)

	fig = plt.figure(figsize=(7,5))
	ax1 = plt.subplot(111)

	# equilivrium
	path = 'level6-intertemporalShiledAllocatin/step4-intertemporalEquilibrium/'
	fileName = 'MalfunctionRatio.txt'
	MeanOfMR_All_t1, MeanOfMR_All_t2, StdDevOfMR_All_t1, StdDevOfMR_All_t2= ReadFile(path, fileName)

	#x		= TotalNumShields
	#y 		= MeanOfMR_All_t1
	#yerr	= StdDevOfMR_All_t1
	#ax1.errorbar(x,y, yerr=yerr, label='IE P1', lw=1, ls='-', marker='o', elinewidth=3, capsize=4)

	#x		= TotalNumShields
	#y 		= MeanOfMR_All_t2
	#yerr	= StdDevOfMR_All_t2
	#ax1.errorbar(x,y, yerr=yerr, label='IE P2', lw=1, ls='--', marker='s', elinewidth=3, capsize=4)

	plt.plot(TotalNumShields, MeanOfMR_All_t1,   label = 'P1-IE', lw=2, ls='-', marker='o')
	plt.plot(TotalNumShields, MeanOfMR_All_t2,   label = 'P2-IE', lw=2, ls='--', marker='s')

	path = 'level6-intertemporalShiledAllocatin/step5-integerLinearProgramming/'
	fileName = 'MalfunctionRatio_ILP.txt'
	MeanOfMR_All_t1, MeanOfMR_All_t2, StdDevOfMR_All_t1, StdDevOfMR_All_t2 = ReadFile(path, fileName)

	#x		= TotalNumShields
	#y 		= MeanOfMR_All_t1
	#yerr	= StdDevOfMR_All_t1
	#ax1.errorbar(x,y, yerr=yerr, label='ILP P1', lw=1, ls='-', marker='o', elinewidth=3, capsize=4)

	#x		= TotalNumShields
	#y 		= MeanOfMR_All_t2
	#yerr	= StdDevOfMR_All_t2
	#ax1.errorbar(x,y, yerr=yerr, label='ILP P2', lw=1, ls='--', marker='s', elinewidth=3, capsize=4)

	plt.plot(TotalNumShields, MeanOfMR_All_t1,   label = 'P1-ILP', lw=1, ls='-', marker='o')
	plt.plot(TotalNumShields, MeanOfMR_All_t2,   label = 'P2-ILP', lw=1, ls='--', marker='s')

	plt.legend(frameon=True)
	plt.xlabel('Total Number of Shields',fontdict={'family' : 'Times New Roman', 'size': 12})
	plt.ylabel('Expectation of Malfunction Ratio (%)',fontdict={'family' : 'Times New Roman', 'size': 12})
	plt.title('Expectation of Malfunction Ratio of The UGV Unit in Two Periods')
	#plt.xlim(-0.5,11)
	#plt.ylim(0,15)
	plt.tight_layout()
	plt.savefig('figure-level6_Comparison_MalfunctionRatios.png',dpi=300)

	plt.show()
