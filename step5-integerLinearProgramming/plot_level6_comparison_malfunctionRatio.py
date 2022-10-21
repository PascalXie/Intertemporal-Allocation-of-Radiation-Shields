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
	IE_MeanOfMR_All_t1, IE_MeanOfMR_All_t2, StdDevOfMR_All_t1, StdDevOfMR_All_t2= ReadFile(path, fileName)


	path = 'level6-intertemporalShiledAllocatin/step5-integerLinearProgramming/'
	fileName = 'MalfunctionRatio_ILP.txt'
	ILP_MeanOfMR_All_t1, ILP_MeanOfMR_All_t2, StdDevOfMR_All_t1, StdDevOfMR_All_t2 = ReadFile(path, fileName)

