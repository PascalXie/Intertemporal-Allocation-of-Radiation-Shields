import numpy as np
import matplotlib.pyplot as plt
import csv
import string
import seaborn as sns
import math
import random
import os 
from scipy.optimize import curve_fit

TotalNumUGVs = 6
delta = 0.8

def ReadFile(path, fileName):
	print("Data File Path : {}".format(path))
	print("File Name : {}".format(fileName))

	# read
	f = open(path+fileName)
	lines = f.readlines()

	Alphas = []

	for line in lines:
		line = line.strip().split()

		alpha = float(line[1])
		Alphas.append(alpha)

	return Alphas

def Equilibrium(TotalNumberShields, Alphas):
	# scenario 2
	
	NumShields_time1 = []
	NumShields_time2 = []

	for i in range(TotalNumUGVs):
		alpha = Alphas[i]
		NumShields_time1.append(0)

		NumShields_A_UGV = alpha*TotalNumberShields
		NumShields_time2.append(NumShields_A_UGV)

	return NumShields_time1, NumShields_time2

def IntertemporalEquilibrium(TotalNumberShields, Alphas):
	# scenario 2
	
	NumShields_time1 = []
	NumShields_time2 = []

	for i in range(TotalNumUGVs):
		alpha = Alphas[i]
		NumShields_A_UGV = alpha/(1.+delta)*TotalNumberShields
		#NumShields_A_UGV = math.floor(NumShields_A_UGV)
		NumShields_A_UGV = round(NumShields_A_UGV,0)
		NumShields_time1.append(NumShields_A_UGV)

		NumShields_A_UGV = delta*alpha/(1.+delta)*TotalNumberShields
		#NumShields_A_UGV = math.floor(NumShields_A_UGV)
		NumShields_A_UGV = round(NumShields_A_UGV,0)
		NumShields_time2.append(NumShields_A_UGV)

	for i,ele1 in enumerate(NumShields_time1):
		if ele1>8:
			NumShields_time1[i] = 8

		ele2 = NumShields_time2[i]
		if ele2>8:
			NumShields_time2[i] = 8

		

	return np.array(NumShields_time1), np.array(NumShields_time2)

if __name__ == '__main__':
	print("hello")

	Alphas = ReadFile("../step3-doseAndShields/","Alphas_IntertemporalEquilibrium.txt")

	#NumShields_time1, NumShields_time2 = IntertemporalEquilibrium(20,Alphas)
	#print(NumShields_time1)
	#print(NumShields_time2)

	for i in range(4):
		num = 12. + float(i)*12
		NumShields_time1, NumShields_time2 = IntertemporalEquilibrium(num,Alphas)
		print('I: ', i, ', Number of shields: ', num)
		print('sum', np.sum(NumShields_time1), NumShields_time1)
		print('sum', np.sum(NumShields_time2), NumShields_time2)

		fo = open("IntertemporalAllocationResults_"+str(int(num))+".txt","w")
		for j,ele1 in enumerate(NumShields_time1):
			ele2 = NumShields_time2[j]
			line = str(int(num)) + '  ' + str(ele1) + '  ' + str(ele2) + '\n'
			fo.write(line)
		fo.close()
