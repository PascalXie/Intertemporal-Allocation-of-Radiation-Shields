import numpy as np
import matplotlib.pyplot as plt
import csv
import string
import seaborn as sns
import math
import random
import os 
from scipy.optimize import curve_fit

# curve fitting
a_ = []
D_0shield_ = []
def dose(x, UGVID):
	return a_[UGVID] * x + D_0shield_[UGVID]

def ReadFile_dose_fitting(path, fileName):
	print("Data File Path : {}".format(path))
	print("File Name : {}".format(fileName))

	# read
	f = open(path+fileName)
	lines = f.readlines()

	for line in lines:
		line = line.strip().split()

		a_.			append(float(line[0]))
		D_0shield_.	append(float(line[1]))

	return 

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

	path = '../step3-doseAndShields/'
	fileName = 'FittingResults_shieldNum_forStep4.txt'
	ReadFile_dose_fitting(path, fileName)
	for i, a in enumerate(a_):
		D_0shield = D_0shield_[i]
		print(i, a, D_0shield)

	#
	# plot
	#
	fig = plt.figure(figsize=(7,5))

	#
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

	# for plot2
	NumShields_Delta_All = []
	# for plot2 !

	# for plot3
	Doses_All_t1 = []
	Doses_All_t2 = []
	# for plot3 !

	#
	time_A_period = 3. # hour
	path = './'
	for i in range(4):
		num = 12 + 12*i
		fileName = 'IntertemporalAllocationResults_'+str(num)+'.txt'
		NumShields_t1, NumShields_t2 = ReadFile(path, fileName)

		# for plot2
		NumShields_Delta = np.array(NumShields_t1)-np.array(NumShields_t2)
		NumShields_Delta_All.append(NumShields_Delta)
		print('NumShields_Delta', NumShields_Delta)
		# for plot2 !

		Doses_t1 = []
		Doses_t2 = [] 
		for j, num_t1 in enumerate(NumShields_t1):
			UGVID = j
			num_t2 = NumShields_t2[j]
			dose_t1 = dose(num_t1,UGVID)*time_A_period
			dose_t2 = dose(num_t2,UGVID)*time_A_period
			print(j, dose_t1, dose_t2, dose_t2-dose_t1)
			Doses_t1.append(dose_t1)
			Doses_t2.append(dose_t2)

		# for plot3
		Doses_All_t1.append(np.array(Doses_t1))
		Doses_All_t2.append(np.array(Doses_t2))
		# for plot3 !

		#plt.plot(UGVID_t1, Doses_t1,  label = r'$P1 N_{Total}$='+str(num), lw=2, ls='-', marker='o')
		#plt.plot(UGVID_t2, Doses_t2,  label = r'$P2 N_{Total}$='+str(num), lw=2, ls='--', marker='s')

#	plt.legend(frameon=True)
#	plt.xlabel('UGV Index in Two Time Periods',fontdict={'family' : 'Times New Roman', 'size': 12})
#	plt.ylabel('Absorption Dose Rate (Gy/hour)',fontdict={'family' : 'Times New Roman', 'size': 12})
#	plt.title('Absorption Doses in Two Periods')
#	#plt.xlim(-0.5,11)
#	#plt.ylim(0,2500)
#	plt.tight_layout()
#	plt.savefig('figure-level6_dosesInTwoPeriods.png',dpi=300)

	#
	# plot 2
	#
	fo = open("Differences_shieldAllocation.txt","w")

	fig = plt.figure(figsize=(7,5))

	TotalNumShields = []
	TotalNumShield_Delta = []
	for i in range(4):
		num = 12 + 12*i
		TotalNumShields.append(num)
		print("Delta,",NumShields_Delta_All[i])
		TotalNumShield_Delta.append(np.sum(NumShields_Delta_All[i])/float(num)*100.)

		line = str(num) + ' '
		for j, delta in enumerate(NumShields_Delta_All[i]):
			line = line + str(int(delta)) + ' '
		line = line + '\n'
		fo.write(line)
	fo.close()

	print(TotalNumShield_Delta)
#	plt.plot(TotalNumShields, TotalNumShield_Delta, lw=2, ls='-', marker='o')
#	#plt.legend(frameon=False)
#	plt.xlabel('Total Number of Shields',fontdict={'family' : 'Times New Roman', 'size': 12})
#	plt.ylabel(r'$\frac{\Delta N}{N_{Total}}$ (%)',fontdict={'family' : 'Times New Roman', 'size': 12})
#	plt.title('Difference of Allocated Shielded in Two Periods')
#	#plt.xlim(-0.5,11)
#	#plt.ylim(0,2500)
#	plt.tight_layout()
#	plt.savefig('figure-level6_NumberShieldsDelta.png',dpi=300)
#
#	#plt.show()

	#
	# plot 3
	#
	MeanOfMR_All_t1 = []
	MeanOfMR_All_t2 = []

	StdDevOfMR_All_t1 = []
	StdDevOfMR_All_t2 = []

	TotalNumShields = []

	for i in range(4):
		num = 12 + 12*i
		TotalNumShields.append(num)
		MeanOfMR_All_t1.append(np.mean(Doses_All_t1[i])/160.*100.)
		MeanOfMR_All_t2.append(np.mean(Doses_All_t2[i])/160.*100.)

		StdDevOfMR_All_t1.append(np.std(Doses_All_t1[i]/160.*100.))
		StdDevOfMR_All_t2.append(np.std(Doses_All_t2[i]/160.*100.))

	print('MeanOfMR_All_t1')
	print(MeanOfMR_All_t1)
	print('MeanOfMR_All_t2')
	print(MeanOfMR_All_t2)

	fo = open("MalfunctionRatio.txt","w")
	for i in range(4):
		num = TotalNumShields[i]
		line = str(num) + ' ' + str(MeanOfMR_All_t1[i]) + ' ' + str(StdDevOfMR_All_t1[i]) + ' ' + str(MeanOfMR_All_t2[i]) + ' ' + str(StdDevOfMR_All_t2[i]) + '\n'
		fo.write(line)
	fo.close()

	#
	# plot 4
	#
	fo = open("DevicesWorkingTime.txt","w")
	Devices_name = np.array(['CPUBoard','Driver(Copley)','BrushedMotor', 'BrushlessMotor', 'CCDCamera(AV,reflective)'])
	Devices_MalDose = np.array([164.07, 217.8, 5089.9, 759.2, 745.7])


	WorkingTimes_All = []

	for i in range(4):
		num = 12 + 12*i
		TotalNumShields.append(num)

		for j, dose_UGV in enumerate(Doses_All_t1[i]):
			# write data
			if i!=3:
				continue

			print('----')
			print('NumID ',i,j, dose_UGV)
			WorkingTimes = Devices_MalDose/dose_UGV*time_A_period
			print(WorkingTimes)

			print(i)
			line = str(j+1) + ' & '

			for k,name in enumerate(Devices_name):
				WorkingTime = round(WorkingTimes[k],1)
				line = line + str(WorkingTime) + ' & '

			dose_UGV_t2 = Doses_All_t2[i][j]
			WorkingTimes = Devices_MalDose/dose_UGV_t2*time_A_period

			for k,name in enumerate(Devices_name):
				WorkingTime = round(WorkingTimes[k],1)
				line = line + str(WorkingTime) + ' & '

			line += '\\\\ \n'
			fo.write(line)
			# write data !

	fo.close()

