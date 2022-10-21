import numpy as np
import matplotlib.pyplot as plt
import csv
import string
import seaborn as sns
import math
import random
import os 
from scipy.optimize import curve_fit

time = 6. # 6 hours is the time that the mission lasts

# curve fitting
def func(x, a, D_0shield):
	return a * x + D_0shield


def ReadFile(path, fileName):
	print("Data File Path : {}".format(path))
	print("File Name : {}".format(fileName))

	# read
	f = open(path+fileName)
	lines = f.readlines()

	Doses = []
	Dists = []
	ShieldNums = []

	for i in range(6):
		Doses.append([])
		dist = 12. + 2.*float(i)
		Dists.append(dist)

	for i in range(9):
		ShieldNums.append(float(i))


	for line in lines:
		line = line.strip().split()

		name = str(line[0])
		ID			= int(name.split('_')[1])
		ShieldNum	= int(name.split('_')[3])
		edep = float(line[1])/1.e3 # MeV


		mass = 7.45 # kg
		dose = edep * 1.6e-13 / mass * 3600./(8.06e-11) # Gy/h
		Doses[ID].append(dose)

	for i in range(6):
		print('Dist ',Dists)
		print('Doses[i] ',Doses[i])

	# curve fitting
	for i in range(6):
		print('-----')
		print('--------------------------')
		print('Distance ', Dists[i])
		print('Doses[i] ',Doses[i])

		popt, pcov = curve_fit(func, ShieldNums, Doses[i])
		perr = np.sqrt(np.diag(pcov))
		print("popt",popt)
		print("pcov",pcov)
		print("perr",perr)

	return Dists, ShieldNums, Doses

if __name__ == '__main__':
	print("hello")

	fo = open("FittingResults_shieldNum.txt","w")
	fo2 = open("Doses_shieldNum_dist.txt","w")
	fo3 = open("FittingResults_shieldNum_forStep4.txt","w")

	alphas = []

	Dists, ShieldNums, Doses = ReadFile("./","output.txt")

	#
	# plot
	#
	fig = plt.figure(figsize=(7,5))

	# shieldNumber - dose
	for i in range(6):
		# fitting results
		popt, pcov = curve_fit(func, ShieldNums, Doses[i])
		perr = np.sqrt(np.diag(pcov))
		print("popt",popt)
		print("pcov",pcov)
		print("perr",perr)

		line =        '$'+ str(round(popt[0],2)) + '\pm' + str(round(perr[0],2)) + '$    ' 
		line = line + '$'+ str(round(popt[1],2)) + '\pm' + str(round(perr[1],2)) + '$\n'
		fo.write(line)

		line =        str(popt[0]) + ' ' + str(popt[1]) + '\n'
		fo3.write(line)

		for j, dose in enumerate(Doses[i]):
			dist = Dists[i]
			num = ShieldNums[j]
			line2 = 'UGV'+str(i+1) + ' ' + str(dist) + ' ' + str(num) + ' ' + str(dose) + '\n'
			fo2.write(line2)
		
		#alpha = -1.*popt[0]*popt[1]
		alpha = popt[1]
		alphas.append(alpha)

	fo.close()
	fo2.close()
	fo3.close()

	# alpha for intertemporal equilibrium
	alphas = np.array(alphas)
	alphas = alphas/np.sum(alphas)
	print(alphas)

	fo = open("Alphas_IntertemporalEquilibrium.txt","w")
	for i,alpha in enumerate(alphas):
		line = str(i) + " " + str(alpha) + '\n'
		print(line)
		fo.write(line)

	fo.close()




	# shieldNumber - dose
	for i in range(4):
		# simulated results
		plt.plot(ShieldNums, Doses[i],  label = 'UGVID'+str(i), lw=2, ls='-', marker='o')
		# fitting results
		popt, pcov = curve_fit(func, ShieldNums, Doses[i])
		x = np.array(ShieldNums)
		y = popt[0]*x + popt[1]
		plt.plot(x,y,  label = r'UGVID'+str(i)+'-Fit', lw=1, ls='--', marker='s')


	plt.legend(frameon=True)
	plt.xlabel('Number of shields',fontdict={'family' : 'Times New Roman', 'size': 12})
	plt.ylabel('Absorption Dose Rate (Gy/h)',fontdict={'family' : 'Times New Roman', 'size': 12})
	plt.title('Shielding Abilities')
	plt.xlim(-0.5,11)
	#plt.ylim(0,2500)
	plt.tight_layout()
	plt.savefig('figure-level6_shieldingAbilities_shields.png',dpi=300)

	plt.show()
