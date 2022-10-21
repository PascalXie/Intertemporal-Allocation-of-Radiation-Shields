import numpy as np
import matplotlib.pyplot as plt
import csv
import string
import seaborn as sns
import math

def Plot(path, fileName):
	print("Data File Path : {}".format(path))
	print("File Name : {}".format(fileName))

	# read
	f = open(path+fileName)
	lines = f.readlines()

	binCentors = []
	contents   = []

	for line in lines:
		line = line.strip().split()

		binCentors	.append(float(line[0])) 
		contents	.append(float(line[1]))

	# uniformed intensity
	intensity = np.array(contents)
	intensity = intensity/np.sum(contents)

	hist = []
	hist.append(binCentors)
	hist.append(intensity)

	print(intensity)

	return hist

if __name__ == '__main__':
	print("hello")

	path = '../step1-Geant4Simulation/shieldingBoard-cityAtmosphereRadiation/data/'

	#
	#Act_grd[0] = 71.7E-6; // MBq/m2, 129I
	#Act_grd[1] = 29.0;    // MBq/m2, 134Cs
	#Act_grd[2] = 28.67;   // MBq/m2, 137Cs
	Act_grd_list = []
	Act_grd_list.append(71.7E-6)
	Act_grd_list.append(29.0)
	Act_grd_list.append(28.67)

	Act_grd = np.array(Act_grd_list)
	print(Act_grd)
	Act_grd = Act_grd/np.sum(Act_grd)
	print(Act_grd, np.sum(Act_grd))

	Intensity_total = 0

	fo = open("GammaIntensity.txt","w")

	# 129I
	filename = '129I.txt'
	hist = Plot(path,filename)
	energy = hist[0]
	intensity = hist[1]*Act_grd[0]
	print("energy",energy)
	print("intensity", intensity)
	Intensity_total = Intensity_total + np.sum(intensity)

	for i,e in enumerate(energy):
		inten = round(float(intensity[i]),8)
		print(i,e,inten)
		line = '129I  & '+str(e) + ' & ' + str(inten) + '\\\\ \n'
		fo.write(line)


	# 134Cs 
	filename = '134Cs.txt'
	hist = Plot(path,filename)
	energy = hist[0]
	intensity = hist[1]*Act_grd[1]
	print("energy",energy)
	print("intensity", intensity)
	Intensity_total = Intensity_total + np.sum(intensity)
	for i,e in enumerate(energy):
		inten = round(float(intensity[i]),3)
		print(i,e,inten)
		line = '134Cs  & '+str(e) + ' & ' + str(inten) + '\\\\ \n'
		fo.write(line)

	# 137Cs 
	filename = '137Cs.txt'
	hist = Plot(path,filename)
	energy = hist[0]
	intensity = hist[1]*Act_grd[2]
	print("energy",energy)
	print("intensity", intensity)
	Intensity_total = Intensity_total + np.sum(intensity)
	for i,e in enumerate(energy):
		inten = round(float(intensity[i]),3)
		print(i,e,inten)
		line = '137Cs  & '+str(e) + ' & ' + str(inten) + '\\\\ \n'
		fo.write(line)

	print("Intensity_total", Intensity_total)
	fo.close()
