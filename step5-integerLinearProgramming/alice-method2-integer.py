import numpy as np
import matplotlib.pyplot as plt
import csv
import string
import seaborn as sns
import math

from scipy.optimize import linprog
import sys
from queue import Queue

class ILP():
    def __init__(self, c, A_ub, b_ub, A_eq, b_eq, bounds):
        """ maximize C^T * X """
        """ s.t. A * X <= B """
		# global variables
        self.LOWER_BOUND = -sys.maxsize
        self.UPPER_BOUND = sys.maxsize
        self.opt_val = None
        self.opt_x = None
        self.Q = Queue()

		# the fixed variables
        self.c = -c
        self.A_eq = A_eq
        self.b_eq = b_eq
        self.bounds = bounds

		# describing the optimized problem
        r = linprog(-c, A_ub, b_ub, A_eq, b_eq, bounds)

		# if the problem is not feasible
        if not r.success:
            raise ValueError('Not a feasible problem!')

		# adding the result and constraints into the list
        self.Q.put((r, A_ub, b_ub))

    def solve(self):
        while not self.Q.empty():
			# get the paramters of the current problem
            res, A_ub, b_ub = self.Q.get(block=False)

			# the value of the optimized cost funtion should be less than the lower bound
            if -res.fun < self.LOWER_BOUND:
                continue

			# if all the elements in x are ingeters
			# or change the lower bound
            if all(list(map(lambda f: f.is_integer(), res.x))):
                if self.LOWER_BOUND < -res.fun:
                    self.LOWER_BOUND = -res.fun

                if self.opt_val is None or self.opt_val < -res.fun:
                    self.opt_val = -res.fun
                    self.opt_x = res.x

                continue

			# branching
            else:
				# find the first non-integer value, record its index, denoted by 'idx'
                idx = 0
                for i, x in enumerate(res.x):
                    if not x.is_integer():
                        break
                    idx += 1

				# construct the new constraints
                new_con1 = np.zeros(A_ub.shape[1])
                new_con1[idx] = -1
                new_con2 = np.zeros(A_ub.shape[1])
                new_con2[idx] = 1
                new_A_ub1 = np.insert(A_ub, A_ub.shape[0], new_con1, axis=0)
                new_A_ub2 = np.insert(A_ub, A_ub.shape[0], new_con2, axis=0)
                new_b_ub1 = np.insert(
                    b_ub, b_ub.shape[0], -math.ceil(res.x[idx]), axis=0)
                new_b_ub2 = np.insert(
                    b_ub, b_ub.shape[0], math.floor(res.x[idx]), axis=0)

				# adding the new constaints into the queue
                r1 = linprog(self.c, new_A_ub1, new_b_ub1, self.A_eq,
                             self.b_eq, self.bounds)
                r2 = linprog(self.c, new_A_ub2, new_b_ub2, self.A_eq,
                             self.b_eq, self.bounds)
                if not r1.success and r2.success:
                    self.Q.put((r2, new_A_ub2, new_b_ub2))
                elif not r2.success and r1.success:
                    self.Q.put((r1, new_A_ub1, new_b_ub1))
                elif r1.success and r2.success:
                    if -r1.fun > -r2.fun:
                        self.Q.put((r1, new_A_ub1, new_b_ub1))
                        self.Q.put((r2, new_A_ub2, new_b_ub2))
                    else:
                        self.Q.put((r2, new_A_ub2, new_b_ub2))
                        self.Q.put((r1, new_A_ub1, new_b_ub1))

def ReadTimes(path, fileName):
	print("Data File Path : {}".format(path))
	print("File Name : {}".format(fileName))

	# read
	f = open(path+fileName)
	lines = f.readlines()

	UGVNames	= []
	times		= []

	for line in lines:
		line = line.strip().split()

		times		.append(float(line[1]))
		UGVNames	.append(str(line[2])) 

	hist = []
	hist.append(UGVNames)
	hist.append(times)

	return hist

def ReadEnvs(path, fileName):
	print("Data File Path : {}".format(path))
	print("File Name : {}".format(fileName))

	# read
	f = open(path+fileName)
	lines = f.readlines()

	UGVName_Env1 = []
	UGVName_Env2 = []

	for line in lines:
		line = line.strip().split()

		ID		= float(line[0])
		Idx_env	= float(line[1])
		UGVName	= str(line[2])
		
		if Idx_env==1:
			UGVName_Env1.append(UGVName)
		elif Idx_env==2:
			UGVName_Env2.append(UGVName)
		else:
			print("!!!!!!!!!!!1")

	#print('UGVName_Env1', len(UGVName_Env1) )
	#print('UGVName_Env2', len(UGVName_Env2) )

	return UGVName_Env1, UGVName_Env2

def matchUGV(UGVTimes, UGVName_Env1, UGVName_Env2):
	UGVInfo = {}
	for i,UGVName in enumerate(UGVTimes[0]):
		time = UGVTimes[1][i]
		info = [0,0]
		if UGVName in UGVName_Env1:
			info = [time, 1]
		elif UGVName in UGVName_Env2:
			info = [time, 2]
		else:
			print("!!!!!")

		UGVInfo[UGVName] = info
			
	#print('UGVInfo', UGVInfo)

	# histogram
	Doses = []

	DoseRate = [5.7, 9.4] # no shielding
	MalfunctionDose = 50.

	for UGVName in UGVInfo:
		time	= UGVInfo[UGVName][0]
		Idx_env = UGVInfo[UGVName][1]

		if Idx_env==1:
			Doses.append(time*DoseRate[0]/MalfunctionDose*100.)
		elif Idx_env==2:
			Doses.append(time*DoseRate[1]/MalfunctionDose*100.)

	hist = np.histogram(Doses, range=(0.,50.))
	print("hist",hist)

	return UGVInfo, hist

def allocation_ILP(Wealth):
	Time_path = '../step1-workingTimes/'
	Env_path = '../step3-workingEnvironments/'

	Time_filename = 'UGVworkingTimes_DZ0.txt'
	Env_filename = 'UGVworkingEnvs_DZ0.txt'
	UGVTimes = ReadTimes(Time_path,Time_filename)
	UGVName_Env1, UGVName_Env2 = ReadEnvs(Env_path, Env_filename)
	UGVInfo_DZ1, hist = matchUGV(UGVTimes, UGVName_Env1, UGVName_Env2)

	Time_filename = 'UGVworkingTimes_DZ1.txt'
	Env_filename = 'UGVworkingEnvs_DZ1.txt'
	UGVTimes = ReadTimes(Time_path,Time_filename)
	UGVName_Env1, UGVName_Env2 = ReadEnvs(Env_path, Env_filename)
	UGVInfo_DZ2, hist = matchUGV(UGVTimes, UGVName_Env1, UGVName_Env2)

	Time_filename = 'UGVworkingTimes_DZ2.txt'
	Env_filename = 'UGVworkingEnvs_DZ2.txt'
	UGVTimes = ReadTimes(Time_path,Time_filename)
	UGVName_Env1, UGVName_Env2 = ReadEnvs(Env_path, Env_filename)
	UGVInfo, hist = matchUGV(UGVTimes, UGVName_Env1, UGVName_Env2)
	UGVInfo_DZ3, hist = matchUGV(UGVTimes, UGVName_Env1, UGVName_Env2)


	C			= [] # vector for the integer linear programming
	P			= [] # price vector for the integer linear programming
	UGVNames	= []

	price = 100.
	time_mission = 6.
	pDRs = [-3.597, -4.4829]
	
#	for UGVName in UGVInfo_DZ1:
#		time	= UGVInfo_DZ1[UGVName][0]
#		Idx_env = UGVInfo_DZ1[UGVName][1]
#		pDR = pDRs[Idx_env-1]
#		#print('UGVName',UGVName,'Idx_env', Idx_env, 'pDR',pDR,'time',time)
#		C.append(pDR*time/time_mission)
#		P.append(price)
#		UGVNames.append(UGVName)
#		
#	for UGVName in UGVInfo_DZ2:
#		time	= UGVInfo_DZ2[UGVName][0]
#		Idx_env = UGVInfo_DZ2[UGVName][1]
#		pDR = pDRs[Idx_env-1]
#		#print('UGVName',UGVName,'Idx_env', Idx_env, 'pDR',pDR,'time',time)
#		C.append(pDR*time/time_mission)
#		P.append(price)
#		UGVNames.append(UGVName)
#
#	for UGVName in UGVInfo_DZ3:
#		time	= UGVInfo_DZ3[UGVName][0]
#		Idx_env = UGVInfo_DZ3[UGVName][1]
#		pDR = pDRs[Idx_env-1]
#		#print('UGVName',UGVName,'Idx_env', Idx_env, 'pDR',pDR,'time',time)
#		C.append(pDR*time/time_mission)
#		P.append(price)
#		UGVNames.append(UGVName)

	#
	# debug
	#
	counter = 0
	for UGVName in UGVInfo_DZ2:
		time	= UGVInfo_DZ2[UGVName][0]
		Idx_env = UGVInfo_DZ2[UGVName][1]
		pDR = pDRs[Idx_env-1]
		#print('UGVName',UGVName,'Idx_env', Idx_env, 'pDR',pDR,'time',time)
		C.append(pDR*time/time_mission)
		P.append(price)
		UGVNames.append(UGVName)

		counter += 1 
		if counter>5:
			break

	c = np.array(C) * -1.
	A = np.array([P])
	b = np.array([Wealth])
	Aeq = None
	beq = None
	bounds = []
	for i in range(len(C)):
		bounds.append((0, 8))
	
	solver = ILP(c, A, b, Aeq, beq, bounds)
	solver.solve()
	
	#print("Allocation result:", solver.opt_val, solver.opt_x)

	N_B = solver.opt_x

	N_B = np.floor(N_B)

	for i in range(N_B.size):
		n = N_B[i]
		if n>=9:
			n = 8
		N_B[i] = n

	allocation = []
	allocation.append(UGVNames)
	allocation.append(N_B)

	return allocation

if __name__ == '__main__':
	print("hello")


	Wealth = 1.e3
	allocation = allocation_ILP(Wealth)
	UGVNames = allocation[0]
	N_B = allocation[1]

	for i,name in enumerate(UGVNames):
		print(name, N_B[i])
