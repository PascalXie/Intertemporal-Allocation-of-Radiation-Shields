import numpy as np
import matplotlib.pyplot as plt
import csv
import string
import seaborn as sns
import math

from scipy.optimize import linprog
import sys
from queue import Queue

delta_ = 0.8
# curve fitting
a_ = [] # pDR
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



def test4():
	""" the result is xxx"""
	c = np.array([0.5,0.1,0.25])
	A = np.array([[5,3,4], [0.6,0.8,1.5]])
	b = np.array([24,4])
	Aeq = None
	beq = None
	bounds = [(0, None), (0, None) , (0, None)]

	solver = ILP(c, A, b, Aeq, beq, bounds)
	solver.solve()

	print("Test 4's result:", solver.opt_val, solver.opt_x)
	return

def IntertemporalAllocation_ILP(N_Total):
	""" the result is xxx"""
	c = []

	# matrix C
	for a in a_:
		a = -1.*a
		c.append(a)
	for a in a_:
		a = -1.*a
		c.append(delta_*a)
	c = np.array(c)
	print('Matrix c',c)

	# matrix A_ub
	A = []
	for a in a_:
		A.append(1.)
	for a in a_:
		A.append(1.)
	A = np.array(A)
	print('Matrix A',A)

	b = np.array([N_Total])
	Aeq = None
	beq = None

	bounds = []
	for a in a_:
		bounds.append((0, 8))
	for a in a_:
		bounds.append((0, 8))
	bounds = np.array(bounds)
	print('Matrix bounds',bounds)

	solver = ILP(c, A, b, Aeq, beq, bounds)
	solver.solve()

	print("Test 4's result:", solver.opt_val, solver.opt_x)

	# get results
	N_H_T = solver.opt_x
	print(N_H_T)

	return N_H_T

def GetDoses(N_H_T):
	print('GetDoses N_H_T ',N_H_T)
	Doses = []

	# time period 1
	for i in range(6):
		ID = i
		NumShield = N_H_T[ID]
		d = dose(NumShield, ID)
		Doses.append(d)

	# time period 2
	for i in range(6):
		ID = i
		NumShield = N_H_T[6+ID]
		d = dose(NumShield, ID)
		Doses.append(d)
	
	print('Doses ',Doses)

	return

if __name__ == '__main__':
	print("hello")

	# step 1 : read fitting results
	path = '../step3-doseAndShields/'
	fileName = 'FittingResults_shieldNum_forStep4.txt'
	ReadFile_dose_fitting(path, fileName)
	for i, a in enumerate(a_):
		D_0shield = D_0shield_[i]
		print(i, a, D_0shield)

	# step 2 : ILP
	#N_Total = 48
	#IntertemporalAllocation_ILP(N_Total)

	for i in range(4):
		N_Total = 12. + float(i)*12.
		print('N_Total ',N_Total)
		N_H_T = IntertemporalAllocation_ILP(N_Total)

		# step 3 : calculate doses
		doses = GetDoses(N_H_T)

		fo = open("IntertemporalAllocationResults_ILP_"+str(int(N_Total))+".txt","w")
		for j in range(6):
			N_h_1 = N_H_T[j]
			N_h_2 = N_H_T[6+j]
			line = str(int(N_Total)) + '  ' + str(int(N_h_1)) + '  ' + str(int(N_h_2)) + '\n'
			fo.write(line)
		fo.close()




