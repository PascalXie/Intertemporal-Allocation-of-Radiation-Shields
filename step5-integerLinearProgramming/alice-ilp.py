from scipy.optimize import linprog
import numpy as np
import math
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


def test1():
    """ the result is [4,2]"""
    c = np.array([40, 90])
    A = np.array([[9, 7], [7, 20]])
    b = np.array([56, 70])
    Aeq = None
    beq = None
    bounds = [(0, None), (0, None)]

    solver = ILP(c, A, b, Aeq, beq, bounds)
    solver.solve()

    print("Test 1's result:", solver.opt_val, solver.opt_x)
    print("Test 1's true optimal x: [4, 2]\n")


def test2():
    """ the result is [2,4]"""
    c = np.array([3, 13])
    A = np.array([[2, 9], [11, -8]])
    b = np.array([40, 82])
    Aeq = None
    beq = None
    bounds = [(0, None), (0, None)]

    solver = ILP(c, A, b, Aeq, beq, bounds)
    solver.solve()

    print("Test 2's result:", solver.opt_val, solver.opt_x)
    print("Test 2's true optimal x: [2, 4]\n")


def test3():
    """ the result is [4,1]"""
    c = np.array([0.5,0.25])
    A = np.array([[5,4], [0.6,1.5]])
    b = np.array([24,4])
    Aeq = None
    beq = None
    bounds = [(0, None), (0, None)]

    solver = ILP(c, A, b, Aeq, beq, bounds)
    solver.solve()

    print("Test 3's result:", solver.opt_val, solver.opt_x)
    print("Test 3's true optimal x: [4,1]\n")

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

if __name__ == '__main__':
    test4()
