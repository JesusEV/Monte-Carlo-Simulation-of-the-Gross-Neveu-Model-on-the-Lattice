import sys
import numpy as np
import time
from matrix import Mmatrix

def function(rows, columns,rand_set):
	M = Mmatrix(rows, columns, rand_set)
	det = np.linalg.det(M)
	return det

def obtain_factors_ids(i, points):
	binary = np.array([int(x) for x in list('{0:0b}'.format(i))])
	zeros = np.zeros(2*points-binary.size)
	return np.append(zeros, binary)

def coef_polinom_calc(points):
	file = open('pol_coefficients.txt', 'a')
	Q_set = np.zeros(points)
	for i in range(2**(2*points)):
		factors_ids = obtain_factors_ids(i, points)
		M = Mmatrix(rows, columns, Q_set)		
		for j, factor in zip(list(range(2*points)),factors_ids):
			if factor == 1:
				M[:, j] = np.zeros(2*points)
				M[j, :] = np.zeros(2*points)
				M[j,j] = 1
		det = np.linalg.det(M)
		if i == 2**(2*points): file.write(str(det))
		else: file.write(str(det) + '\n')
	file.close()

def obtain_coef_from_file():
	file = open('pol_coefficients.txt', 'r')
	content = file.read()
	file.close()
	list_content = content.split('\n')
	list_content.pop() # last line is blanck always
	return np.array([float(x) for x in list_content])
	
	
def obtain_literals(points, rand_set):
	literals = []
	for i in range(2**(2*points)):
		factors_ids = obtain_factors_ids(i, points)
		literal = 1
		for j, factor in zip(list(range(2*points)),factors_ids):
			if factor == 1: literal *= rand_set[j//2]
		literals.append(literal)

	return np.array(literals)

if __name__ == '__main__':

	coef_calc = sys.argv[1]
	rows = 3
	columns = 3
	points = rows*columns

	if coef_calc == '1': coef_polinom_calc(points)

	if coef_calc == '0': 
		rand_set = [np.random.standard_normal(1) for j in range(0,points)]

		start = time.time()
		coef = obtain_coef_from_file()

		i=0
		for c in coef: 
			if c == 0.0: i+=1
		
		print(i)

		"""
		literals = obtain_literals(points, rand_set)
		det_pol = np.dot(coef, literals)
		end = time.time()
		print('\n')
		print('Polinomical determinant calculation: ', det_pol)
		print('time: ', end - start)

		M = Mmatrix(rows, columns, rand_set)
		start = time.time()		
		det_common = np.linalg.det(M)
		end = time.time()
		print('Standard determinant calculation: ', det_common)
		print('time: ', end - start)
		print('\n')
		"""
		

