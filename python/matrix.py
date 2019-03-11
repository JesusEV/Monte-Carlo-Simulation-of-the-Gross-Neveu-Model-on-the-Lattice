import numpy as np

def matrix0_writer(matrix0,N,M, Q_set):
	matrix0[2*N][2*M] = Q_set[N]
	matrix0[2*N][2*M+1] = 0
	matrix0[2*N+1][2*M] = 0
	matrix0[2*N+1][2*M+1] = Q_set[N]

def matrix1_writer(matrix1,N,M, a_p_b_c):
	a = 1
	if a_p_b_c: a = -1
	matrix1[2*N][2*M] = a*1/2
	matrix1[2*N][2*M+1] = a*1/2
	matrix1[2*N+1][2*M] = a*1/2
	matrix1[2*N+1][2*M+1] = a*1/2

def matrix2_writer(matrix2,N,M, a_p_b_c):
	a = 1
	if a_p_b_c: a = -1 
	matrix2[2*N][2*M] = a*1/2
	matrix2[2*N][2*M+1] = a*-1/2
	matrix2[2*N+1][2*M] = a*-1/2
	matrix2[2*N+1][2*M+1] = a*1/2

def matrix3_writer(matrix3,N,M):
	matrix3[2*N][2*M] = 1/2
	matrix3[2*N][2*M+1] = -2
	matrix3[2*N+1][2*M] = 2
	matrix3[2*N+1][2*M+1] = 1/2

def matrix4_writer(matrix4,N,M):
	matrix4[2*N][2*M] = 1/2
	matrix4[2*N][2*M+1] = 2
	matrix4[2*N+1][2*M] = -2
	matrix4[2*N+1][2*M+1] = 1/2

def Mmatrix(rows, columns, Q_set):
	matrix0 = np.zeros((columns*rows*2, columns*rows*2))
	matrix1 = np.zeros((columns*rows*2, columns*rows*2))
	matrix2 = np.zeros((columns*rows*2, columns*rows*2))
	matrix3 = np.zeros((columns*rows*2, columns*rows*2))
	matrix4 = np.zeros((columns*rows*2, columns*rows*2))
	Mmatrix = np.zeros((columns*rows*2, columns*rows*2))

	for n_1 in range(0,columns):
		for n_2 in range(0,rows):
			for m_1 in range(0,columns):
				for m_2 in range(0,rows):

					N = n_2*columns+n_1
					M = m_2*columns+m_1

					#  % operators to take into account boundary conds
					r_0_1 = (n_1)%columns
					r_0_2 = (n_2)%rows

					r_1_1 = (n_1+1)%columns
					r_1_2 = (n_2)%rows

					r_2_1 = (n_1-1)%columns
					r_2_2 = (n_2)%rows

					r_3_1 = (n_1)%columns
					r_3_2 = (n_2+1)%rows

					r_4_1 = (n_1)%columns
					r_4_2 = (n_2-1)%rows

					P0 = r_0_2*columns+r_0_1
					Q0 = m_2*columns+m_1
					
					P1 = r_1_2*columns+r_1_1
					Q1 = m_2*columns+m_1

					P2 = r_2_2*columns+r_2_1
					Q2 = m_2*columns+m_1

					P3 = r_3_2*columns+r_3_1
					Q3 = m_2*columns+m_1

					P4 = r_4_2*columns+r_4_1
					Q4 = m_2*columns+m_1

					
	 
					if P0 == Q0:
						matrix0_writer(matrix0, N, M, Q_set)
					if P1 == Q1:
						a_p_b_c = False
						if n_1 == columns-1: a_p_b_c = True 
						matrix1_writer(matrix1, N, M, a_p_b_c)

					if P2 == Q2:
						a_p_b_c = False
						if n_1 == 0: a_p_b_c = True
						matrix2_writer(matrix2, N, M, a_p_b_c)

					if P3 == Q3:
						matrix3_writer(matrix3, N, M)

					if P4 == Q4:
						matrix4_writer(matrix4, N, M)

					else:
						continue

	Mmatrix = matrix0 - matrix1 - matrix2 - matrix3 - matrix4

	return Mmatrix


if __name__ == '__main__':

	columns = 2
	rows = 2
	Q_set = np.ones(2*columns*rows)
	
	M = Mmatrix(rows, columns, Q_set)
	det = np.linalg.det(M)
	print("Dirac Matrix: ")
	print(M)
	print("Determinant: ", det)