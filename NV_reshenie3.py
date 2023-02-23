from qutip import Qobj, Options, mesolve
import qutip as q
import numpy as np
import scipy
from math import *
from math import pi, exp
import matplotlib.pyplot as plt
hamiltonian1 = np.zeros((3, 3))
hamiltonian_spin = np.zeros((3, 3))
hamiltonian = np.kron(hamiltonian1, hamiltonian_spin)

H = Qobj(hamiltonian)



def reshenie3(s, index1, index2, j, endtime):

	g0 = (1./(13*10**(-9)))**(0.5)
	gl = (s/(13*10**(-9)))**(0.5)
	g1 = (1./(170*10**(-9)))**(0.5)
	g_prod_relax_3A2 = (2/(4*10**(-3)))**(0.5)
	g_T2 = (2/(1*10**(-6)))**(0.5)
	
	matr_size = 9
	
	def operator(n, m, g):
		A = np.zeros((matr_size, matr_size)) 
		A[n][m] = g
		A = Qobj(A)
		return A
		
	def operator1(g):
		A = np.zeros((matr_size, matr_size)) 
		A[0][0] = g
		A[2][2] = -g
		A = Qobj(A)
		return A
	
	L1 = operator(3, 0, gl)		#3A2 ms=0 -> 3E ms=0
	L2 = operator(4, 1, gl)		#3A2 ms=-1 -> 3E ms=-1	
	L3 = operator(5, 2, gl)		#3A2 ms=+1 -> 3E ms=+1
	L4 = operator(0, 3, g0)		#3E ms=0 -> 3A2 ms=0
	L5 = operator(1, 4, g0)		#3E ms=-1 -> 3A2 ms=-1
	L6 = operator(2, 5, g0)		#3E ms=+1 -> 3A2 ms=+1
	L7 = operator(6, 4, g0)		#3E ms=-1 -> 1E
	L8 = operator(6, 5, g0)		#3E ms=+1 -> 1E
	L9 = operator(0, 6, g1*0.8)	#1E -> 3A2 ms=0
	L91 = operator(1, 6, g1*0.1)	#1E -> 3A2 ms=-1
	L92 = operator(2, 6, g1*0.1)	#1E -> 3A2 ms=+1
	
	#Перекрёстные переходы 3A2 ->3E
	
	L10 = operator(4, 0, gl/100.)		#3A2 ms=0 -> 3E ms=-1
	L11 = operator(5, 0, gl/100.)		#3A2 ms=0 -> 3E ms=+1
	#
	L12 = operator(3, 1, 0)		#3A2 ms=-1 -> 3E ms=0
	L13 = operator(5, 1, 0)		#3A2 ms=-1 -> 3E ms=+1
	L14 = operator(3, 2, 0)		#3A2 ms=+1 -> 3E ms=0
	L15 = operator(4, 2, 0)		#3A2 ms=+1 -> 3E ms=-1
	
	#Перекрёстные переходы 3E -> 3A2
	
	L16 = operator(1, 3, g0/100.)		#3E ms=0 -> 3A2 ms=-1
	L17 = operator(2, 3, g0/100.)		#3E ms=0 -> 3A2 ms=+1
	#
	L18 = operator(0, 4, 0)		#3E ms=-1 -> 3A2 ms=0
	L19 = operator(2, 4, 0)		#3E ms=-1 -> 3A2 ms=+1
	L20 = operator(1, 5, 0)		#3E ms=+1 -> 3A2 ms=-1
	L21 = operator(0, 5, 0)		#3E ms=+1 -> 3A2 ms=0
	
	
	#Продольная релаксация на 3А2
	
	L22 = operator(1, 0, g_prod_relax_3A2)		#ms=0 -> ms=-1
	L23 = operator(2, 0, g_prod_relax_3A2)		#ms=0 -> ms=+1
	L24 = operator(0, 1, g_prod_relax_3A2)		#ms=-1 -> ms=0
	L25 = operator(2, 1, g_prod_relax_3A2)		#ms=-1 -> ms=+1
	L26 = operator(0, 2, g_prod_relax_3A2)		#ms=+1 -> ms=0
	L27 = operator(1, 2, g_prod_relax_3A2)		#ms=+1 -> ms=-1
	
	#Поперечная релаксация на 3А2
	
	L28 = operator1( g_T2) #12, 21
	L29 = operator(0, 2, g_T2) #13, 31
	L30 = operator(1, 2, g_T2) #23, 32
	
	

	options = Options(nsteps=1e6, atol=1e-5)

	def bra(x):
		mas = [[0, 0, 0, 0, 0, 0, 0, 0, 0]]
		mas[0][x] = 1
		mas = Qobj(mas)
		return mas

	def ket(x):
		mas = [[0], [0], [0], [0], [0], [0], [0], [0], [0]]
		mas[x][0] = 1
		mas = Qobj(mas)
		return mas


	starttime = 0
	#endtime = 1e-5
	num_intermediate_state = 100


	state_evaluation_times = np.linspace(
		starttime,
		endtime,
		num_intermediate_state
	)


	#__________________________________________-enter density matrix here!-________________________
	'''
	rho0 = np.zeros((matr_size, matr_size))
	rho0[j][j] = 1.0
	print("Начальное состояние")
	print(rho0)
	rho0=Qobj(rho0)
	'''
	A= [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
	psi_ket = [[1.0/sqrt(2)], [1.0/sqrt(2)], [0]]
	psi_bra = [[1.0/sqrt(2), 1.0/sqrt(2), 0]]
	psi_proektor = np.dot(psi_ket, psi_bra)
	print(psi_proektor)
	rho0= np.kron(A, psi_proektor)
	print("Начальное состояние")
	print(rho0)
	rho0=Qobj(rho0)
	
	

	#______________________________________________________________________________________________

	result = mesolve(
		H,
		rho0,
		state_evaluation_times,
		#[L1, L2, L3, L4, L5, L6, L7, L8, L9, L91, L92, L10, L11, L12, L13, L14, L15, L16, L17, L18, L19, L20, L21, L22, L23, L24, L25, L26, L27, L28],
		[L22, L23, L24, L25, L26, L27, L28],
		[],
		options=options
	)




	solve = bra(index2) * (result.states * ket(index1))

	points_to_plot = []

	for s in solve:
		if s == solve[0]:
			points_to_plot.append(0.5)
		else:
			points_to_plot.append(s.data.data.real[0])
		

	return points_to_plot







