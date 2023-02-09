from qutip import Qobj, Options, mesolve
import qutip as q
import numpy as np
import scipy
from math import *
from math import pi, exp
import matplotlib.pyplot as plt
from NV_reshenie1 import reshenie1
from NV_reshenie2 import reshenie2
from NV_reshenie3 import reshenie3
hamiltonian1 = np.zeros((3, 3))
hamiltonian_spin = np.zeros((3, 3))
hamiltonian = np.kron(hamiltonian1, hamiltonian_spin)

H = Qobj(hamiltonian)




starttime = 0
endtime = 1e-6
num_intermediate_state = 100
state_evaluation_times = np.linspace(starttime, endtime, num_intermediate_state )


fig, ax = plt.subplots()

#ax.plot(state_evaluation_times, points_to_plot, label=text)
start_naselennost = 0

#ax.plot(state_evaluation_times, reshenie2(1, 2, start_naselennost, endtime), label="3A2: ms = 0")
#ax.plot(state_evaluation_times, reshenie2(1, 1, start_naselennost, endtime), label="3A2: ms = -1")
#ax.plot(state_evaluation_times, reshenie2(1, 2, start_naselennost, endtime), label="3A2: ms = +1")
#ax.plot(state_evaluation_times, reshenie2(10, 3, start_naselennost, endtime), label="3E: ms = 0")
#ax.plot(state_evaluation_times, reshenie2(10, 4, start_naselennost, endtime), label="3E: ms = -1")
#ax.plot(state_evaluation_times, reshenie2(10, 5, start_naselennost, endtime), label="3E: ms = +1")
#ax.plot(state_evaluation_times, reshenie2(10, 6, start_naselennost, endtime), label="1E")


#ax.plot(state_evaluation_times, reshenie2(10, 6, 0, endtime), label="1E")
#ax.plot(state_evaluation_times, reshenie2(10, 0, 0, endtime), label="3A2 ms=0")
#ax.plot(state_evaluation_times, reshenie1(10, 3, 0, endtime), label="ms=0, s=10")
#ax.plot(state_evaluation_times, reshenie1(10, 3, 1, endtime), label="ms=-1, s=10")
#ax.plot(state_evaluation_times, reshenie1(1, 3, 0, endtime), label="ms=0, s=1")
#ax.plot(state_evaluation_times, reshenie1(1, 3, 1, endtime), label="ms=-1, s=1")

#ax.plot(state_evaluation_times, reshenie1(10, 3, 0, endtime), label="ms=0, s=0.1")
#ax.plot(state_evaluation_times, reshenie1(10, 3, 1, endtime), label="ms=-1, s=0.1")

ax.plot(state_evaluation_times, reshenie3(1, 0, 0, start_naselennost, endtime))

'''
y = np.zeros(num_intermediate_state)
index = 0
for time in state_evaluation_times:
	y[index] = 0.33 + 0.66*np.exp(-time/(0.3333*4*1e-3))
	index+=1
ax.plot(state_evaluation_times, y, label = "exp(-3t/T1)")
'''
'''
s_value = np.linspace(1e-5, 1, 100)
emmiting = np.zeros(100)
index = 0
for s in s_value:
	emmiting[index] = reshenie2(s, 0, 0, endtime)
	index += 1
	print(index)

ax.plot(s_value, emmiting)	
'''
ax.set_xlabel("time, c")
ax.set_ylabel("number")
#ax.set_title("")
plt.legend()
plt.show()
exit()
