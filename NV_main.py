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
endtime = 1e-5
num_intermediate_state = 100
state_evaluation_times = np.linspace(starttime, endtime, num_intermediate_state )


fig, ax = plt.subplots()

#ax.plot(state_evaluation_times, points_to_plot, label=text)
start_naselennost = 0

#ax.plot(state_evaluation_times, reshenie2(1, 0, start_naselennost, endtime), label="3A2: ms = 0")
#ax.plot(state_evaluation_times, reshenie2(1, 1, start_naselennost, endtime), label="3A2: ms = -1")
#ax.plot(state_evaluation_times, reshenie2(1, 2, start_naselennost, endtime), label="3A2: ms = +1")
#ax.plot(state_evaluation_times, reshenie2(1, 3, start_naselennost, endtime), label="3E: ms = 0")
#ax.plot(state_evaluation_times, reshenie2(1, 4, start_naselennost, endtime), label="3E: ms = -1")
#ax.plot(state_evaluation_times, reshenie2(1, 5, start_naselennost, endtime), label="3E: ms = +1")
#ax.plot(state_evaluation_times, reshenie2(1, 6, start_naselennost, endtime), label="1E")


#ax.plot(state_evaluation_times, reshenie2(10, 6, 0, endtime), label="1E")
#ax.plot(state_evaluation_times, reshenie2(10, 0, 0, endtime), label="3A2 ms=0")
#ax.plot(state_evaluation_times, reshenie1(10, 3, 0, endtime), label="ms=0, s=10")
#ax.plot(state_evaluation_times, reshenie1(10, 3, 1, endtime), label="ms=-1, s=10")
#ax.plot(state_evaluation_times, reshenie1(1, 3, 0, endtime), label="ms=0, s=1")
#ax.plot(state_evaluation_times, reshenie1(1, 3, 1, endtime), label="ms=-1, s=1")

#ax.plot(state_evaluation_times, reshenie1(10, 3, 0, endtime), label="ms=0, s=0.1")
#ax.plot(state_evaluation_times, reshenie1(10, 3, 1, endtime), label="ms=-1, s=0.1")

ax.plot(state_evaluation_times, reshenie3(1, 0, 1, start_naselennost, endtime), label="cross-relax for rho_01")
#ax.plot(state_evaluation_times, reshenie3(1, 0, 0, start_naselennost, endtime), label="3A2: ms = 0")
#ax.plot(state_evaluation_times, reshenie3(1, 2, 2, start_naselennost, endtime), label="3A2: ms = +1")

'''
y = np.zeros(num_intermediate_state)
index = 0
for time in state_evaluation_times:
	y[index] = 0.33 + (0.5-0.33)*np.exp(-3*time/(4*1e-3))
	index+=1
ax.plot(state_evaluation_times, y, label = "exp(-3t/T1)")
'''


y = np.zeros(num_intermediate_state)
index = 0
for time in state_evaluation_times:
	y[index] = 0.5*np.exp(-time/(1.0*1e-6))
	index+=1
ax.plot(state_evaluation_times, y, label = "0.5*exp(-t/T2)")

ax.set_xlabel("time, c")
ax.set_ylabel("number")
#ax.set_title("Продольная релаксация(T1 = 4мс, Т2=1мкс )")
ax.set_title("Поперечная релаксация релаксация(T1 = 4мс, Т2=1мкс )")
plt.legend()
plt.show()
exit()




