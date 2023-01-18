import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

##########################################################################################
############################## Depleted Boltzmann ########################################
##########################################################################################

def equations(u0,F0):
	u,F = u0,F0

	f = []
	for i in range(len(u)):
		f.append(u[i]-np.exp((F[i]-2.*gamma*N*u[i])/T)/np.sum(np.exp((F-2.*gamma*N*u)/T)))	#Depleted Boltzmann Distribution for E=F-gammaNu
		#f.append(u[i]-np.exp((F[i]-1.*gamma*N*u[i])/T)/np.sum(np.exp((F-1.*gamma*N*u)/T))) #Boltzmann Distribution for E=F-gammaNu
	return f

##########################################################################################
######################### Parameters and Inicialization ##################################
##########################################################################################

num_patches = 100
N = 20.
T = .5
gamma = 1.

x,F0,u0 = [],[],[]
for i in range(num_patches):
	F0.append(float(i)/50.+0)
	x.append(i)

F0 = np.array(F0)
for i in range(num_patches):
	u0.append(np.exp(F0[i]/T)/np.sum(np.exp(F0/T)))

u_1 =  fsolve(equations, u0,F0)

##########################################################################################
##################################### Plot ###############################################
##########################################################################################

plt.plot(x,u_1,c="g",linewidth=2.,alpha=1.,label=r"$\beta =$"+str(T))

plt.legend(loc=2)
plt.xlabel("Patch quality")
plt.ylabel("Density foragers")
plt.savefig("Depleted_Boltzmann_energies.png",dpi=300)

plt.show()
