import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

plt.rcParams.update({'font.size': 14})

##########################################################################################
##################################### OPTIMAL ############################################
##########################################################################################


# def S_function(rho,i,gamma):
# 	#S = A[i] - gamma*rho[i] - beta*np.sqrt((i-x_nest)**2)										#PROPORTIONAL DEPLETION
# 	#S = A[i] - 2.*gamma*rho[i]																	#PROPORTIONAL DEPLETION Depleted optimization
# 	S = A[i] - A[i]*gamma*rho[i]/(A[i]+gamma*rho[i]+eps) - (A[i]**2)*gamma*rho[i]/(A[i]+gamma*rho[i]+eps)**2
# 	#S = A[i] - np.minimum(A[i],(A[i]-gamma*rho[i])*gamma*rho[i]) - beta*np.sqrt((i-x_nest)**2)
# 	#S = A[i] - A[i]*gamma*rho[i]/(A[i]+gamma*rho[i]+0.00001) - beta*np.sqrt((i-x_nest)**2)		#SATURATING-COMPETITION DEPLETION
# 	#S = A[i] - A[i]*gamma*rho[i]/(gamma*rho[i]+1.) - beta*np.sqrt((i-x_nest)**2)				#ERIC-IFD DEPLETION
# 	#S = A[i] - A[i]*gamma*rho[i]/(A[i]+1.) - beta*np.sqrt((i-x_nest)**2)						#SATURATING-REBECCA DEPLETION
# 	return S
# 
# num_patches = 100 
# N = 20
# gamma = 1.
# beta = 0.
# x_nest = int(num_patches/2)
# eps = 0.000000000000001
# 
# t_total = 10000							#t_total should make the S at the final stage pretty uniform! Otherwise, increase t_total!
# dt = float(N)/float(t_total)
# x = np.linspace(0,num_patches,num=num_patches)
# 
# A,rho = [],[]
# for i in range(num_patches):
# 	A.append(float(i)/50.)
# 	rho.append(0.)
# 	
# print "Food",np.sum(A),"Consumption",N*gamma
# if np.sum(A)<N*gamma:
# 	print np.sum(A),N*gamma
# 	print "ERROR!!!! too much consumption"
# ###################################### ALGORITHM #########################################	
# 
# for t in range(t_total):
# 	S = []
# 	for i in range(num_patches):
# 		S.append(S_function(rho,i,gamma))
# 	
# 	best_x = np.argmax(S)
# 	rho[best_x] += dt
	    
##########################################################################################
####################################### Boltzmann ########################################
##########################################################################################

def equations(u0,F0):
	u,F = u0,F0

	f = []
	for i in range(len(u)):
		f.append(u[i]-np.exp((F[i]-1.*gamma*N*u[i])/T)/np.sum(np.exp((F-1.*gamma*N*u)/T)))
		#f.append(u[i]-np.exp((F[i]-F[i]*gamma*N*u[i]/(F[i]+gamma*N*u[i]+eps))/T)/np.sum(np.exp((F-F*gamma*N*u/(F+gamma*N*u+eps))/T)))
		
		#f.append(u[i]-np.exp(np.minimum(F[i],gamma*N*u[i])/T)/np.sum(np.exp(np.minimum(F,gamma*N*u)/T)))
		#f.append(u[i]-np.exp(np.minimum(F[i],F[i]*gamma*N*u[i])/T)/np.sum(np.exp(np.minimum(F,F*gamma*N*u)/T)))
		#f.append(u[i]-np.exp(np.minimum(F[i],(F[i]-gamma*N*u[i])*gamma*N*u[i])/T)/np.sum(np.exp(np.minimum(F,(F-gamma*N*u)*gamma*N*u)/T)))
		
		#f.append( u[i]-np.exp(( F[i] - F[i]*gamma*N*u[i]/(F[i]+gamma*N*u[i]) )/T)/np.sum(np.exp(( F - F*gamma*N*u/(F+gamma*N*u) )/T)) )

	return f

##################################### Parameters #########################################

num_patches = 100
N = 20.
TT = [5.,2.,.5,0.2,.1]
gamma = 1.
beta = 0.
T,dt = 0.05,0.05

eps = 0.000000000000001

x,F0,u0 = [],[],[]
for i in range(num_patches):
	F0.append(float(i)/50.+0)
	x.append(i)

F0 = np.array(F0)
for i in range(num_patches):
	u0.append(np.exp(F0[i]/TT[-1])/np.sum(np.exp(F0/TT[-1])))
	#u0.append(rho[i]/N)
	#u0.append(1./num_patches)
#print u0

total_food,total_consum = np.sum(F0),N*gamma
if total_food<total_consum:
	print "RESOURCE FINISH!!!"
print total_food,total_consum

# u_beta,energy,entropy,optimality = [],[],[],[]
# for i in range(174):
# 	if i<101:
# 		T += dt
# 	elif i<141:
# 		T += dt*10
# 	elif i<156:
# 		T += dt*100
# 	else:
# 		T += dt*1000
# 	u_beta.append(fsolve(equations, u0,F0))
# 	#energy.append(sum(u_beta[i]*(F0-gamma*N*u_beta[i])))
# 	energy.append(sum(u_beta[i]*(F0-F0*gamma*N*u_beta[i]/(F0+gamma*N*u_beta[i]))))
# 	entropy.append(-sum(u_beta[i]*np.log(u_beta[i])))
# 	optimality.append(1./T)
# 
# print "OPTIMALITY:",optimality
# print "BOLTZMANN ENERGY:",energy
# print "BOLTZMANN ENTROPY:",entropy
	
T=TT[0]
u_1 =  fsolve(equations, u0,F0)
T=TT[1]
u_2 =  fsolve(equations, u0,F0)
T=TT[2]
u_3 =  fsolve(equations, u0,F0)
T=TT[3]
u_4 =  fsolve(equations, u0,F0)
T=TT[4]
u_5 =  fsolve(equations, u0,F0)

##########################################################################################
############################## Modified Boltzmann ########################################
##########################################################################################

def equations(u0,F0):
	u,F = u0,F0

	f = []
	for i in range(len(u)):
		f.append(u[i]-np.exp((F[i]-2.*gamma*N*u[i])/T)/np.sum(np.exp((F-2.*gamma*N*u)/T)))
		#f.append(u[i]-np.exp((F[i]-(F[i]*gamma*N*u[i]/(F[i]+gamma*N*u[i]+eps))*(1.+F[i]/(F[i]+gamma*N*u[i]+eps)))/T)/np.sum(np.exp((F-(F*gamma*N*u/(F+gamma*N*u+eps))*(1.+F/(F+gamma*N*u+eps)))/T)))
		
		#f.append(u[i]-np.exp((np.minimum(F[i],gamma*N*u[i])+gamma*N*u[i])/T)/np.sum(np.exp((np.minimum(F,gamma*N*u)+gamma*N*u)/T)))
		#f.append(u[i]-np.exp((np.minimum(F[i],F[i]*gamma*N*u[i])+gamma*N*u[i])/T)/np.sum(np.exp((np.minimum(F,F*gamma*N*u)+gamma*N*u)/T)))
		#f.append(u[i]-np.exp((np.minimum(F[i],(F[i]-gamma*N*u[i])*gamma*N*u[i])+gamma*N*u[i])/T)/np.sum(np.exp((np.minimum(F,(F-gamma*N*u)*gamma*N*u)+gamma*N*u)/T)))
		
		#f.append( u[i]-np.exp(( F[i] - F[i]*gamma*N*u[i]/(F[i]+gamma*N*u[i]) - (F[i]**2)*gamma*N*u[i]/(F[i]+gamma*N*u[i])**2 )/T)/np.sum(np.exp(( F - F*gamma*N*u/(F+gamma*N*u) - (F**2)*gamma*N*u/(F+gamma*N*u)**2 )/T)) )
	return f

##################################### Parameters #########################################

num_patches = 100
#N = 20.
#T = .1
#gamma = 1.
beta = 0.
T,dt = 0.05,0.05

#eps = 0.00001

x,F0,u0 = [],[],[]
for i in range(num_patches):
	F0.append(float(i)/50.+0)
	x.append(i)

F0 = np.array(F0)
for i in range(num_patches):
	u0.append(np.exp(F0[i]/TT[-1])/np.sum(np.exp(F0/TT[-1])))
	#u0.append(rho[i]/N)
	#u0.append(1./num_patches)
#print u0

total_food,total_consum = np.sum(F0),N*gamma
if total_food<total_consum:
	print "RESOURCE FINISH!!!"
print total_food,total_consum

# u_beta,energy,entropy,optimality = [],[],[],[]
# for i in range(174):
# 	if i<101:
# 		T += dt
# 	elif i<141:
# 		T += dt*10
# 	elif i<156:
# 		T += dt*100
# 	else:
# 		T += dt*1000
# 	u_beta.append(fsolve(equations, u0,F0))
# 	#energy.append(sum(u_beta[i]*(F0-gamma*N*u_beta[i])))
# 	energy.append(sum(u_beta[i]*(F0-F0*gamma*N*u_beta[i]/(F0+gamma*N*u_beta[i]))))
# 	entropy.append(-sum(u_beta[i]*np.log(u_beta[i])))
# 	optimality.append(1./T)
# 
# print "OPTIMALITY:",optimality
# print "MODY BOLTZMANN ENERGY:",energy
# print "MODY BOLTZMANN ENTROPY:",entropy

T=TT[0]
u_1_m =  fsolve(equations, u0,F0)
T=TT[1]
u_2_m =  fsolve(equations, u0,F0)
T=TT[2]
u_3_m =  fsolve(equations, u0,F0)
T=TT[3]
u_4_m =  fsolve(equations, u0,F0)
T=TT[4]
u_5_m =  fsolve(equations, u0,F0)



###################################### IFD ###############################################
IFD =  []
for i in range(len(F0)):
	IFD.append(F0[i]/total_food)
	
##########################################################################################
##################################### OPTIMAL ############################################
##########################################################################################


def S_function(rho,i,gamma):
	S = A[i] - gamma*rho[i] - beta*np.sqrt((i-x_nest)**2)										#PROPORTIONAL DEPLETION
	#S = A[i] - 2.*gamma*rho[i]																	#PROPORTIONAL DEPLETION Depleted optimization
	#S = A[i] - A[i]*gamma*rho[i]/(A[i]+gamma*rho[i]+eps) - (A[i]**2)*gamma*rho[i]/(A[i]+gamma*rho[i]+eps)**2
	#S = A[i] - np.minimum(A[i],(A[i]-gamma*rho[i])*gamma*rho[i]) - beta*np.sqrt((i-x_nest)**2)
	#S = A[i] - A[i]*gamma*rho[i]/(A[i]+gamma*rho[i]+0.00001) - beta*np.sqrt((i-x_nest)**2)		#SATURATING-COMPETITION DEPLETION
	#S = A[i] - A[i]*gamma*rho[i]/(gamma*rho[i]+1.) - beta*np.sqrt((i-x_nest)**2)				#ERIC-IFD DEPLETION
	#S = A[i] - A[i]*gamma*rho[i]/(A[i]+1.) - beta*np.sqrt((i-x_nest)**2)						#SATURATING-REBECCA DEPLETION
	return S

num_patches = 100 
#N = 20
#gamma = 1.
beta = 0.
x_nest = int(num_patches/2)

t_total = 10000							#t_total should make the S at the final stage pretty uniform! Otherwise, increase t_total!
dt = float(N)/float(t_total)
x = np.linspace(0,num_patches,num=num_patches)

A,rho,rho2 = [],[],[]
for i in range(num_patches):
	A.append(float(i)/50.)
	rho.append(0.)
	rho2.append(0.)
	
print "Food",np.sum(A),"Consumption",N*gamma
if np.sum(A)<N*gamma:
	print np.sum(A),N*gamma
	print "ERROR!!!! too much consumption"
###################################### ALGORITHM #########################################	

for t in range(t_total):
	S = []
	for i in range(num_patches):
		S.append(S_function(rho,i,gamma))
	best_x = np.argmax(S)
	rho[best_x] += dt
	
for t in range(t_total):
	S2 = []
	for i in range(num_patches):
		S2.append(S_function(rho2,i,2*gamma))
	best_x = np.argmax(S2)
	rho2[best_x] += dt

unif=[]
for i in range(len(x)):
	unif.append(1/float(num_patches))

##########################################################################################
######################################## PLOT ############################################
##########################################################################################

################### OPTIMAL

S = np.multiply(S, 1/np.sum(S))
rho = np.multiply(rho, 1/np.sum(rho))
A = np.multiply(A, 1/np.sum(A))
S2 = np.multiply(S2, 1/np.sum(S2))
rho2 = np.multiply(rho2, 1/np.sum(rho2))

#plt.plot(x,A,label="A")
#plt.plot(x,S,label="S")
plt.plot(x,rho,label="Max (Min($\epsilon _i$))",linewidth=2.,c="b",linestyle="--")
plt.plot(x,rho2,label="Max $\epsilon$",linewidth=2.,c="g",linestyle="--")
plt.plot(x,unif,label="Random",linewidth=2.,c="r",linestyle="--")
#plt.scatter(x,rho)
plt.legend(loc=2)
#plt.show()

#print "OPTIMAL ENERGY:", sum(rho), sum(rho*(F0-gamma*N*rho)), "OPTIMAL ENTROPY:",-sum(rho*np.log(rho+0.0000000000000001))
print "OPTIMAL ENERGY:", sum(rho), sum(rho*(F0-F0*gamma*N*rho/(F0+gamma*N*rho+eps))), "OPTIMAL ENTROPY:",-sum(rho*np.log(rho+0.0000000000000001))
rho_rand = np.full(num_patches,1/float(num_patches))
print "RANDOM ENERGY:", sum(rho_rand*(F0-gamma*N*rho_rand)), "OPTIMAL ENTROPY:",-sum(rho_rand*np.log(rho_rand))
#print "OPTIMAL ENERGY:", sum(rho), sum((F0-gamma*N*rho))
#print "OPTIMAL ENERGY:", sum(rho), sum((F0))
#print sum(rho), "OPTIMAL ENERGY:", sum(rho*(F0-F0*gamma*rho/(F0+gamma*rho+0.00001))),"OPTIMAL ENTROPY 1:",-sum((rho/float(N))*np.log((rho+0.000000000000001)/float(N)))

# print x.tolist()
#print A.tolist()
# print S.tolist()
#print rho.tolist()
	
################### Boltzmann

plt.plot(x,u_1,c="b",linewidth=2.,alpha=0.2,label=r"$\beta = 0.2$")
plt.plot(x,u_2,c="b",linewidth=2.,alpha=0.4,label=r"$\beta = 0.5$")
plt.plot(x,u_3,c="b",linewidth=2.,alpha=0.6,label=r"$\beta = 2.0$")
plt.plot(x,u_4,c="b",linewidth=2.,alpha=0.8,label=r"$\beta = 5.0$")
plt.plot(x,u_5,c="b",linewidth=2.,alpha=1.,label=r"$\beta = 10.$")
#plt.plot(x,u_5,c="b",linewidth=2.,alpha=1.,label="Boltzmann")
#plt.plot(x,IFD,label="IFD",linewidth=5.,c="black",linestyle="--")
plt.legend(loc=2)
plt.xlabel("Patch quality")
plt.ylabel("Density foragers")
plt.savefig("Boltzmann_energies.png",dpi=300)
plt.show()

print "BOLTZMANN ENERGY 1:", sum(u_1), sum(u_1*(F0-gamma*N*u_1))
print "BOLTZMANN ENERGY 2:", sum(u_2), sum(u_2*(F0-gamma*N*u_2))
print "BOLTZMANN ENERGY 3:", sum(u_3), sum(u_3*(F0-gamma*N*u_3))
# print "BOLTZMANN ENERGY 1:", sum(u_1), sum((F0-gamma*N*u_1))
# print "BOLTZMANN ENERGY 2:", sum(u_2), sum((F0-gamma*N*u_2))
# print "BOLTZMANN ENERGY 3:", sum(u_3), sum((F0-gamma*N*u_3))
#print "BOLTZMANN ENERGY 1:", sum(u_1), sum(u_1*(F0-F0*gamma*N*u_1/(F0+gamma*N*u_1+0.00001))),"BOLTZMANN ENTROPY 1:",-sum(u_1*np.log(u_1))
#print "BOLTZMANN ENERGY 2:", sum(u_2), sum(u_2*(F0-F0*gamma*N*u_2/(F0+gamma*N*u_2+0.00001))),"BOLTZMANN ENTROPY 2:",-sum(u_2*np.log(u_2))
#print "BOLTZMANN ENERGY 3:", sum(u_3), sum(u_3*(F0-F0*gamma*N*u_3/(F0+gamma*N*u_3+0.00001))),"BOLTZMANN ENTROPY 3:",-sum(u_3*np.log(u_3))

################### Modified Boltzmann

plt.plot(x,rho,label="Optimal",linewidth=2.,c="b",linestyle="--")
plt.plot(x,rho2,label="Optimal 2",linewidth=2.,c="g",linestyle="--")
plt.plot(x,unif,label="Random",linewidth=2.,c="r",linestyle="--")

plt.plot(x,u_1_m,c="g",linewidth=2.,alpha=0.2,label=r"$\beta = 0.2$")
plt.plot(x,u_2_m,c="g",linewidth=2.,alpha=0.4,label=r"$\beta = 0.5$")
plt.plot(x,u_3_m,c="g",linewidth=2.,alpha=0.6,label=r"$\beta = 2.0$")
plt.plot(x,u_4_m,c="g",linewidth=2.,alpha=0.8,label=r"$\beta = 5.0$")
plt.plot(x,u_5_m,c="g",linewidth=2.,alpha=1.,label=r"$\beta = 10.$")
#plt.plot(x,u_5_m,c="g",linewidth=2.,alpha=1.,label="Depleted Boltzmann")

#plt.plot(x,IFD,label="IFD",linewidth=5.,c="black",linestyle="--")
plt.legend(loc=2)
plt.xlabel("Patch quality")
plt.ylabel("Density foragers")
plt.savefig("Depleted_Boltzmann_energies.png",dpi=300)

print "MODIFIED BOLTZMANN ENERGY 1:", sum(u_1_m), sum(u_1_m*(F0-gamma*N*u_1_m))
print "MODIFIED BOLTZMANN ENERGY 2:", sum(u_2_m), sum(u_2_m*(F0-gamma*N*u_2_m))
print "MODIFIED BOLTZMANN ENERGY 3:", sum(u_3_m), sum(u_3_m*(F0-gamma*N*u_3_m))
# print "BOLTZMANN ENERGY 1:", sum(u_1), sum((F0-gamma*N*u_1))
# print "BOLTZMANN ENERGY 2:", sum(u_2), sum((F0-gamma*N*u_2))
# print "BOLTZMANN ENERGY 3:", sum(u_3), sum((F0-gamma*N*u_3))
#print "MODIFIED BOLTZMANN ENERGY 1:", sum(u_1_m), sum(u_1_m*(F0-F0*gamma*N*u_1_m/(F0+gamma*N*u_1_m+0.00001))),"MODIFIED BOLTZMANN ENTROPY 1:",-sum(u_1_m*np.log(u_1_m))
#print "MODIFIED BOLTZMANN ENERGY 2:", sum(u_2_m), sum(u_2_m*(F0-F0*gamma*N*u_2_m/(F0+gamma*N*u_2_m+0.00001))),"MODIFIED BOLTZMANN ENTROPY 2:",-sum(u_2_m*np.log(u_2_m))
#print "MODIFIED BOLTZMANN ENERGY 3:", sum(u_3_m), sum(u_3_m*(F0-F0*gamma*N*u_3_m/(F0+gamma*N*u_3_m+0.00001))),"MODIFIED BOLTZMANN ENTROPY 3:",-sum(u_3_m*np.log(u_3_m))

plt.show()

#print "Solution:"
#print u.tolist()
#print IFD
#print F0.tolist()
#print "Error",equations((x, y))

