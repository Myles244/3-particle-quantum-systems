import numpy as np
import scipy.linalg as la
import plotter


def integrate(function, N, width):
	random_inputs=(np.random.random((N,6))-0.5)*2*width

	return (64*width**6)*np.sum(function(random_inputs))/N

def acurate_integrate(function, N_integrals, N_points,width):
	#repeatedly do the integral and return the mean and standard error

	values=np.empty((N_integrals))
	for i in range(N_integrals):
		values[i]=integrate(function, N_points, width)
		print(i,values[i],np.mean(values[0:i+1]),np.std(values[0:i+1])/np.sqrt(i))
	return np.mean(values),np.std(values)/np.sqrt(N_integrals-1)

#class State:
#	def __init__(self, wavefunction):
#		self.wavefunction = wavefunction

	#def norm(self):
	#	return acurate_integrate(prod,100,100000, 8, 0, 1)

#for an array of values of kappa
k=1

#a) calculate matrix N

#define wavefunctions phi
Z=2
r0=0.52917#anstroms
ijms=np.array([[0,0,0]])#,[0,0,1],[0,1,0],[1,0,0]])
def phi(ijm,kappa,R):
	r1_norm=np.sqrt(R[:,0]**2+R[:,1]**2+R[:,2]**2)
	r2_norm=np.sqrt(R[:,3]**2+R[:,4]**2+R[:,5]**2)
	r12_norm=np.sqrt((R[:,0]-R[:,3])**2+(R[:,1]-R[:,4])**2+(R[:,2]-R[:,5])**2)
	return (r1_norm+r2_norm)**ijm[0]*(r1_norm-r2_norm)**ijm[1]*r12_norm**ijm[2]*np.exp(-Z*(r1_norm+r2_norm)/(kappa*r0))

#create states with the wavefunctions phi_wavefuctions
#phis=State(lambda R:phi_wavefunction(ijms,1,R))

#calculate the matrix N, where N_nn'=<phi_n|phi_n'>

#start by generating every combination of ns and n_primes in the matrix
ns=np.arange(0,np.shape(ijms)[0])
n_primes=ns

#calculate each element of the array by doing the integral associated with each 
#using a loop here to avoid excessive memory usage since each element requires an integral with 1000000s of paramiters
N=np.empty((np.shape(ijms)[0],np.shape(ijms)[0],2))

for n in ns:
	for n_prime in n_primes:
		N[n,n_prime]=acurate_integrate(
			lambda R: phi(ijms[n],k,R)*phi(ijms[n_prime],k,R),
			N_integrals=100,
			N_points=10000000, #10million is a sensibe max
			width=3
		)

print(N)

#N=acurate_integrat e(lambda R: phi_wavefunction())

#b) calculate matices beta sqrt_beta inverse_sqrt_beta Y Y_T

#find eigenvalues and eigen vectors of N

#construct diagonal matrix +-sqrt beta from eigen values of N

#construct Y from eigen vectors of N

#c) calculate matrix H_tilde and then P

#calculate elements of H_tilde

#calculate P using definition

#c) Find the eigenvalues of P

#for the ground state plot the energy against k