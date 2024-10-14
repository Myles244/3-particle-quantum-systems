import numpy as np
import scipy.linalg as la
import plotter

def integrate(function, N, width, min, max):
	#intergates a real valued function over 4d space using the monte carlo method with N random points
	#each input dimenstion extends from -width to width
	#the output dimesnion extends from -max to max

	#could using gaussian random values reduce uncertant, uthis would eeffect the formlas tho
	#their might be a way to estemate the error on this integral without doing it a bunch check wikipedia

	#genrate N random input vectors
	random_inputs=(np.random.random((N,4))-0.5)*2*width

	#generate N random outputs
	random_outputs=min+(np.random.random((N)))*(max-min)

	#calculate the function at the random inputs and count it if the random point lies inside the function
	N_inside=np.count_nonzero(function(random_inputs)>random_outputs)

	#return the volume
	#which is the probability that it was inside the function times the toatl area then adding the offset from the min 
	return ((N_inside/N)*(max-min)+min)*16*width**4

def acurate_integrate(function, N_integrals, N_points,width, min, max):
	#repeatedly do the integral and return the mean and standard error

	values=np.empty((N_integrals))
	for i in range(N_integrals):
		values[i]=integrate(function, N_points, width,min, max)
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
ijms=np.array([[0,0,0],[0,0,1],[0,1,0],[1,0,0]])
def phi(ijm,kappa,R):
	r1_norm=np.sqrt(R[...,0]**2+R[...,1]**2)
	r2_norm=np.sqrt(R[...,2]**2+R[...,3]**2)
	r12_norm=np.sqrt((R[...,0]-R[...,2])**2+(R[...,1]-R[...,3])**2)
	return (r1_norm+r2_norm)**ijm[0]*(r1_norm-r2_norm)**ijm[1]*r12_norm**ijm[2]*np.exp(-Z*(r1_norm+r2_norm)/(kappa*r0))

#create states with the wavefunctions phi_wavefuctions
#phis=State(lambda R:phi_wavefunction(ijms,1,R))

#calculate the matrix N, where N_nn'=<phi_n|phi_n'>

#start by generating every combination of aloowed ijm values and storing them in an array where the index is the n/nprime
ijm_ns,ijm_n_primes = np.meshgrid(np.arange(0,np.shape(ijms)[0]),np.arange(0,np.shape(ijms)[0]))

#calculate each element of the array by doing the integral associated with each 
#using a loop here to avoid excessive memory usage since each element requires an integral with 1000000s of paramiters
N=np.empty((np.shape(ijms)[0],np.shape(ijms)[0],2))
for n,ijm_n in enumerate(ijm_ns):
	for n_prime,ijm_n_prime in enumerate(ijm_n_primes):
		N[n,n_prime]=acurate_integrate(
			lambda R: phi(ijm_n,k,R)*phi(ijm_n_prime,k,R),
			N_integrals=100,
			N_points=1000000,
			width=4,
			min=-1,
			max=1,
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