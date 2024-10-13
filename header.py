import numpy as np
import matplotlib.pyplot as plt
def integrate(function, N, width, max):
	#intergates a real valued function over 4d space using the monte carlo method with N random points
	#each input dimenstion extends from -width to width
	#the output dimesnion extends from -max to max

	#could using gaussian random values reduce uncertant, uthis would eeffect the formlas tho
	#their might be a way to estemate the error on this integral without doing it a bunch check wikipedia

	#genrate N random input vectors
	random_inputs=(np.random.random((N,4))-0.5)*2*width

	#generate N random outputs
	random_outputs=(np.random.random((N))-0.5)*2*max

	#calculate the function at the random inputs and count it if the random point lies inside the function
	N_inside=np.count_nonzero(function(random_inputs)>random_outputs)

	#return the volume
	#which is the probability that it was inside the function take 0.5 (because anything below the origin is negative) times the toatl area
	return (N_inside/N-0.5)*32*max*width**4

def acurate_integrate(function, N_integrals, N_points,width, max):
	#repeatedly do the integral and return the mean and standard error

	values=np.empty((N_integrals))
	for i in range(N_integrals):
		values[i]=integrate(function, N_points, width, max)
	return np.mean(values),np.std(values)/np.sqrt(N_integrals-1)

class State:
	def __innit__(self, wavefunction):
		self.wavefunction = wavefunction

	def norm(self,R):
		def prod(R):
			return np.linalg.abs(self.wavefunction(R))
		return acurate_integrate(prod,100,1000000, 4, 1)


class PseudoBasis:
	def __innit__(self, basis_states):
		self.basis_states = basis_states


