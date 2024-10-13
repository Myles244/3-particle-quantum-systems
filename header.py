import numpy as np

def intergate(function, N, range, max):
	#intergates a real valued function over 4d space using the monte carlo method with N random points
	

	#genrate N random input vectors
	random_inputs=(np.random.random((N,4))-0.5)*2*range

	#generate N random outputs
	random_outputs=(np.random.random((N))-0.5)*2*max

	#calculate the function at the random inputs and count it if the random point lies inside the function
	N_inside=np.count_nonzero(function(random_inputs)>random_outputs)
	
	#return the volume, which is the probability of the point being inside, times the total volume that the points could be in
	return (N/N_inside)*16*max*range**4


class State:
	def __innit__(self, wavefunction):
		self.wavefunction = wavefunction

	def norm(self,R):
		phi=self.wavefunction(R)


class PseudoBasis:
	def __innit__(self, basis_states):
		self.basis_states = basis_states