import numpy as np
import math
import matplotlib.pyplot as plt

# --- constants ---
#elementry charge squared
e2=14.4#eV Angstroms
#bohr radius
r0=0.53917#Angstroms
#charge on the nuclius is emelentry charge units
Z=2

#plot the expected curve when single basis state is used
kappas=np.linspace(0.5,3,1000)
def test_E0(kappa):
	return (e2/r0)*(4/(kappa**2)-27/(4*kappa))

#the ijm value associated with each n value
jkms=np.array([[0,0,0]])

#the dimension of the vector space created by the phi_n vectors
dim=np.shape(jkms)[0]

def P(kappa):
	#the N matrix
	N=np.zeros((dim,dim))
	#calculating elements
	for n in range(dim):
		for n_prime in range(dim):
			j,k,m=jkms[n]
			j_prime,k_prime,m_prime=jkms[n_prime]
			J=j+j_prime
			K=k+k_prime
			M=m+m_prime
			if M%2==0:
				N[n,n_prime]= 2*np.pi**2 * math.factorial(J+K+M+5) * (kappa*r0/(2*Z))**(J+K+M+6) * (1/(M+2)) * (1/(K+1) - 1/(K+3) - 1/(K+M+3) + (1/(K+M+5)))

	print(N)

P(1)

plt.plot(kappas,test_E0(kappas))
plt.show()
