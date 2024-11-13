import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
import scipy.constants as cnst
# --- constants ---
#elementry charge squared
e2=(cnst.e**2/(4*np.pi*cnst.epsilon_0))/(cnst.eV*cnst.angstrom)
#bohr radius
r0=(4*np.pi*cnst.epsilon_0*cnst.hbar**2)/(cnst.angstrom*cnst.m_e*cnst.e**2)
#charge on the nuclius is emelentry charge units
Z=2
#reduced plancs constant
hbar=cnst.hbar/cnst.eV
#mass of electorn
me=cnst.m_e*cnst.angstrom**2/cnst.eV

#the dimension of the vector space created by the phi_n vectors

#generator function of element of N given JKM
def fN(J,K,M,ootl):
	if M==-2:
		return 0
	elif K%2 ==0:
		return (
			2*np.pi**2 
			* math.factorial(J+K+M+5) 
			* ootl**(J+K+M+6) 
			* (1/(M+2)) 
			* (1/(K+1) - 1/(K+3) - 1/(K+M+3) + (1/(K+M+5)))
		)
	else:
		return 0

#same for C
def fC(J, K, M, ootl):
	if M==-2:
		return 0
	elif K%2==0:
		return (
			8*np.pi**2 
			* math.factorial(J+K+M+4) 
			* ootl**(J+K+M+5) 
			* (1/(M+2)) 
			* (1/(1+K)-1/(K+M+3))
		)
	else:
		return 0

#same for T
def fT(j,j_prime,k,k_prime,m,m_prime,J,K,M,l,ootl):
	return (
		(hbar**2)/(2*me) * ( 
			2 * (
				(l**2)*fN(J,K,M,ootl) 
				- J*l*fN(J-1,K,M,ootl) 
				+ j*j_prime*fN(J-2,K,M,ootl) 
				+ k*k_prime*fN(J,K-2,M,ootl) 
				+ m* m_prime*fN(J,K,M-2,ootl)
				) 
			+(1/2)*(
				-M*l*(fC(J,K,M,ootl)-fC(J,K+2,M-2,ootl))
				+ (m*j_prime + m_prime*j)*(fC(J-1,K,M,ootl)-fC(J-1,K+2,M-2,ootl)) 
				+ (m*k_prime + m_prime*k)*(fC(J+1,K,M-2,ootl)-fC(J-1,K,M,ootl))
			)
		)
	)

def find_E(kappas,jkms):

	dim=np.shape(jkms)[0]

	#the N matrix
	N=np.zeros((dim,dim,np.size(kappas)))

	#the H_tilde matrix
	H_tilde=np.zeros((dim,dim,np.size(kappas)))

	#lambda
	l=Z/(kappas*r0)

	#usefull constand 1/2lambda ie 'one over two lambda'
	ootl=kappas*r0/(2*Z)

	
	#calculating elements
	for n in range(dim):
		for n_prime in range(dim):
			j,k,m=jkms[n]
			j_prime,k_prime,m_prime=jkms[n_prime]
			J=j+j_prime
			K=k+k_prime
			M=m+m_prime
			N[n,n_prime]=fN(J,K,M, ootl)
			#print((n,n_prime))
			#print("\tC:",fC(J,K,M,ootl))
			#print("\tW:",fW(J,K,M,ootl))
			#print("\tT:",fT(j,j_prime,k,k_prime,m,m_prime,J,K,M,l,ootl))
			H_tilde[n,n_prime]=-Z*e2*fC(J,K,M, ootl)+e2*fN(J,K,M-1, ootl)+fT(j,j_prime,k,k_prime,m,m_prime,J,K,M, l,ootl)
			#print(fT(j,j_prime,k,k_prime,m,m_prime,J,K,M,l,ootl))

	#changing the N and H_tilde form matrices with lists as element to lists of matrices
	N=np.transpose(N)
	H_tilde=np.transpose(H_tilde)

	#print(H_tilde)
	#print(N)

	#finding the eigen values of N
	N_eigenvalues, N_eigenvectors=np.linalg.eig(N)

	#calculating invs sqrt beta
	invs_sqrt_beta=np.zeros((dim,dim,np.size(kappas)))
	for n in range(dim):
		invs_sqrt_beta[n,n]=1/np.sqrt(N_eigenvalues[:,n])
	invs_sqrt_beta=np.transpose(invs_sqrt_beta)

	#calculating Y
	Y=N_eigenvectors
	Y_T=np.transpose(Y,axes=[0,2,1])

	#P	
	P=invs_sqrt_beta@Y_T@H_tilde@Y@invs_sqrt_beta

	#get energy eigen values from Ps
	energy_eigenvalues,zs=np.linalg.eig(P)	

	#sort the energy eigenvalues 
	order=np.argsort(energy_eigenvalues)
	sorted_energy_eigenvalues=np.zeros((np.size(kappas),dim))
	sorted_zs=np.zeros(zs.shape)
	for i in range(np.size(kappas)):
		sorted_energy_eigenvalues[i]=energy_eigenvalues[i,order[i]]
		sorted_zs[i]=zs[i,order[i]]
	
	return sorted_energy_eigenvalues