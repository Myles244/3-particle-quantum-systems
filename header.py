import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
import scipy.constants as cnst

# --- constants ---
#in units of angstroms and electron volts, where the prefactor 1/(4 pi epsilon0) has been absorbed into the definition of electric charge

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

#generator function of elements of N for hylleraas' basis states given J K M and one over two lambda or 'ootl'
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

#generator function of elements of C for hylleraas' basis states given J K M and one over two lambda or 'ootl'
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

#generator function of elements of T for hylleraas' basis states given J K M and one over two lambda or 'ootl'
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

#finds energy egien values (and their corosponding eigen vectors in the projected hilber space if return_as=True) given a kapp and a let of the combitnations of jkm
def find_E(kappas,jkms,return_as=False):

	#the dimension of the sub space
	dim=np.shape(jkms)[0]

	#the N matrix defined as N_nn'=<psi_n|psi_n'>
	N=np.zeros((dim,dim,np.size(kappas)))

	#the H_tilde matrix defined as H_nn'=<>
	H_tilde=np.zeros((dim,dim,np.size(kappas)))

	#lambda
	l=Z/(kappas*r0)

	#usefull constand 1/2lambda ie 'one over two lambda'
	ootl=kappas*r0/(2*Z)

	
	#calculating elements of the N and H_tilde matrices
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

	#changing the N and H_tilde form matrices with lists as element to lists of matrices so that we can use broadcasting of matrix fucntions
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

	#construct the P matrix
	P=invs_sqrt_beta@Y_T@H_tilde@Y@invs_sqrt_beta

	#get energy eigen values from Ps
	energy_eigenvalues,zs=np.linalg.eig(P)	

	#sort the energy eigenvalues 
	order=np.argsort(energy_eigenvalues)
	sorted_energy_eigenvalues=np.zeros((np.size(kappas),dim))
	sorted_As=np.zeros(zs.shape)
	for i in range(np.size(kappas)):
		sorted_energy_eigenvalues[i]=energy_eigenvalues[i,order[i]]
	
	#return the aplitudes if return_as==true else just return the energy eigen values
	if return_as==False:
		return sorted_energy_eigenvalues
	
	else:
		As=Y@invs_sqrt_beta@zs
		for i in range(np.size(kappas)):
			sorted_As[i]=As[i,order[i]]
		return sorted_energy_eigenvalues,sorted_As

#use binary search like algorythem to find the optimal kappa
#we can do this because we know the funciton has a single global minimum and is smooth
def binary_search(f):

	#staring information
	#the minimum must be between furthest left and furthest right
	current_best=np.array([1,f(1)])
	furthest_left=np.array([0.8,f(0.1)])
	furthest_right=np.array([1.2,f(2)])

	#the uncertanty on kappa, at which the program will stop
	min_error=0.00001

	#whist we havent yet satisfied our minimum error
	while (np.abs(furthest_right[0]-furthest_left[0])>min_error):

		#find the value in between our furthest left and best guess
		#then by coparing it with the values we already know we can find a tighter constraint on kappa
		#ie a better furthest left/furthest right/best guess

		value_to_decide=(current_best+furthest_left)/2
		value_to_decide[1]=f(value_to_decide[0])

		#if it is a new minumum use it has the current best and the current best is the new furthest right
		if value_to_decide[1]<current_best[1]:
			furthest_right=current_best
			current_best=value_to_decide

		#if it isint a new minmum then the minmum cant be to the left of it so we can set it to the new furthest left
		elif value_to_decide[1]>=current_best[1]:
			furthest_left=value_to_decide


		# now we can do the same for the furthest right
		value_to_decide=(current_best+furthest_right)/2
		value_to_decide[1]=f(value_to_decide[0])

		if value_to_decide[1]<current_best[1]:
			furthest_left=current_best
			current_best=value_to_decide

		elif value_to_decide[1]>=current_best[1]:
			furthest_right=value_to_decide

		print(current_best)
	return current_best

#these are the basis funcitons but their not actually being used for anything atm
def phi(jkl,kappa,r1,r2):
	r1_norm=np.linalg.norm(r1,axis=-1)
	r2_norm=np.linalg.norm(r2,axis=-1)
	r12_norm=np.linalg.norm(r2-r1,axis=-1)
	return (
		(r1_norm+r2_norm)**jkl[0]
		*(r1_norm+r2_norm)**jkl[1]
		*(r12_norm)**jkl[2]
		*np.exp((-Z/(kappa*r0))*(r1_norm+r2_norm))
		)
