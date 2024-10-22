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
#reduced plancs constant
hbar=6.5821220*10**(-16)#eV
#mass of electorn
me=9.11*10**(-31)#kg

#plot the expected curve when single basis state is used
kappas=np.linspace(0.5,3,1000)
def test_E0(kappa):
	return (e2/r0)*(4/(kappa**2)-27/(4*kappa))

#the ijm value associated with each n value
jkms=np.array([[0,0,0]])

#the dimension of the vector space created by the phi_n vectors
dim=np.shape(jkms)[0]

#generator function of element of N given JKM
def fN(J,K,M,ootl):
	if M==-2:
		return 0
	elif K%2 ==0:
		return 2*np.pi**2 * math.factorial(J+K+M+5) * ootl**(J+K+M+6) * (1/(M+2)) * (1/(K+1) - 1/(K+3) - 1/(K+M+3) + (1/(K+M+5)))
	else:
		return 0

#same for C
def fC(J, K, M, ootl):
	if M==-2:
		return 0
	elif K%2==0:
		return -Z*e2 * 8*np.pi**2 * math.factorial(J+K+M+4) * ootl**(J+K+M+5) * (1/(M+2)) * (1/(1+K)-1/(K+M+3))
	else:
		return 0

#same for W
def fW(J,K,M, ootl):
	return e2*fN(J,K,M-1,ootl)

#same for T
def fT(j,j_prime,k,k_prime,m,m_prime,J,K,M,l,ootl):
	return (hbar**2)/(2*me) * ( 2 * ((l**2)*fN(J,K,M,ootl) - J*l*fN(J-1,K,M,ootl) + j*j_prime*fN(J-2,K,M,ootl) + k*k_prime*fN(J,K-2,M,ootl) + m* m_prime*fN(J,K,M-2,ootl)) +(1/2)*((-M*l*(fC(J,K,M,ootl)-fC(J,K+2,M-2,ootl))) + (m*j_prime + m_prime*j)*(fC(J-1,K,M,ootl)-fC(J-1,K+2,M-2,ootl)) + (m*k_prime + m_prime*k)*(fC(J+1,K,M-2,ootl)-fC(J-1,K,M,ootl))))

def P(kappa):

	#the N matrix
	N=np.zeros((dim,dim))

	#the H_tilde matrix
	H_tilde=np.zeros((dim,dim))

	#lambda
	l=Z/(kappa*r0)

	#usefull constand 1/2lambda ie 'one over two lambda'
	ootl=kappa*r0/(2*Z)

	
	#calculating elements
	for n in range(dim):
		for n_prime in range(dim):
			j,k,m=jkms[n]
			j_prime,k_prime,m_prime=jkms[n_prime]
			J=j+j_prime
			K=k+k_prime
			M=m+m_prime
			N[n,n_prime]=fN(J,K,M, ootl)
			H_tilde[n,n_prime]=fC(J,K,M, ootl)+fW(J,K,M, ootl)+fT(j,j_prime,k,k_prime,m,m_prime,J,K,M, l,ootl)
			#print(fT(j,j_prime,k,k_prime,m,m_prime,J,K,M,l,ootl))

	#finding the eigen values of N
	N_eigenvalues, N_eigenvectors=np.linalg.eig(N)

	#calculating invs sqrt beta
	invs_sqrt_beta=np.zeros((dim,dim))
	for n in range(dim):
		invs_sqrt_beta[n,n]=1/np.sqrt(N_eigenvalues[n])
	
	#calculating Y
	Y=N_eigenvectors
	print(N)
	print(H_tilde)
	return invs_sqrt_beta@np.transpose(Y)@H_tilde@Y@invs_sqrt_beta
myP=P(1)
print(myP)
E,_=np.linalg.eig(myP)

print(E)
plt.scatter([1],E/2)
plt.plot(kappas,test_E0(kappas))
plt.show()
