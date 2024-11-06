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
hbar=6.5821220*10**(-16)#eVs
#mass of electorn
me=hbar**2/(r0*e2)

#the ijm value associated with each n value
crd=2 #cube root dim
pjkmv=np.arange(0,crd) # posible jkm values
jkms=np.zeros((crd**3,3),dtype=int)
for j in range(crd):
	for k in range(crd):
		for m in range(crd):
			jkms[(crd**2)*j+crd*k+m]=[j,k,m]

jkms=[[0,0,0],[0,2,0],[0,0,1]]

#the dimension of the vector space created by the phi_n vectors
dim=np.shape(jkms)[0]

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

def P(kappa):

	#the N matrix
	N=np.zeros((dim,dim,np.size(kappa)))

	#the H_tilde matrix
	H_tilde=np.zeros((dim,dim,np.size(kappa)))

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
			#print((n,n_prime))
			#print("\tC:",fC(J,K,M,ootl))
			#print("\tW:",fW(J,K,M,ootl))
			#print("\tT:",fT(j,j_prime,k,k_prime,m,m_prime,J,K,M,l,ootl))
			H_tilde[n,n_prime]=-Z*e2*fC(J,K,M, ootl)+e2*fN(J,K,M-1, ootl)+fT(j,j_prime,k,k_prime,m,m_prime,J,K,M, l,ootl)
			#print(fT(j,j_prime,k,k_prime,m,m_prime,J,K,M,l,ootl))

	#changing the N and H_tilde form matrices with lists as element to lists of matrices
	N=np.transpose(N)
	H_tilde=np.transpose(H_tilde)

	print()

	#print(H_tilde)
	#print(N)

	#finding the eigen values of N
	N_eigenvalues, N_eigenvectors=np.linalg.eig(N)

	#calculating invs sqrt beta
	invs_sqrt_beta=np.zeros((dim,dim,np.size(kappa)))
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
	sorted_energy_eigenvalues=np.zeros((kappas.size,dim))
	sorted_zs=np.zeros(np.shape())
	for i in range(kappas.size):
		sorted_energy_eigenvalues[i]=energy_eigenvalues[i,order[i]]
		sorted_zs[i]=zs[i,order[i]]
	

#this is a change blah blah blah
#this is another one

kappas=np.linspace(0.8,1.5,10000)




plt.scatter(kappas,sorted_energy_eigenvalues[:,0],color='C'+str(i),marker='.')

#find optimal k
kappas_order=np.argsort(sorted_energy_eigenvalues[:,0])
optimal_kappa=kappas[kappas_order[0]]
optimal_energy_eigenvalues=sorted_energy_eigenvalues[kappas_order[0]]
#plt.plot([optimal_kappa,optimal_kappa],[-100,150],'--')

#plt.scatter(kappas,sorted_test_energy_eigenvalues[:,1])



#plot the expected curve when single basis state is used
def test_E0(kappa):
	return (e2/r0)*(4/(kappa**2)-27/(4*kappa))

#plt.plot(kappas,test_E0(kappas))

print("Optimal value for kappa:",optimal_kappa)
print("Coresponing energy eigenvalues:",optimal_energy_eigenvalues)

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

def phi_0(r):
	value=0
	for n in range(dim):
		value+=

plt.show()
