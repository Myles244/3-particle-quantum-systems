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

	print()

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



x1, x2, y1, y2 = 0.8, 1.5, -80, -70 
subkappas=np.linspace(0.8,1.5,10000)
jkms=[[0,0,0],[0,0,1],[0,2,0]]

subE=find_E(subkappas,jkms)

kappas=np.linspace(0.5,20,10000)

E=find_E(kappas,jkms)

#find optimal k
kappas_order=np.argsort(subE[:,0])
optimal_kappa=subkappas[kappas_order[0]]
optimal_energy_eigenvalues=subE[kappas_order[0]]
#plt.plot([optimal_kappa,optimal_kappa],[-100,150],'--')

#plt.scatter(kappas,sorted_test_energy_eigenvalues[:,1])
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





fig, ax = plt.subplots()

ax.set_yticks([0,-20,-40,-60,-80])
ax.set_xticks([0,5,10,15,20])

ax.plot(kappas,E[:,0])

# inset Axes....
 # subregion of the original image
axins = ax.inset_axes(
    [0.5, 0.1, 0.4, 0.4],
    xlim=(x1, x2), ylim=(y1, y2), xticks=[optimal_kappa], yticks=[optimal_energy_eigenvalues[0]])
axins.set_yticklabels(["$-79.0\,eV$"])
axins.set_xticklabels(["$1.10$"])
axins.plot(subkappas,subE[:,0])
axins.scatter([optimal_kappa],[optimal_energy_eigenvalues[0]],marker='o')
axins.plot([0,optimal_kappa],[optimal_energy_eigenvalues[0],optimal_energy_eigenvalues[0]],color='C0',linestyle='--')
axins.plot([optimal_kappa,optimal_kappa],[-80,optimal_energy_eigenvalues[0]],color='C0',linestyle='--')

ax.indicate_inset((x1,y1,x2-x1,y2-y1), edgecolor="black")
ax.add_patch(ConnectionPatch(xyA=(x2, y1), xyB=(0, 0), axesA=ax, axesB=axins, coordsA="data", coordsB="axes fraction", lw=0.7,ls="--"))
ax.add_patch(ConnectionPatch(xyA=(x2, y2), xyB=(0, 1), axesA=ax, axesB=axins, coordsA="data", coordsB="axes fraction", lw=0.7,ls="--"))


ax.set_xlabel("$\kappa$")
ax.set_ylabel("$E/eV$")

plt.savefig("variation.png",dpi=300)


#the ijm value associated with each n value
crd=3 #cube root dim
pjkmv=np.arange(0,crd) # posible jkm values
jkms=np.zeros((crd**3,3),dtype=int)
for j in range(crd):
	for k in range(crd):
		for m in range(crd):
			jkms[(crd**2)*j+crd*k+m]=[j,2*k,m]

#print(jkms)

tempf=np.linalg.norm(jkms*[1,0.5,1],axis=-1)
#print(tempf)
jkms_order=np.argsort(tempf)
#print(jkms_order)
jkms=jkms[jkms_order]
#print(jkms)

kappas=np.linspace(0.2,4,10000)
jkmss=[[[0,0,0],[0,0,1],[0,2,0]],jkms[:4],jkms[:9],jkms[:81]]

Es=np.zeros((4,kappas.size))
for r in range(4):
	Es[r]=find_E(kappas,jkmss[r])[:,0]

	
fig, ax = plt.subplots()

ax.set_ylim((-85,5))
ax.set_xlim((0,4))
for r in range(4):
	ax.plot(kappas,Es[r])

ax.set_yticks([0,-20,-40,-60,-80])
ax.set_xticks([0,1,2,3])


ax.set_xlabel("$\kappa$")
ax.set_ylabel("$E/eV$")

axins = ax.inset_axes(
    [0.5, 0.5, 0.4, 0.4],
    xlim=(x1, x2), ylim=(y1, y2), xticks=[optimal_kappa], yticks=[optimal_energy_eigenvalues[0]])
axins.set_yticklabels(["$-79.0\,eV$"])
axins.set_xticklabels(["$1.10$"])

subEs=np.zeros((4,subkappas.size))
for r in range(4):
	subEs[r]=find_E(subkappas,jkmss[r])[:,0]

for r in range(4):
	axins.plot(subkappas,subEs[r])
#axins.scatter([optimal_kappa],[optimal_energy_eigenvalues[0]],marker='o')
#axins.plot([0,optimal_kappa],[optimal_energy_eigenvalues[0],optimal_energy_eigenvalues[0]],color='C0',linestyle='--')
#axins.plot([optimal_kappa,optimal_kappa],[-80,optimal_energy_eigenvalues[0]],color='C0',linestyle='--')



kappas_order=np.argsort(subEs[-1])
optimal_kappa=subkappas[kappas_order[0]]
optimal_energy_eigenvalues=subEs[-1,kappas_order[0]]

def binary_search(f):
	current_best=np.array([1,f(1)])
	furthest_left=np.array([0.1,f(0.1)])
	furthest_right=np.array([2,f(2)])
	min_error=0.00001
	while (np.abs(furthest_right[0]-furthest_left[0])>min_error):
		value_to_decide=(current_best+furthest_left)/2
		value_to_decide[1]=f(value_to_decide[0])

		if value_to_decide[1]<current_best[1]:
			furthest_right=current_best
			current_best=value_to_decide

		elif value_to_decide[1]>=current_best[1]:
			furthest_left=value_to_decide

		value_to_decide=(current_best+furthest_right)/2
		value_to_decide[1]=f(value_to_decide[0])

		if value_to_decide[1]<current_best[1]:
			furthest_left=current_best
			current_best=value_to_decide

		elif value_to_decide[1]>=current_best[1]:
			furthest_right=value_to_decide

		print(current_best)
	return current_best

best_kappa=binary_search(lambda val:find_E(val,jkmss[-1])[0,0])

axins.scatter(best_kappa[0],best_kappa[1],marker='o',color='C3')
#axins.plot([0,optimal_kappa],[optimal_energy_eigenvalues,optimal_energy_eigenvalues],color='C0',linestyle='--')
#axins.plot([optimal_kappa,optimal_kappa],[-80,optimal_energy_eigenvalues],color='C0',linestyle='--')

print("Optimal value for kappa:",optimal_kappa)
print("Coresponing energy eigenvalues:",optimal_energy_eigenvalues)

print("best_kappa",best_kappa)

ax.indicate_inset((x1,y1,x2-x1,y2-y1), edgecolor="black")
ax.add_patch(ConnectionPatch(xyA=(x1, y2), xyB=(0, 1), axesA=ax, axesB=axins, coordsA="data", coordsB="axes fraction", lw=0.7,ls="--"))
ax.add_patch(ConnectionPatch(xyA=(x2, y2), xyB=(0, 0), axesA=ax, axesB=axins, coordsA="data", coordsB="axes fraction", lw=0.7,ls="--"))



fig.savefig("convergence.png",dpi=300)

print(binary_search(lambda val:find_E(val,jkmss[0])[0,0]))