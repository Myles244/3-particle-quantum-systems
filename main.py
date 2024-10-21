import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt

#monte carlo integrator
def integrate(function, N, stddev):
	#using a normal distrobtution in stead of uniform one has many advantages.
	#It increases the acuracy as the function is being integrated over the domain
	#It increases the precision since the are more points towards the center where most of the variance of the wavefunction is
	random_inputs=(np.random.normal(0,stddev,(N,6)))
	
	#The actual value is calculated as Sum(f(X)/p(X))/N. the 1/p(x) removes the bias on the values of f(x) from the distrobution of p(x)
	return np.sum(function(random_inputs)*((stddev**6)*8*np.pi**3)*np.exp(((np.linalg.norm(random_inputs,axis=-1)/stddev)**2)/2))/N

#does multiple integrations and returns the mean and standard error
def acurate_integrate(function, N_integrals, N_points, stddev):
	values=np.empty((N_integrals))
	for i in range(N_integrals):
		values[i]=integrate(function, N_points, stddev)
		#print(i,values[i],np.mean(values[0:i+1]),np.std(values[0:i+1])/np.sqrt(i))
	return np.mean(values),np.std(values)/np.sqrt(N_integrals-1)

def P(functions,hamiltonian):

	dim=np.shape(functions)[0]

	N=np.empty((dim,dim))

	H=np.empty((dim,dim))

	N_points=100000
	np.random.seed(0)
	random_inputs=np.random.normal(0,r0,(N_points,6))

	for n in range(dim):
		for n_prime in range(dim):
			F1=lambda R:functions[n](R)*functions[n_prime](R)*np.exp(0.5*(np.linalg.norm(R,axis=-1)/r0)**2)
			N[n,n_prime]=(8*r0**6*np.pi**3)*np.sum(F1(random_inputs))/N_points
			F2=lambda R:functions[n](R)*hamiltonian(lambda S:functions[n_prime](S),R)*np.exp(0.5*(np.linalg.norm(R,axis=-1)/r0)**2)
			H[n,n_prime]=(8*r0**6*np.pi**3)*np.sum(F2(random_inputs))/N_points
			#print("\tN%s%s = <Φ%s|Φ%s> = %s" %(n,n_prime,n,n_prime,N[n,n_prime]))

	#find eigenvalues and eigen vectors of N
	#print("Finding the eigen values and eigen vectors of N")
	N_eigenvalues,N_eigenvectors = np.linalg.eig(N)

	#print("Constructing beta and Y")
	#construct diagonal matrix +-sqrt beta from eigen values of N
	invs_sqrt_beta=np.zeros((dim,dim))
	beta=np.zeros((dim,dim))
	for i in range(dim):
		invs_sqrt_beta[i,i]=1/(N_eigenvalues[i]**(0.5))
		beta[i,i]=N_eigenvalues[i]

	#construct Y from eigen vectors of N
	Y=N_eigenvectors

	#c) calculate matrix H_tilde and then P

	#calculate P using definition
	#print("Calculating P")
	return invs_sqrt_beta@np.transpose(Y)@H@Y@invs_sqrt_beta


#define wavefunctions phi
#R is the combination of both r1 and r2
def phi(ijm,kappa,R):
	r1_norm=np.sqrt(R[:,0]**2+R[:,1]**2+R[:,2]**2)
	r2_norm=np.sqrt(R[:,3]**2+R[:,4]**2+R[:,5]**2)
	r12_norm=np.sqrt((R[:,0]-R[:,3])**2+(R[:,1]-R[:,4])**2+(R[:,2]-R[:,5])**2)
	return (r1_norm+r2_norm)**ijm[0]*(r1_norm-r2_norm)**ijm[1]*r12_norm**ijm[2]*np.exp(-Z*(r1_norm+r2_norm)/(kappa*r0))


def laplaccian(function, R, h):
	return (1/h**2)*(function(R+2*h)-2*function(R+h)+function(R))

def H_opperator(function, R):
	r1_norm=np.sqrt(R[:,0]**2+R[:,1]**2+R[:,2]**2)
	r2_norm=np.sqrt(R[:,3]**2+R[:,4]**2+R[:,5]**2)
	r12_norm=np.sqrt((R[:,0]-R[:,3])**2+(R[:,1]-R[:,4])**2+(R[:,2]-R[:,5])**2)

	return laplaccian(function,R,0.00001)-Z/r1_norm-Z/r2_norm+1/r12_norm

#kappa is the maramiter to vary, in order to minimise the energy levels
kappas=np.linspace(0.5,2,10)
Z=2#e
r0=0.52917#anstroms
ijms=np.array([[0,0,0]])#,[0,0,1],[0,1,0],[1,0,0]]
dim=np.shape(ijms)[0]
phis=[]
E=np.zeros((10,1,1))
for i in range(100):
	for j in range(dim):
		phis.append(lambda R,ijm=ijms[j],kappa=kappas[i]:phi(ijm,kappa,R))
	print(phis)
	Eigvalue,_=np.linalg.eig(P(phis,H_opperator))
	print(Eigvalue)
E=E.flatten()
plt.plot(kappas,E)
plt.show()
#for the ground state plot the energy against k