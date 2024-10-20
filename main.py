import numpy as np
import scipy.linalg as la
import plotter


def integrate(function, N, width):
	random_inputs=(np.random.random((N,6))-0.5)*2*width

	return (64*width**6)*np.sum(function(random_inputs))/N

def gaussian_integrate(function, N, stddev):

	random_inputs=(np.random.normal(0,stddev,(N,6)))

	return np.sum(function(random_inputs)*((stddev**6)*8*np.pi**3)*np.exp(((np.linalg.norm(random_inputs,axis=-1)/stddev)**2)/2))/N

def acurate_integrate(function, N_integrals, N_points,width):
	#repeatedly do the integral and return the mean and standard error

	values=np.empty((N_integrals))
	for i in range(N_integrals):
		values[i]=integrate(function, N_points, width)
		print(i,values[i],np.mean(values[0:i+1]),np.std(values[0:i+1])/np.sqrt(i))
	return np.mean(values),np.std(values)/np.sqrt(N_integrals-1)

def acurate_gaussian_integrate(function, N_integrals, N_points, stddev):

	values=np.empty((N_integrals))
	for i in range(N_integrals):
		values[i]=gaussian_integrate(function, N_points, stddev)
		#print(i,values[i],np.mean(values[0:i+1]),np.std(values[0:i+1])/np.sqrt(i))
	return np.mean(values),np.std(values)/np.sqrt(N_integrals-1)

#for an array of values of kappa
k=1

#a) calculate matrix N

#define wavefunctions phi
Z=2
r0=0.52917#anstroms
ijms=np.array([[0,0,0]])#,[0,0,1],[0,1,0],[1,0,0]])
dim=np.shape(ijms)[0]
def phi(ijm,kappa,R):
	r1_norm=np.sqrt(R[:,0]**2+R[:,1]**2+R[:,2]**2)
	r2_norm=np.sqrt(R[:,3]**2+R[:,4]**2+R[:,5]**2)
	r12_norm=np.sqrt((R[:,0]-R[:,3])**2+(R[:,1]-R[:,4])**2+(R[:,2]-R[:,5])**2)
	return (r1_norm+r2_norm)**ijm[0]*(r1_norm-r2_norm)**ijm[1]*r12_norm**ijm[2]*np.exp(-Z*(r1_norm+r2_norm)/(kappa*r0))

#calculate the matrix N, where N_nn'=<phi_n|phi_n'>

#start by generating every combination of ns and n_primes in the matrix
ns=np.arange(0,dim)
n_primes=ns

#calculate each element of the array by doing the integral associated with each 
#using a loop here to avoid excessive memory usage since each element requires an integral with 1000000s of paramiters
print("Calculating elements of N:")
N=np.empty((dim,dim,2))

for n in ns:
	for n_prime in n_primes:
		N[n,n_prime]=acurate_gaussian_integrate(
			lambda R: phi(ijms[n],k,R)*phi(ijms[n_prime],k,R),
			N_integrals=10,
			N_points=100000, #10million is a sensibe max
			stddev=r0
		)
		print("\tN%s%s = <Φ%s|Φ%s> = %s ± %s" %(n,n_prime,n,n_prime,N[n,n_prime,0],N[n,n_prime,1]))


#b) calculate matices beta sqrt_beta inverse_sqrt_beta Y Y_T

#find eigenvalues and eigen vectors of N
print("Finding the eigen values and eigen vectors of N")
N_eigenvalues,N_eigenvectors = np.linalg.eig(N[:,:,0])

print("Constructing beta and Y")
#construct diagonal matrix +-sqrt beta from eigen values of N
invs_sqrt_beta=np.zeros((dim,dim))
beta=np.zeros((dim,dim))
for i in range(dim):
	invs_sqrt_beta[i,i]=1/(N_eigenvalues[i]**(0.5))
	beta[i,i]=N_eigenvalues[i]

#construct Y from eigen vectors of N
Y=N_eigenvectors

#c) calculate matrix H_tilde and then P

def laplaccian(function, R, h):
	return (1/h**2)*(function(R+2*h)-2*function(R+h)+function(R))

def H_opperator(function, R):
	r1_norm=np.sqrt(R[:,0]**2+R[:,1]**2+R[:,2]**2)
	r2_norm=np.sqrt(R[:,3]**2+R[:,4]**2+R[:,5]**2)
	r12_norm=np.sqrt((R[:,0]-R[:,3])**2+(R[:,1]-R[:,4])**2+(R[:,2]-R[:,5])**2)

	return laplaccian(function,R,0.00001)-Z/r1_norm-Z/r2_norm+1/r12_norm

#calculate elements of H_tilde
print("Calculating elements of H tilde:")
H=np.zeros((dim,dim,2))

for n in ns:
	for n_prime in n_primes:
		H[n,n_prime]=acurate_gaussian_integrate(
			lambda R: phi(ijms[n],k,R)*H_opperator(lambda R_prime:phi(ijms[n_prime],k,R_prime),R),
			N_integrals=10,
			N_points=100000, #10million is a sensibe max
			stddev=r0
		)
		print("\tH%s%s = <Φ%s|H|Φ%s> = %s ± %s" %(n,n_prime,n,n_prime,H[n,n_prime,0],H[n,n_prime,1]))


#calculate P using definition
print("Calculating P")
P=invs_sqrt_beta@np.transpose(Y)@H[:,:,0]@Y@invs_sqrt_beta

#c) Find the eigenvalues of P
energy_eigenvalues,_=np.linalg.eig(P)
print(energy_eigenvalues)

#for the ground state plot the energy against k