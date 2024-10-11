import numpy as np
import scipy.linalg as la
import plotter
#for an array of values of kappa

#a) calculate matrix N

#define wavefunctions phi
Z=2
r0=0.52917#anstroms
def phi(i,j,m,kappa,R):
	r1_norm=np.linalg.norm(r2,axis=-1)
	r2_norm=np.linalg.norm(r2,axis=-1)
	return (r1_norm+r2_norm)**i*(r1_norm-r2_norm)**j*np.linalg.norm(r1-r2,axis=-1)**m*np.exp(-Z*(r1_norm+r2_norm)/(kappa*r0))

#test phi by plotting it
x=np.linspace(-2,2,100)
R=np.transpose(np.meshgrid(x,x,x,x))
Phi=10*phi(0,0,0,1,r)
#print(np.shape(Phi))
#plotter.plot(Phi)

#plt.show()
#calculate norms and save in matrix N

#b) calculate matices beta sqrt_beta inverse_sqrt_beta Y Y_T

#find eigenvalues and eigen vectors of N

#construct diagonal matrix +-sqrt beta from eigen values of N

#construct Y from eigen vectors of N

#c) calculate matrix H_tilde and then P

#calculate elements of H_tilde

#calculate P using definition

#c) Find the eigenvalues of P

#for the ground state plot the energy against k