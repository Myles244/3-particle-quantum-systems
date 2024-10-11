import numpy as np
import scipy.linalg as la
import plotter
#for an array of values of kappa

#a) calculate matrix N

#define wavefunctions phi
Z=2
r0=0.52917#anstroms
def phi(i,j,m,kappa,r1,r2):
	r1_norm=la.norm(r1)
	r2_norm=la.norm(r2)
	return (r1_norm+r2_norm)**i*(r1_norm-r2_norm)**j*la.norm(r1-r2)**m*np.exp(-Z*(r1_norm+r2_norm)/(kappa*r0))

#test phi by plotting it


x,y=np.meshgrid(np.linspace(-1,1,1000),np.linspace(-1,1,1000))
r=np.sqrt(x**2+y**2)
#Phi=r*np.exp(-r/(2*r0))*np.exp(1j*np.angle(np.stack(x+y*1j,axis=-1)))
Phi=x+1j*y
plotter.plot(Phi)

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