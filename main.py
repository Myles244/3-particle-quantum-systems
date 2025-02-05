from header import *
import scipy.constants as cnst
import matplotlib.pyplot as plt
import math

# --- constants --- 

#masses of nuclear particles:
        
#mass of the alpha particle
m1=cnst.physical_constants['alpha particle mass'][0]

#mass of the muon
m2=cnst.physical_constants['muon mass'][0]

#mass of the electron
m3=cnst.physical_constants['electron mass'][0]

#hbar
hbar=cnst.physical_constants['reduced Planck constant'][0]

#elementry charge
e=cnst.physical_constants['elementry charge'][0]

#vacume permiativity
epsilon0=cnst.physical_constants['vacuum electric permittivity'][0]


subspace=Subspace(10)

def N_func(i,j,alphas,betas,gammas):
        
    return 512*np.pi**3/((np.conj(alphas[i])+alphas[j])**3*(np.conj(betas[i])+betas[j])**3*(np.conj(gammas[i])+gammas[j])**3)
	
subspace.set_N_func(N_func)

def H_func(i,j,alphas,betas,gammas):

	T=-0.5*hbar**2*((alphas[i]**2+betas[i]**2)/m1+(alphas[i]**2+gammas[i]**2)/m1+(betas[i]**2+gammas[i]**2)/m1)
	V=(e**2/(4*np.pi*epsilon0))*(-(np.conj(alphas[i])+alphas[j])-(np.conj(betas[i])+betas[j])+0.5*(np.conj(gammas[i])+gammas[j]))
	
	return (T+V)*N_func(i,j,alphas,betas,gammas)
	
subspace.set_H_func(H_func)



subspace.set_params(1000,np.meshgrid(np.linspace(0.5,1.6,1000))

subspace.make_N_mats()
subspace.make_H_mats()
subspace.find_N_eigens()
subspace.make_Y_mats()
subspace.make_invs_sqrt_beta_mats()
subspace.make_P_mats()
subspace.find_P_eigens()
subspace.find_energy_levels()

plt.plot(subspace.params[0],subspace.energy_levels[:,0])
plt.show()

