from header import *
import scipy.constants as cnst

# --- constants --- 

#masses of nuclear particles:
        
#mass of the alpha particle
m1=(cnst.physical_constants['alpha particle mass'][0])*(10**(-12))**2/cnst.eV

#mass of the muon
m2=(cnst.physical_constants['muon mass'][0])*(10**(-12))**2/cnst.eV

#mass of the electron
m3=cnst.m_e*(10**(-12))**2/cnst.eV

#hbar
hbar=cnst.hbar/cnst.eV

#elementry charge
e=cnst.e

#vacume permiativity
epsilon0=(cnst.physical_constants['vacuum electric permittivity'][0])*(10**(-12))*cnst.eV


subspace=Subspace(1000)

def N_func(i,j,alphas,betas,gammas):
        
    return 512*(np.abs(alphas[i])*np.abs(betas[i])*np.abs(gammas[i])*np.abs(alphas[j])*np.abs(betas[j])*np.abs(gammas[j]))**1.5/((np.conj(alphas[i])+alphas[j])**3*(np.conj(betas[i])+betas[j])**3*(np.conj(gammas[i])+gammas[j])**3)

subspace.set_N_func(N_func)

def twocostheta(A,B,C):
	
	return 3*B/(2*A)+3*A/(2*B)-3*A*B/(C**2)
		

def H_func(i,j,alphas,betas,gammas):

	T=-0.5*hbar**2*(
            (alphas[j]**2+betas[j]**2)/m1  +  (alphas[j]**2+gammas[j]**2)/m2  +  (betas[j]**2+gammas[j]**2)/m3
            
            -2*alphas[j]*(1/m1+1/m2)*(np.conj(alphas[i])+alphas[j])/2
			-2*betas[j]*(1/m1+1/m3)*(np.conj(betas[i])+betas[j])/2	  
			-2*gammas[j]*(1/m2+1/m3)*(np.conj(gammas[i])+gammas[j])/2	  
                  
			+alphas[j]*betas[j]/m1*twocostheta(np.conj(alphas[i])+alphas[j],np.conj(betas[i])+betas[j],np.conj(gammas[i])+gammas[j])
			+alphas[j]*gammas[j]/m2*twocostheta(np.conj(betas[i])+betas[j],np.conj(gammas[i])+gammas[j],np.conj(alphas[i])+alphas[j])
			+betas[j]*gammas[j]/m3*twocostheta(np.conj(alphas[i])+alphas[j],np.conj(gammas[i])+gammas[j],np.conj(betas[i])+betas[j])

			)
      
	
	V=(e**2/(4*np.pi*epsilon0))*(-(np.conj(alphas[i])+alphas[j])-(np.conj(betas[i])+betas[j])+0.5*(np.conj(gammas[i])+gammas[j]))

	return (T+V)*N_func(i,j,alphas,betas,gammas)

subspace.set_H_func(H_func)

param_range=np.geomspace(0.01,100,10)
params=np.expand_dims(np.array(np.meshgrid(param_range,param_range,param_range,indexing='ij')).reshape(3,-1),-1)

subspace.set_params(1,params)

subspace.make_N_mats_vectorized()
subspace.make_H_mats_vectorized()
subspace.find_N_eigens()
subspace.make_Y_mats()
subspace.make_invs_sqrt_beta_mats()
subspace.make_P_mats()
subspace.find_P_eigens()
subspace.find_energy_levels()

print()