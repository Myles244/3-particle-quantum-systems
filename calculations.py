from exponential import *
from uncertainties import ufloat

mp.dps=100

#load in good params
params=np.load("data/best100params.npy",allow_pickle=True)

#generate random params
#randomparams=np.abs(np.diag([400,1,0.02])@np.exp(np.random.normal(0,5,(3,50))))*mp.mpf(1)

#theseparams=np.append(params,randomparams,axis=1)

#calculate the energy levels
subspace=Subspace(params.shape[1], verbose=True)

subspace.set_N_func(N_func)
subspace.set_H_func(H_func)

subspace.set_params(params)

subspace.make_N_mat()
subspace.make_H_mat()
subspace.find_N_eigens()
subspace.make_Y_mat()
subspace.make_invs_sqrt_beta_mats()
subspace.make_P_mats()
subspace.find_P_eigens()
subspace.find_energy_levels()
subspace.find_energy_eigenstates()

expdelta=delta_r23(subspace.energy_eigenstates[0],params)

#define masses and gs with uncertainty
def ufloat_from_const(name):
    temp=cnst.physical_constants[name]
    return ufloat(temp[0],temp[2])

energyunit=ufloat_from_const("Hartree energy in eV")
MHz_conversion_factor=ufloat_from_const("hartree-hertz relationship")/1000000
alpha=ufloat_from_const("fine-structure constant")
g2=ufloat_from_const("muon g factor")
g3=ufloat_from_const("electron g factor")
m3=1
m2=ufloat_from_const("muon-electron mass ratio")

def HFS(expdelta):
    prefactor=(
        (2*np.pi/(3))
        *alpha**2
        *(g2*g3/(m2*m3))
        )*MHz_conversion_factor
    return prefactor*expdelta

def energy_Hartrees():
    return np.float64(subspace.energy_levels[0])

def energy_KeV():
    return np.float64(subspace.energy_levels[0])*energyunit/1000

def expectation_of_the_delta():
    return "$"+str(ufloat_from_sample(expdeltas))+"$"

def hyperfine_splitting():
    return "$"+str(HFS(ufloat_from_sample(expdeltas)))+"$"