from exponential import *
from uncertainties import ufloat

mp.dps=200

#load in good params
params=np.load("data/bestparams.npy",allow_pickle=True)[:,:50]

ground_state_energies=[]
expdeltas=[]

for i in range(100):
    print(i)
    #generate random params
    randomparams=np.abs(np.diag([400,1,0.02])@np.exp(np.random.normal(0,5,(3,50))))*mp.mpf(1)

    theseparams=np.append(params,randomparams,axis=1)

    #calculate the energy levels
    subspace=Subspace(theseparams.shape[1])

    subspace.set_N_func(N_func)
    subspace.set_H_func(H_func)

    subspace.set_params(theseparams)

    subspace.make_N_mat()
    subspace.make_H_mat()
    subspace.find_N_eigens()
    subspace.make_Y_mat()
    subspace.make_invs_sqrt_beta_mats()
    subspace.make_P_mats()
    subspace.find_P_eigens()
    subspace.find_energy_levels()
    subspace.find_energy_eigenstates()

    expdelta=delta(subspace.energy_eigenstates[0],theseparams)

    ground_state_energies.append(np.float64(subspace.energy_levels[0]))
    expdeltas.append(np.float64(expdelta))

#define masses and gs with uncertainty

def ufloat_from_sample(sample):
    return ufloat(np.mean(sample),np.std(sample,ddof=1)/np.sqrt(len(sample)))

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

import matplotlib.pyplot as plt
plt.hist(np.float64(ground_state_energies))
plt.show()
plt.hist(np.float64(expdeltas))
plt.show()

def ground_state_energy():
    return "$"+str(ufloat_from_sample(ground_state_energies)*energyunit/1000)+"$ KeV"

def expectation_of_the_delta():
    return "$"+str(ufloat_from_sample(expdeltas))+"$"

def hyperfine_splitting():
    return "$"+str(HFS(ufloat_from_sample(expdeltas)))+"$"