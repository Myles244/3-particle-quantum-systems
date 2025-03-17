from exponential import *
from time import time
print()

mp.dps=200

#generate random params
params=np.load("data/bestparams copy.npy",allow_pickle=True)

#calculate the energy levels

starttime=time()
subspace=Subspace(params.shape[1],verbose=True)

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

#print the ground energy level
print("\nGround State energy level:",subspace.energy_levels[0])

#calculate the exectation of the delta 
expdelta=delta_r23(subspace.energy_eigenstates[0],params)

print("\nthe expectation of the delta23:",expdelta)
from nuclear_structure_corection import *
print("\nthe expectation of the delta13:",delta_r13(subspace.energy_eigenstates[0],params))
print("\nthe expectation of the delta12:",delta_r12(subspace.energy_eigenstates[0],params))

#calculate the hyperfine splitting
print("\nThe hyperfine splitting:",HFS(expdelta),"\n")

print("time",starttime-time())