from exponential import *
from scipy.optimize import minimize

params=np.load("data/bestparams.npy",allow_pickle=True)
bestparams=params

x0=params.flatten()
bestE=0

def f(x):
    global bestE
    global bestparams
    theseparams=x.reshape(3,-1)*mp.mpf(1)

    for i in range(theseparams.shape[1]):
        a,b,c=theseparams[:,i].flatten()
        if a<0 or b<0 or c<0 or (a==0 and b==0) or (a==0 and c==0) or (b==0 and c==0):
            return 1000
        
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
    
    if subspace.energy_levels[0]<bestE:
        print("Success, new lowest energy level found, saving paramiters")
        print(np.float64(subspace.energy_levels[0]))
        bestE=subspace.energy_levels[0]
        bestparams=theseparams
        #np.save("data/bestparams.npy",bestparams,allow_pickle=True)
    return np.float64((subspace.energy_levels[0]+402.637302)*1*10**6)

print(params.shape[1])

x0=params.flatten()

f(x0)

res=minimize(f,x0=x0)
print(res)
