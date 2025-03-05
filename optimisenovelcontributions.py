from exponential import *
from scipy.optimize import minimize
from time import time

params=np.load('data/altparams.npy',allow_pickle=True)

def f(x,oldN_eigenvalues,oldN_eigenvectors):
    if x[0]<=0 or x[1]<=0 or x[2]<=0:
        return 10000

    newparams=np.append(params,np.expand_dims(np.array(x),1)*mp.mpf(1),axis=1)
    subspace=Subspace(newparams.shape[1])

    subspace.set_N_func(N_func)
    subspace.set_H_func(H_func)

    subspace.set_params(newparams)

    subspace.make_N_mat()

    subspace.make_H_mat()
    subspace.find_N_eigens()
    subspace.make_Y_mat()
    subspace.make_invs_sqrt_beta_mats()
    subspace.make_P_mats()
    subspace.find_P_eigens()
    subspace.find_energy_eigenstates()

    novel_contribution=np.sum(
        subspace.energy_eigenstates[0]**2
        *(
            np.array(subspace.N_mat.tolist())[:,-1]
            -np.sum(
                np.array((oldN_eigenvectors.T*subspace.N_mat[:subspace.dim-1,:subspace.dim]).tolist()).transpose()*(
                np.array((oldN_eigenvectors.T*subspace.N_mat[:subspace.dim-1,subspace.dim-1]).tolist()).flatten()/
                np.array(oldN_eigenvalues.tolist()).flatten())
            ,axis=1)
        )**2
    )
    return -novel_contribution


for i in range(200):
    starttime=time()
    print()
    print(params.shape[1],"+1")
    x0=np.abs(np.array([400,1,0.02])*np.exp(np.random.normal(0,3,3)))
    print(x0)

    #first calculate N matrix without the new basis to get an orthogonal basis
    subspace=Subspace(params.shape[1])
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
    subspace.find_energy_eigenstates()
    subspace.find_energy_levels()
    print("Ground State energy level:",subspace.energy_levels[0])
    #calculate the exectation of the delta 
    expdelta=delta(subspace.energy_eigenstates[0],params)
    print("the expectation of the delta:",expdelta)

    res=minimize(f,x0=x0,method='Nelder-Mead',args=(subspace.N_eigenvalues,subspace.N_eigenvectors),options={'xatol':0.01})
    params=np.append(params,np.expand_dims(np.array(res.x),1)*mp.mpf(1),axis=1)
    print(res)
    print("time",time()-starttime)
    np.save("data/altparams",params,allow_pickle=True)
