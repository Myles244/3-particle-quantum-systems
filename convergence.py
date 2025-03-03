from exponential import *
from multiprocessing import Pool

N=[]
E=[]
D=[]
params=np.append(np.load("altparams.npy",allow_pickle=True),np.load("altaltparams.npy",allow_pickle=True),axis=1)


def f(i):
    subspace=Subspace(i)

    subspace.set_N_func(N_func)
    subspace.set_H_func(H_func)

    subspace.set_params(params[:,:i])

    subspace.make_N_mat()
    subspace.make_H_mat()
    subspace.find_N_eigens()
    subspace.make_Y_mat()
    subspace.make_invs_sqrt_beta_mats()
    subspace.make_P_mats()
    subspace.find_P_eigens()
    subspace.find_energy_levels()
    subspace.find_energy_eigenstates()
    expdelta=delta(subspace.energy_eigenstates[0],params[:,:i])
    print(i)
    return([i,subspace.energy_levels[0],expdelta])

ilist=range(1,params.shape[1])

data=[]

if __name__ == '__main__':
    with Pool() as p:
        data=(p.map(f, ilist))
        p.close()
        p.join()

        data=np.array(data)

        np.save("convergence",data,allow_pickle=True)




        







