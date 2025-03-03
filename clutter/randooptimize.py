from exponential import *

mp.dps=50

#generate random params
params=np.load("params.npy",allow_pickle=True)[:,:200]

champion=(params,mp.mpf("-402.63730185402387328557354382533412485104065903481"))

E=[]

def f(x):

    global champion

    randoms=np.abs(np.diag([400,1,0.02])@np.exp(np.random.normal(0,2,(3,10))))*mp.mpf(1)

    myparams=np.append(champion[0][:,:190],randoms,axis=1)

    #calculate the energy levels
    subspace=Subspace(myparams.shape[1],verbose=True)

    subspace.set_N_func(N_func)
    subspace.set_H_func(H_func)

    subspace.set_params(myparams)

    subspace.make_N_mat()
    subspace.make_H_mat()
    subspace.find_N_eigens()
    subspace.make_Y_mat()
    subspace.make_invs_sqrt_beta_mats()
    subspace.make_P_mats()
    subspace.find_P_eigens()
    subspace.find_energy_levels()
    subspace.find_energy_eigenstates()
    expdelta=delta(subspace.energy_eigenstates[0],myparams)

    if subspace.energy_levels[0]<champion[1]:
        print(x,"has been victorious,\n\t E:",subspace.energy_levels[0],"\n\t D:",expdelta)
        order=np.argsort(-subspace.energy_eigenstates[0]**2)
        newparams=myparams[:,order]
        champion=(newparams,subspace.energy_levels[0])
        np.save("bestfirst200params",champion[0],allow_pickle=True)
        E.append(champion[1])

    else:
        print(x,'has failed')


for i in range(1000):
    f(i+1)

import matplotlib as plt

plt.plot(E)

plt.show()