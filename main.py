from exponential import *

params=[]
N=[]
E=[]
A=[]

print()

dim=100

#generate random params
params=np.abs(np.diag([mp.mpf(400),mp.mpf(1),mp.mpf(0.02)])@np.random.normal(1,1,(3,dim)))

#calculate the energy levels
subspace=Subspace(dim,verbose=True)

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
def term(ai,bi,ci,aj,bj,cj):
    Ai=2*np.real(ai)+2*np.real(bi)
    Bi=2*np.real(ai)+2*np.real(ci)
    Ci=2*np.real(bi)+2*np.real(ci)
    Aj=2*np.real(aj)+2*np.real(bj)
    Bj=2*np.real(aj)+2*np.real(cj)
    Cj=2*np.real(bj)+2*np.real(cj)

    return (
        np.sqrt(Ai*Aj)/(np.conjugate(ai+bi)+aj+bj)
        *np.sqrt(Bi*Bj)/(np.conjugate(ai+bi)+aj+bj)
        *np.sqrt(Ci*Cj)/(np.conjugate(ai+bi)+aj+bj)
        *Ai*Bi*Ci/np.sqrt( Ai**2 * (Bi + Ci) + Bi**2 * (Ai + Ci) + Ci**2 * (Ai + Bi) + Ai*Bi*Ci)
        *Aj*Bj*Cj/np.sqrt( Aj**2 * (Bj + Cj) + Bj**2 * (Aj + Cj) + Cj**2 * (Aj + Bj) + Aj*Bj*Cj)
    )

def delta(Amplitudes,params):
    ais, ajs=np.meshgrid(params[0],params[0])
    bis, bjs=np.meshgrid(params[1],params[1])
    cis, cjs=np.meshgrid(params[2],params[2])
    Ampis,Ampjs=np.meshgrid(Amplitudes,Amplitudes)
    terms=term(ais,bis,cis,ajs,bjs,cjs)
    sum=np.sum(np.conjugate(Ampis)*Ampjs*terms)
    return sum/(4*np.pi)

expdelta=delta(subspace.energy_eigenstates[0],params)

print("\nthe expectation of the delta:",expdelta)

#calculate the hyperfine splitting
alpha=mp.mpf(cnst.alpha)

g2=mp.mpf(2.002331846)

g3=mp.mpf(2.002319304386)

MHz_conversion_factor=mp.mpf(6.57968392061)*(mp.mpf(10)**mp.mpf(9))

def HFS(expdelta):
    prefactor=(
        (2*np.pi/(3))
        *(four_pi_epsilon0*hbar**4/(e**2))
        *alpha**2
        *(g2*g3/(m2*m3))
        )*MHz_conversion_factor
    return prefactor*expdelta

print("\nThe hyperfine splitting:",HFS(expdelta),"\n")
