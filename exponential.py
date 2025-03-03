from header import *
import scipy.constants as cnst


# --- constants --- 

#masses of nuclear particles:
        
#mass of the alpha particle
m1=mp.mpf(cnst.physical_constants['alpha particle mass'][0])/mp.mpf(cnst.m_e)

#mass of the muon
m2=mp.mpf(cnst.physical_constants['muon mass'][0])/mp.mpf(cnst.m_e)

#mass of the electron
m3=mp.mpf(1)

#hbar
hbar=mp.mpf(1)

#elementry charge
e=mp.mpf(1)

#vacume permiativity
four_pi_epsilon0=mp.mpf(1)



#prefactor
def pref(lambdas):
    return (
        lambdas[0]/np.sqrt(lambdas[3]*lambdas[6])
        *lambdas[1]/np.sqrt(lambdas[4]*lambdas[7])
        *lambdas[2]/np.sqrt(lambdas[5]*lambdas[8])
        )

#messy norm term
def mnt(x,y,z): 
    return (x**2*(y+z)+y**2*(x+z)+z**2*(x+y)+x*y*z)

def I110(lambdas):

    I=(
        0.5
        *pref(lambdas)
        *((2*lambdas[0]**2 + lambdas[0]*lambdas[1] + lambdas[0]*lambdas[2] + lambdas[1]*lambdas[2])
        /(
            np.sqrt(mnt(lambdas[3],lambdas[4],lambdas[5]))
            *np.sqrt(mnt(lambdas[6],lambdas[7],lambdas[8])) 
        ))
        )
    
    return I

def I111(lambdas):

    I=(
        pref(lambdas)
        *(mnt(lambdas[0],lambdas[1],lambdas[2])
        /(
            np.sqrt(mnt(lambdas[3],lambdas[4],lambdas[5]))
            *np.sqrt(mnt(lambdas[6],lambdas[7],lambdas[8])) 
        ))
        )
    return I

def I210(lambdas):
        
    I=(
        pref(lambdas)
        *((lambdas[0]**2*(3*lambdas[0]+2*lambdas[1]+lambdas[2])+lambdas[1]**2*(lambdas[0]+lambdas[2]) + lambdas[0]*lambdas[1]*lambdas[2])
        /(
            np.sqrt(mnt(lambdas[3],lambdas[4],lambdas[5]))
            *np.sqrt(mnt(lambdas[6],lambdas[7],lambdas[8])) 
        ))
        )
    
    return I

def I300(lambdas):

    I=(
        3
        *pref(lambdas)
        *(lambdas[0]**2+lambdas[1]**2)/np.sqrt(mnt(lambdas[3],lambdas[4],lambdas[5]))
        *(lambdas[0]+lambdas[1])/np.sqrt(mnt(lambdas[6],lambdas[7],lambdas[8]))
        )

    return I

def N_func(i,j,alphas,betas,gammas):
    lambdas=[
        1/(np.conjugate(alphas[i])+alphas[j]+np.conjugate(betas[i])+betas[j]),1/(np.conjugate(alphas[i])+alphas[j]+np.conjugate(gammas[i])+gammas[j]),  1/(np.conjugate(betas[i])+betas[j]+np.conjugate(gammas[i])+gammas[j]),  
        0.5/(np.real(alphas[i])+np.real(betas[i])), 0.5/(np.real(alphas[i])+np.real(gammas[i])), 0.5/(np.real(betas[i])+np.real(gammas[i])),
         0.5/(np.real(alphas[j])+np.real(betas[j])), 0.5/(np.real(alphas[j])+np.real(gammas[j])), 0.5/(np.real(betas[j])+np.real(gammas[j]))
        ]
    I=I111(lambdas)
    
    return I

def order(first,second,third,lambdas):
    return [lambdas[first+second-1],lambdas[first+third-1],lambdas[second+third-1],
            lambdas[3],lambdas[4],lambdas[5],
            lambdas[6],lambdas[7],lambdas[8]
            ]


def H_func(i,j,alphas,betas,gammas):

	lambdas=[
		1/(np.conjugate(alphas[i])+alphas[j]+np.conjugate(betas[i])+betas[j]),1/(np.conjugate(alphas[i])+alphas[j]+np.conjugate(gammas[i])+gammas[j]),  1/(np.conjugate(betas[i])+betas[j]+np.conjugate(gammas[i])+gammas[j]),  
		0.5/(np.real(alphas[i])+np.real(betas[i])), 0.5/(np.real(alphas[i])+np.real(gammas[i])), 0.5/(np.real(betas[i])+np.real(gammas[i])),
		0.5/(np.real(alphas[j])+np.real(betas[j])), 0.5/(np.real(alphas[j])+np.real(gammas[j])), 0.5/(np.real(betas[j])+np.real(gammas[j]))
		]

	T=-0.5*hbar**2*(
			((alphas[j]**2+betas[j]**2)/m1  +  (alphas[j]**2+gammas[j]**2)/m2  +  (betas[j]**2+gammas[j]**2)/m3)*I111(lambdas)
			
			-2*alphas[j]*(1/m1+1/m2)*I110(order(1,2,0,lambdas))
			-2*betas[j]*(1/m1+1/m3)*I110(order(0,2,1,lambdas))
			-2*gammas[j]*(1/m2+1/m3)*I110(lambdas)
				
			+alphas[j]*betas[j]/m1*(I210(order(0,2,1,lambdas))+I210(order(1,2,0,lambdas))-I300(order(2,0,1,lambdas)))
			+alphas[j]*gammas[j]/m2*(I210(lambdas)+I210(order(2,1,0,lambdas))-I300(order(1,0,2,lambdas)))
			+betas[j]*gammas[j]/m3*(I210(order(1,0,2,lambdas))+I210(order(2,0,1,lambdas))-I300(lambdas))
			)
	
	V=(e**2/(four_pi_epsilon0))*(-2*I110(order(1,2,0,lambdas))-2*I110(order(0,2,1,lambdas))+I110(lambdas))

	return T+V

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

alpha=mp.mpf(cnst.alpha)

g2=mp.mpf(cnst.physical_constants["muon g factor"][0])

g3=mp.mpf(cnst.physical_constants["electron g factor"][0])

MHz_conversion_factor=mp.mpf(cnst.physical_constants["atomic unit of energy"][0])/(mp.mpf(cnst.h)*mp.mpf(1000000))

def HFS(expdelta):
    prefactor=(
        (2*np.pi/(3))
        *(four_pi_epsilon0*hbar**4/(e**2))
        *alpha**2
        *(g2*g3/(m2*m3))
        )*MHz_conversion_factor
    return prefactor*expdelta