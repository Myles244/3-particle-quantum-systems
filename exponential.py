from header import *
import scipy.constants as cnst

#set the precision
mp.dps = 100

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
