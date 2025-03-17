from exponential import *

def delta_r12(Amplitudes,params):
    ais, ajs=np.meshgrid(params[0],params[0])
    bis, bjs=np.meshgrid(params[1],params[1])
    cis, cjs=np.meshgrid(params[2],params[2])
    Ampis,Ampjs=np.meshgrid(Amplitudes,Amplitudes)
    #permute the a,b, cs to get the correct formula 
    terms=term(cis,bis,ais,cjs,bjs,ajs)
    sum=np.sum(np.conjugate(Ampis)*Ampjs*terms)
    return sum/(4*np.pi)

def delta_r13(Amplitudes,params):
    ais, ajs=np.meshgrid(params[0],params[0])
    bis, bjs=np.meshgrid(params[1],params[1])
    cis, cjs=np.meshgrid(params[2],params[2])
    Ampis,Ampjs=np.meshgrid(Amplitudes,Amplitudes)
    terms=term(ais,cis,bis,ajs,cjs,bjs)
    sum=np.sum(np.conjugate(Ampis)*Ampjs*terms)
    return sum/(4*np.pi)

Rrms=mp.mpf(cnst.physical_constants["alpha particle rms charge radius"][0])/cnst.physical_constants["atomic unit of length"][0]

def E_ns(ai,bi,ci,aj,bj,cj):
    return (mp.mpf(32)/mp.mpf(3))*mp.pi**2*Rrms**2*(term(ci,bi,ai,cj,bj,aj)+term(ci,bi,ai,cj,bj,aj))