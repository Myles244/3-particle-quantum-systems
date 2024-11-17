import numpy as np
from header import *


crd=4 #cube root dim
pjkmv=np.arange(0,crd) # posible jkm values
jkms=np.zeros((crd**3,3),dtype=int)
for j in range(crd):
	for k in range(crd):
		for m in range(crd):
			jkms[(crd**2)*j+crd*k+m]=[j,2*k,m]

tempf=np.linalg.norm(jkms*[1,0.5,1],axis=-1)
#print(tempf)
jkms_order=np.argsort(tempf)
#print(jkms_order)
jkms=jkms[jkms_order]
#print(jkms)

best=binary_search(lambda val:find_E(val,jkms[:3])[0,0])

print((best[1]+79.005168)/79.005168)

kappa=best[0]

Es,As=find_E(best[0],jkms[:3],return_as=True)
E=Es[0,0]
A=As[0,:,0]
print(E)
print(A)

A=A/np.linalg.norm(A)

def P000(r):
	return (Z**3/(np.pi*r0**3*kappa**3))*np.exp(-(2*Z*r)/(r0*kappa))

def P020(r):
	return (2*Z**3/(81*np.pi*r0**7*kappa**7))*(r**2*Z**2-6*r*r0*Z*kappa+12*r0**2*kappa**2)**2*np.exp(-(2*Z*r)/(r0*kappa))

def P001(r):
	return (8*Z**3/(np.pi*r0**3*kappa**3))*np.exp(-(4*Z*r)/(r0*kappa))

rs=np.linspace(0,1,1000)

Ps=A[0]**2*P000(rs)+A[1]**2*P001(rs)+A[2]**2*P020(rs)

print(A)

#plt.plot(rs,Ps)
plt.plot(rs,np.sqrt(Ps))
plt.plot(rs,np.exp(-rs/r0)/(np.sqrt(np.pi)*r0**(3/2)))
plt.plot(rs,2*np.sqrt(2)*np.exp(-2*rs/r0)/(np.sqrt(np.pi)*r0**(3/2)))
plt.xticks([0,0.2,0.4,0.6,0.8,1])
plt.ylabel("amplitude")
plt.xlabel("r/r_0")
plt.savefig("wavefunction.png",dpi=300)
