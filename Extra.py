from header import *

crd=3 #cube root dim
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

best=binary_search(lambda val:find_E(val,jkms)[0,0])

print((best[1]+79.005168)/79.005168)