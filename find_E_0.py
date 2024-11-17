from header import *

#finds optimal kappa and calculates the error from experimental the experimental value

#generates all the jkm values used for the basis vectors
crd=3 #cube root dim
pjkmv=np.arange(0,crd) # posible jkm values
jkms=np.zeros((crd**3,3),dtype=int)
for j in range(crd):
	for k in range(crd):
		for m in range(crd):
			jkms[(crd**2)*j+crd*k+m]=[j,2*k,m]

#finds the best value of kappa for the given jkm values
best=binary_search(lambda val:find_E(val,jkms)[0,0])

#prints out optimal E
print("best estemate of the energy eigen value is:",best[1],"eV")

#prints out the relative error
print("which has a relative error compared to experimental values of",np.abs((best[1]+79.005168)/79.005168))