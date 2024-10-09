import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.axes
#for an array of values of kappa

#a) calculate matrix N

#define wavefunctions phi
Z=2
r0=0.52917#anstroms
def phi(i,j,m,kappa,r1,r2):
	r1_norm=la.norm(r1)
	r2_norm=la.norm(r2)
	return (r1_norm+r2_norm)**i*(r1_norm-r2_norm)**j*la.norm(r1-r2)**m*np.exp(-Z*(r1_norm+r2_norm)/(kappa*r0))

#test phi by plotting it

#turns a complex number into an rgb colour code where the hue represents the phase and the amplitude is the luminosity
def complex_number_to_colour(z):
	offset=2/3*np.pi
	phase=((np.angle(z)-offset)%(2*np.pi))/(2*np.pi)
	modulus=abs(z)
	hsv=np.stack((phase,modulus,np.full(np.shape(modulus),1)),axis=-1)
	return matplotlib.colors.hsv_to_rgb(hsv)

x,y=np.meshgrid(np.linspace(-5,5,1000),np.linspace(-5,5,1000))
r=np.sqrt(x**2+y**2)
#Phi=r*np.exp(-r/(2*r0))*np.exp(1j*np.angle(np.stack(x+y*1j,axis=-1)))
Phi=np.exp(-r)

fig = plt.figure(figsize=[8,6])
ax1 = fig.add_axes([0.1,0.1,0.6,0.8])
ax1.pcolorfast(complex_number_to_colour(Phi))
ax1.set_xticks([0,500,1000])
ax1.set_xticklabels(["-5A","0","5A"])
ax1.set_yticks([0,500,1000])
ax1.set_yticklabels(["-5A","0","5A"])
thetas=np.linspace(0,2*np.pi,100)
radii=np.linspace(0,1,50)
tr, rr=np.meshgrid(thetas,radii)
hsv=np.stack((((tr/(2*np.pi))-1/3)%1,rr,np.full(np.shape(tr),1)),axis=-1)
rgb=matplotlib.colors.hsv_to_rgb(hsv)

ax2 = fig.add_axes([0.75,0.1,0.2,0.2],polar='true')
ax2.pcolorfast((tr-np.pi*rr)%(2*np.pi),rr,rgb)
ax2.set_yticks([0])
ax2.set_xticks([0,np.pi/2,np.pi,3/2*np.pi])
ax2.set_xticklabels(['1','i','-1','-i'])
ax2.grid(False)
plt.savefig("ground state")



#plt.show()
#calculate norms and save in matrix N

#b) calculate matices beta sqrt_beta inverse_sqrt_beta Y Y_T

#find eigenvalues and eigen vectors of N

#construct diagonal matrix +-sqrt beta from eigen values of N

#construct Y from eigen vectors of N

#c) calculate matrix H_tilde and then P

#calculate elements of H_tilde

#calculate P using definition

#c) Find the eigenvalues of P

#for the ground state plot the energy against k