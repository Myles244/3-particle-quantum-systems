import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.axes
import numpy as np

#turns a complex number into an rgb colour code where the hue represents the phase and the amplitude is the luminosity
def complex_number_to_colour(z):
	offset=1/3
	phase=(np.angle(z)/(2*np.pi)-offset)%1
	modulus=abs(z)#*np.exp(-0.1*(1+np.cos(3*np.angle(z))))
	luminosity=(modulus<=1.000001)
	hsv=np.stack((phase,modulus,luminosity),axis=-1)
	rgb = matplotlib.colors.hsv_to_rgb(hsv)
	return rgb

#plot a 2d image of the complex vatues where the luminosity reperesents the aplitude (max=1) and the hue represents the phase,
def plot(complex_values):
	fig = plt.figure(figsize=[8,6])

	#plot complex values
	ax1 = fig.add_axes([0.1,0.1,0.6,0.8])
	ax1.pcolorfast(complex_number_to_colour(complex_values))
	ax1.set_xticks([])
	ax1.set_yticks([])

	#calculate colours for legend
	thetas=np.linspace(0,2*np.pi,100)
	radii=np.linspace(0,1,50)
	tr, rr=np.meshgrid(thetas,radii)
	hsv=np.stack((((tr/(2*np.pi))-1/3)%1,rr,np.full(np.shape(tr),1)),axis=-1)
	rgb=matplotlib.colors.hsv_to_rgb(hsv)

	#plot legend
	ax2 = fig.add_axes([0.75,0.1,0.2,0.2],polar='true')
	ax2.pcolorfast((tr-np.pi*rr)%(2*np.pi),rr,rgb)
	ax2.set_yticks([0])
	ax2.set_xticks([0,np.pi/2,np.pi,3/2*np.pi])
	ax2.set_xticklabels(['1','i','-1','-i'])
	ax2.grid(False)
	#plt.show()
	plt.show()

