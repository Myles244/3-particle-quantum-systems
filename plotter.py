import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.axes
import matplotlib.animation as animation
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

def plot3d(, filename="animation.gif"):

	fig = plt.figure()
	ax = fig.add_subplot(projection='3d')

	def initiate():
		# Make data
		theta = np.linspace(0, np.pi, 50)
		phi = np.linspace(0, 2*np.pi, 50)
		phitheta=np.meshgrid(phi,theta)
		print(phitheta[1])
		print(np.outer(theta,np.ones(np.size(phi))))
		radius=display_radius(phitheta)


		x = radius * np.outer(np.sin(theta),np.cos(phi))
		y = radius * np.outer(np.sin(theta),np.sin(phi))
		z = radius * np.outer(np.cos(theta),np.ones(np.size(phi)))

		colours=np.ones((np.shape(z)[0],np.shape(z)[1],3))
		colours[:,:,0]=display_phase(phitheta)
		#print(colours)
		colours=matplotlib.colors.hsv_to_rgb(colours)

		# Plot the surface
		ax.plot_surface(x, y, z,rstride=1,cstride=1,antialiased=True, facecolors = colours)

		# Set an equal aspect ratio
		ax.set_aspect('equal')

		ax.set_axis_off()

	def rotate(frame_data):
		print(frame_data)

		azim=frame_data
		elev=30*np.sin(azim*(np.pi/(270)))
		roll=0

		# Update the axis view and title
		ax.view_init(elev, azim, roll)

	# Only save last 100 frames, but run forever
	ani = animation.FuncAnimation(fig, rotate,frames=1080,init_func=initiate)

	#ani.save(filename, fps=60, dpi=100)

	plt.show()

plot3d(
		display_radius=lambda phitheta:np.abs(np.cos(2*phitheta[1])),
		display_phase=lambda phitheta:phitheta[0]/(2*np.pi)
	)