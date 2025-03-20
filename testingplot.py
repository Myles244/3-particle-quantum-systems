import matplotlib.pyplot as plt
import numpy as np
from mpmath import mp
from mpl_toolkits.mplot3d import axes3d
from matplotlib import colormaps
from exponential import *

#load in data
params=np.load("data/best100params.npy",allow_pickle=True)
params=np.append(params,np.load("data/params.npy",allow_pickle=True),axis=1)
params=np.append(params,np.load("data/altparams.npy",allow_pickle=True),axis=1)
params=np.append(params,np.load("data/altaltparams.npy",allow_pickle=True),axis=1)

subspace=Subspace(params.shape[1],verbose=True)

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
expdelta=delta_r23(subspace.energy_eigenstates[0],params)

print("\nthe expectation of the delta23:",expdelta)

As=np.log(np.float64(np.abs(subspace.energy_eigenstates[0])**2))

xs,ys,zs=np.log(np.float64(params))

#sort out colours
cmap=colormaps["plasma"]
colours=cmap((zs-np.min(zs))/(np.max(zs)-np.min(zs)),(As-np.min(As))/(np.max(As)-np.min(As)))

#creat 3d figure
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

#plot data
ax.scatter(xs, ys, zs,color=colours)

# Set the axis labels
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

# Rotate the axes and update
while True:
    for angle in range(0, 360*4 + 1):
        # Normalize the angle to the range [-180, 180] for display
        angle_norm = (angle + 180) % 360 - 180

        # Cycle through a full rotation of elevation, then azimuth, roll, and all
        elev = azim = roll = 0
        if angle <= 360:
            elev = angle_norm
        elif angle <= 360*2:
            azim = angle_norm
        elif angle <= 360*3:
            roll = angle_norm
        else:
            elev = azim = roll = angle_norm

        # Update the axis view and title
        ax.view_init(elev, azim, roll)
        plt.title('Elevation: %d°, Azimuth: %d°, Roll: %d°' % (elev, azim, roll))

        plt.draw()
        plt.pause(.001) 