from header import *



x1, x2, y1, y2 = 0.8, 1.5, -80, -70 
subkappas=np.linspace(0.8,1.5,10000)
jkms=[[0,0,0],[0,0,1],[0,2,0]]

subE=find_E(subkappas,jkms)

kappas=np.linspace(0.5,20,10000)

E=find_E(kappas,jkms)

#find optimal k
kappas_order=np.argsort(subE[:,0])
optimal_kappa=subkappas[kappas_order[0]]
optimal_energy_eigenvalues=subE[kappas_order[0]]
#plt.plot([optimal_kappa,optimal_kappa],[-100,150],'--')

#plt.scatter(kappas,sorted_test_energy_eigenvalues[:,1])
#plt.plot(kappas,test_E0(kappas))

print("Optimal value for kappa:",optimal_kappa)
print("Coresponing energy eigenvalues:",optimal_energy_eigenvalues)

def phi(jkl,kappa,r1,r2):
	r1_norm=np.linalg.norm(r1,axis=-1)
	r2_norm=np.linalg.norm(r2,axis=-1)
	r12_norm=np.linalg.norm(r2-r1,axis=-1)
	return (
		(r1_norm+r2_norm)**jkl[0]
		*(r1_norm+r2_norm)**jkl[1]
		*(r12_norm)**jkl[2]
		*np.exp((-Z/(kappa*r0))*(r1_norm+r2_norm))
		)





fig, ax = plt.subplots()

ax.set_yticks([0,-20,-40,-60,-80])
ax.set_xticks([0,5,10,15,20])

ax.plot(kappas,E[:,0])

# inset Axes....
 # subregion of the original image
axins = ax.inset_axes(
    [0.5, 0.1, 0.4, 0.4],
    xlim=(x1, x2), ylim=(y1, y2), xticks=[optimal_kappa], yticks=[optimal_energy_eigenvalues[0]])
axins.set_yticklabels(["$-79.0\,eV$"])
axins.set_xticklabels(["$1.10$"])
axins.plot(subkappas,subE[:,0])
axins.scatter([optimal_kappa],[optimal_energy_eigenvalues[0]],marker='o')
axins.plot([0,optimal_kappa],[optimal_energy_eigenvalues[0],optimal_energy_eigenvalues[0]],color='C0',linestyle='--')
axins.plot([optimal_kappa,optimal_kappa],[-80,optimal_energy_eigenvalues[0]],color='C0',linestyle='--')

ax.indicate_inset((x1,y1,x2-x1,y2-y1), edgecolor="black")
ax.add_patch(ConnectionPatch(xyA=(x2, y1), xyB=(0, 0), axesA=ax, axesB=axins, coordsA="data", coordsB="axes fraction", lw=0.7,ls="--"))
ax.add_patch(ConnectionPatch(xyA=(x2, y2), xyB=(0, 1), axesA=ax, axesB=axins, coordsA="data", coordsB="axes fraction", lw=0.7,ls="--"))


ax.set_xlabel("$\kappa$")
ax.set_ylabel("$E/eV$")

plt.savefig("variation.png",dpi=300)


#the ijm value associated with each n value
crd=3 #cube root dim
pjkmv=np.arange(0,crd) # posible jkm values
jkms=np.zeros((crd**3,3),dtype=int)
for j in range(crd):
	for k in range(crd):
		for m in range(crd):
			jkms[(crd**2)*j+crd*k+m]=[j,2*k,m]

#print(jkms)

tempf=np.linalg.norm(jkms*[1,0.5,1],axis=-1)
#print(tempf)
jkms_order=np.argsort(tempf)
#print(jkms_order)
jkms=jkms[jkms_order]
#print(jkms)

kappas=np.linspace(0.2,4,10000)
jkmss=[[[0,0,0],[0,0,1],[0,2,0]],jkms[:4],jkms[:9],jkms[:81]]

Es=np.zeros((4,kappas.size))
for r in range(4):
	Es[r]=find_E(kappas,jkmss[r])[:,0]

	
fig, ax = plt.subplots()

ax.set_ylim((-85,5))
ax.set_xlim((0,4))
for r in range(4):
	ax.plot(kappas,Es[r])

ax.set_yticks([0,-20,-40,-60,-80])
ax.set_xticks([0,1,2,3])


ax.set_xlabel("$\kappa$")
ax.set_ylabel("$E/eV$")

axins = ax.inset_axes(
    [0.5, 0.5, 0.4, 0.4],
    xlim=(x1, x2), ylim=(y1, y2), xticks=[optimal_kappa], yticks=[optimal_energy_eigenvalues[0]])
axins.set_yticklabels(["$-79.0\,eV$"])
axins.set_xticklabels(["$1.10$"])

subEs=np.zeros((4,subkappas.size))
for r in range(4):
	subEs[r]=find_E(subkappas,jkmss[r])[:,0]

for r in range(4):
	axins.plot(subkappas,subEs[r])
#axins.scatter([optimal_kappa],[optimal_energy_eigenvalues[0]],marker='o')
#axins.plot([0,optimal_kappa],[optimal_energy_eigenvalues[0],optimal_energy_eigenvalues[0]],color='C0',linestyle='--')
#axins.plot([optimal_kappa,optimal_kappa],[-80,optimal_energy_eigenvalues[0]],color='C0',linestyle='--')



kappas_order=np.argsort(subEs[-1])
optimal_kappa=subkappas[kappas_order[0]]
optimal_energy_eigenvalues=subEs[-1,kappas_order[0]]

best_kappa=binary_search(lambda val:find_E(val,jkmss[-1])[0,0])

axins.scatter(best_kappa[0],best_kappa[1],marker='o',color='C3')
#axins.plot([0,optimal_kappa],[optimal_energy_eigenvalues,optimal_energy_eigenvalues],color='C0',linestyle='--')
#axins.plot([optimal_kappa,optimal_kappa],[-80,optimal_energy_eigenvalues],color='C0',linestyle='--')

print("Optimal value for kappa:",optimal_kappa)
print("Coresponing energy eigenvalues:",optimal_energy_eigenvalues)

print("best_kappa",best_kappa)

ax.indicate_inset((x1,y1,x2-x1,y2-y1), edgecolor="black")
ax.add_patch(ConnectionPatch(xyA=(x1, y2), xyB=(0, 1), axesA=ax, axesB=axins, coordsA="data", coordsB="axes fraction", lw=0.7,ls="--"))
ax.add_patch(ConnectionPatch(xyA=(x2, y2), xyB=(0, 0), axesA=ax, axesB=axins, coordsA="data", coordsB="axes fraction", lw=0.7,ls="--"))



fig.savefig("convergence.png",dpi=300)

print(binary_search(lambda val:find_E(val,jkmss[0])[0,0]))