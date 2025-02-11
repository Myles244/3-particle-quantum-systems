import numpy as np

class Subspace:

  def __init__(self, dim):

    #the dimension of the subspace
    self.dim = dim

  def set_N_func(self,f):
    #a funciton that given two indecies and some other paramiters returns the specified element of the N ie the inner product
    self.N_func=f

  def set_H_func(self,f):
    #a funciton that given two indecies and some other paramiters returns the specified element of the projecitno of the hamiltonian
    self.H_func=f

  def set_params(self,num,params):
    #the paramiters that are used to calculate H and N matrices
    self.num=num
    self.params=params

  def make_N_mats(self):
    #generate N matrices from N_func
    print("Constructing the N matrices.")

    mats=np.empty((self.dim,self.dim,self.num))

    for i in range(self.dim):
      for j in range(self.dim):
        mats[i,j]=self.N_func(i,j,*self.params)

    self.N_mats=np.transpose(mats,[2,0,1])
    
  def make_H_mats(self):
    #generate H matrices from H_func
    print("Constructing the H matrices.")

    mats=np.empty((self.dim,self.dim,self.num))

    for i in range(self.dim):
      for j in range(self.dim):
        mats[i,j]=self.H_func(i,j,*self.params)

    self.H_mats=np.transpose(mats,[2,0,1])

  def make_N_mats_vectorized(self):
    #generate N matrices from N_func
    print("Constructing the N matrices.")
    indecies=np.array(np.meshgrid(np.arange(self.dim),np.arange(self.dim),indexing='ij')).transpose([1,2,0])
    self.N_mats=np.transpose(self.N_func(indecies[:,:,0],indecies[:,:,1],*self.params),[2,0,1])

  def make_H_mats_vectorized(self):
    #generate H matrices from H_func
    print("Constructing the H matrices.")
    indecies=np.array(np.meshgrid(np.arange(self.dim),np.arange(self.dim),indexing='ij')).transpose([1,2,0])
    self.H_mats=np.transpose(self.H_func(indecies[:,:,0],indecies[:,:,1],*self.params),[2,0,1])

  def find_N_eigens(self):
    print("Finding the eigenvectors and eigenvalues of the N matrices.")
    self.N_eigenvalues, self.N_eigenvectors=np.linalg.eigh(self.N_mats)
  
  def make_invs_sqrt_beta_mats(self):

    print("Constructing the inverse square root beta matrices.")

    invs_sqrt_beta=np.zeros((self.num,self.dim,self.dim))
    for i in range(self.dim):
      invs_sqrt_beta[:,i,i]=1/np.sqrt(self.N_eigenvalues[:,i])

    self.invs_sqrt_beta_mats=invs_sqrt_beta
    
  def make_Y_mats(self):

    print("Constructing the Y matrices.")

    self.Y=self.N_eigenvectors

  def make_P_mats(self):

    print("Constructing the P matrices.")

    Y_T=np.transpose(self.Y, axes=[0,2,1])

    self.P_mats=self.invs_sqrt_beta_mats@Y_T@self.H_mats@self.Y@self.invs_sqrt_beta_mats

  def find_P_eigens(self):

    print("Finding P eigenvectors and eigenvalues.")

    self.P_eigenvalues, self.P_eigenvectors=np.linalg.eigh(self.P_mats)

  def find_energy_levels(self):
    
    print("Calculating the energy levels.")

    energy_levels=np.sort(self.P_eigenvalues,axis=1)

    self.energy_levels=energy_levels

  def find_energy_eigenstates(self):

    print("Calculating the components of the energy eigenstates.")

    order=np.argsort(self.P_eigenvalues,axis=1)

    unorderd_energy_eigenstates=self.Y@self.invs_sqrt_beta_mats@self.P_eigenvectors

    energy_eigenstates=np.zeros((self.num,self.dim,self.dim))

    for i in range(self.num):
      for j in range(self.dim):
        energy_eigenstates[i,j]=unorderd_energy_eigenstates[i,:,order[i,j]]/np.linalg.norm(unorderd_energy_eigenstates[i,:,order[i,j]]) #normalise vectors
        energy_eigenstates[i,j]=energy_eigenstates[i,j]*np.sign(energy_eigenstates[i,j,0]) #set first element as posetive to avoid sign flipping

    self.energy_eigenstates=energy_eigenstates

    
    
    

