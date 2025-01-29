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

    self.N_mats=np.transpose(mats)
    
  def make_H_mats(self):
    #generate H matrices from H_func
    print("Constructing the H matrices.")

    mats=np.empty((self.dim,self.dim,self.num))

    for i in range(self.dim):
      for j in range(self.dim):
        mats[i,j]=self.H_func(i,j,*self.params)

    self.H_mats=np.transpose(mats)

  def find_N_eigens(self):
    print("Finding the eigenvectors and eigenvalues of the N matrices.")
    self.N_eigenvalues, self.N_eigenvectors=np.linalg.eig(self.N_mats)
  
  def make_P_mats(self):

    print("Constructing the P matrices.")
    
    #calculating invs sqrt beta
    invs_sqrt_beta=np.zeros((self.dim,self.dim,self.num))
    for i in range(self.dim):
      invs_sqrt_beta[i,i]=1/np.sqrt(self.N_eigenvalues[:,i])
    invs_sqrt_beta=np.transpose(invs_sqrt_beta)

    #calculating Y
    Y=self.N_eigenvectors
    Y_T=np.transpose(Y,axes=[0,2,1])

    #construct the P matrix
    self.P_mats=invs_sqrt_beta@Y_T@self.H_mats@Y@invs_sqrt_beta

  def find_P_eigens(self):

    print("Finding P eigenvectors and eigenvalues.")

    self.P_eigenvalues, self.P_eigenvectors=np.linalg.eig(self.P_mats)

  def find_energy_levels(self):
    
    print("Calculating the energy levels.")

    energy_levels=np.sort(self.P_eigenvalues,axis=1)

    self.energy_levels=energy_levels
    
    

