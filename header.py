from mpmath import mp
import numpy as np

class Subspace:

  def __init__(self, dim , verbose=False):

    #the dimension of the subspace
    self.dim = dim

    self.verbose=verbose

  def set_N_func(self,f):
    #a funciton that given two indecies and some other paramiters returns the specified element of the N ie the inner product
    self.N_func=f

  def set_H_func(self,f):
    #a funciton that given two indecies and some other paramiters returns the specified element of the projecitno of the hamiltonian
    self.H_func=f

  def set_params(self,params):
    #the paramiters that are used to calculate H and N matrices
    self.params=params

  def make_N_mat(self):
    #generate N matrices from N_func
    if self.verbose:
      print("Constructing the N matrices.")
    indecies=np.array(np.meshgrid(np.arange(self.dim),np.arange(self.dim),indexing='ij')).transpose([1,2,0])
    #first calculate matrix in numpy way then convert tp mpmath matrix, this is mostly to avoid rewriting code
    npmatrix=self.N_func(indecies[:,:,0],indecies[:,:,1],*self.params)
    self.N_mat=mp.matrix(npmatrix.tolist())

  def make_H_mat(self):
    #generate H matrices from H_func
    if self.verbose:
      print("Constructing the H matrices.")
    indecies=np.array(np.meshgrid(np.arange(self.dim),np.arange(self.dim),indexing='ij')).transpose([1,2,0])
    npmatrix=self.H_func(indecies[:,:,0],indecies[:,:,1],*self.params)
    self.H_mat=mp.matrix(npmatrix.tolist())

  def find_N_eigens(self):
    if self.verbose:
      print("Finding the eigenvectors and eigenvalues of the N matrices.")
    self.N_eigenvalues, self.N_eigenvectors=mp.eig(self.N_mat)
  
  def make_invs_sqrt_beta_mats(self):
    if self.verbose:
      print("Constructing the inverse square root beta matrices.")

    self.invs_sqrt_beta_mat=mp.diag((1/np.sqrt(self.N_eigenvalues)).tolist())
    
  def make_Y_mat(self):
    if self.verbose:
      print("Constructing the Y matrices.")

    self.Y=self.N_eigenvectors

  def make_P_mats(self):
    if self.verbose:
      print("Constructing the P matrices.")

    self.P_mats=self.invs_sqrt_beta_mat*self.Y.T*self.H_mat*self.Y*self.invs_sqrt_beta_mat

  def find_P_eigens(self):
    if self.verbose:
      print("Finding P eigenvectors and eigenvalues.")

    self.P_eigenvalues, self.P_eigenvectors=mp.eig(self.P_mats)

  def find_energy_levels(self):
    if self.verbose:
      print("Calculating the energy levels.")

    energy_levels=np.sort(self.P_eigenvalues)

    self.energy_levels=energy_levels

  def find_energy_eigenstates(self):
    if self.verbose:
      print("Calculating the components of the energy eigenstates.")

    order=np.argsort(self.P_eigenvalues)

    unorderd_energy_eigenstates=self.Y*self.invs_sqrt_beta_mat*self.P_eigenvectors

    energy_eigenstates=np.empty((self.dim,self.dim),dtype=mp.mpf)

    for i in range(self.dim):
        energy_eigenstates[i]=unorderd_energy_eigenstates[:,order[i]]
    self.energy_eigenstates=energy_eigenstates

    
    
    

