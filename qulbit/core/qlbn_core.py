import itertools as it
import numpy as np

import pyAgrum as gum
import pyAgrum.lib.notebook as gnb

from qutip import Qobj
from qutip import tensor

# Computes the conjugate tranpose of a quantum state
def conj_transpose( quantum_state ):
	"""
	quantum_state: a quantum state object
	"""
	return Qobj(np.conj(np.transpose(quantum_state)))
   

# Computes the full joint probability distribution of the Bayesian network
def computeFullJoint( bn, var="" ):
	"""
	bn: a classical Bayesian network representing a decision problem
	var: a string containing the variable name to write the distribution, Pr( var | X1, ...)
	"""
	classical_joint = bn.cpt(0)
	for node_id in range(1, len(bn.nodes())):
		classical_joint = classical_joint*bn.cpt(node_id)
	
	if len(var) != 0:
		classical_joint = classical_joint.putFirst(var)
	
	return classical_joint

# Takes a vector of real numbers and converts them to complex amplitudes
def convert_to_complex_amplitude( lst_real_nums ):
	"""
	lst_real_nums: a list of real numbers
	"""
	conjTrans = np.sqrt(np.array(lst_real_nums))
	return conjTrans

# Creates the quantum superposition state corresponding to a 
# classical full joint distribution of a Bayesian Network
def createSuperpositionState( joint ):
	"""
	joint: a multiarray containing the classical full joint distribution of a BN
	"""
	N = len(joint.toarray().shape)
	joint_flat = joint.tolist() 								# convert multiarray to a list of lists
	for i in range(N - 1):
		joint_flat = list(it.chain.from_iterable(joint_flat))	# flatten the list of lists
	
	joint_array = convert_to_complex_amplitude(joint_flat)	# convert to complex amplitudes
	superposition = Qobj(joint_array) 						# convert to a quantum state
	return superposition

# Sets the values of the evidence vars to zero in the full joint probability distribution
def observeEvidence( joint, evidence_vars ):
	"""
	joint: a multiarray containing the classical full joint distribution of a BN
	evidence_vars: a dictionary containing paris of {var_name : value}
	"""		
	return joint
	
# Creates a density matrix, which is the quantum version of a full joint probability distribution
# from a quantum superposition state
def computeDensityMatrix( superposition ):
	"""
	superposition: a quantum object representing the quantum superposition state of the BN
	"""	
	rho = tensor(superposition, conj_transpose(superposition))
	return rho

# Takes partial trace over the subsystem defined by 'axis', which corresponds to the 
# marginal probability distribution
# the main diagonal contains the marginal probabilities of a classical distribution
def computePartialTrace(rho, dims, axis=0):
	"""
	rho: a density matrix representing the full joint of the network
	dims: a list containing the dimension of each subsystem
	axis: the index of the subsystem to be traced out
	"""
	dims_ = np.array(dims)
	# Reshape the matrix into a tensor with the following shape:
	# [dim_0, dim_1, ..., dim_n, dim_0, dim_1, ..., dim_n]
	# Each subsystem gets one index for its row and another one for its column
	reshaped_rho = rho.reshape(np.concatenate((dims_, dims_), axis=None))

	# Move the subsystems to be traced towards the end
	reshaped_rho = np.moveaxis(reshaped_rho, axis, -1)
	reshaped_rho = np.moveaxis(reshaped_rho, len(dims)+axis-1, -1)

	# Trace over the very last row and column indices
	traced_out_rho = np.trace(reshaped_rho, axis1=-2, axis2=-1)

	# traced_out_rho is still in the shape of a tensor
	# Reshape back to a matrix
	dims_untraced = np.delete(dims_, axis)
	rho_dim = np.prod(dims_untraced)
	return Qobj(traced_out_rho.reshape([rho_dim, rho_dim]))

# performs quantum-like probabilistic inferences out of a classical Bayesian Network
def performQuantumLikeInference( bn, query_var, evidence_vars ):
	"""
	bn: a classical Bayesian network representing a decision problem
	query_var: a string containing the variable name to query
	evidence_vars: a dictionary with pairs { var_name : var_val } with observed evidences
	"""
	joint = computeFullJoint( bn, var=query_var )  			 # compute classical full joint distribution
	joint = observeEvidence( joint, evidence_vars )			 # observe evidence
	superposition = createSuperpositionState( joint )		 # compute superposition state
	density_matrix = computeDensityMatrix( superposition )	 # compute the correspondent density matrix
	trace = computePartialTrace(rho, dims, axis= bn.idFromName(query_var))	 # compute the partial trace over the query var
	
	
	return density_matrix, trace

# checks if the density operator is correct by checking if the main diagonal (classical full joint)
# sums to 1	
def checkFullJoint( rho ):
	"""
	rho: a density matrix representing the full joint of the network
	"""
	classical_joint = np.diag( rho )
	isJoint = np.round(sum(classical_joint),8) == 1.0 
	
	return isJoint 

	
	
	
	



