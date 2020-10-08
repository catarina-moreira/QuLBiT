import pyAgrum as gum
import pyAgrum.lib.notebook as gnb
import numpy as np

from qutip import *

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm

import seaborn as sns
sns.set()

# plots the quantum bloch sphere associated to a node in a Bayesian Network
def visualiseQuantumState( bn,  var, experiment_title, path, col = "b", save=False, size = (3, 3) ):

	vec = np.sqrt(bn.cpt(bn.idFromName(var)).tolist())
	
	b = Bloch()
	b.fig = plt.figure( figsize=size )

	if( len( vec.shape) == 1):
		plot_single_state( b, vec, col[0], 0)
	else:
		plot_double_state( b, vec[0,:], vec[1,:], col, 0 )

	if save:
		b.save(path + "bloch/BLOCH_" + var + "_" + experiment_title + ".png")

	return b 

# plots the quantum bloch sphere associated to a node in a Bayesian Network
def plot_single_state( b, vec, color, n ):
	
	psi1 = [ vec[0]*(np.cos(th1) + 1j*np.sin(th1))*basis(2,0) for th1 in np.linspace(0, 2*np.pi, 20) ]
	psi2 = [ vec[1]*(np.cos(th2) + 1j*np.sin(th2))*basis(2,1) for th2 in np.linspace(0, 2*np.pi, 20) ]

	colors = [color for x in range(20)]
	b.vector_color = colors
	for i in range(0, len(psi1)):
		for j in range(0, len(psi2)):
			b.add_states(psi1[i] + psi2[j])
	return b

# plots the quantum bloch sphere associated to a node in a Bayesian Network
def plot_double_state( b, vec1, vec2, color, n ):

	psi11 = [ vec1[0]*(np.cos(th1) + 1j*np.sin(th1))*basis(2,0) for th1 in np.linspace(0, 2*np.pi, 20) ]
	psi12 = [ vec1[1]*(np.cos(th2) + 1j*np.sin(th2))*basis(2,1) for th2 in np.linspace(0, 2*np.pi, 20) ]

	psi21 = [ vec2[0]*(np.cos(th1) + 1j*np.sin(th1))*basis(2,0) for th1 in np.linspace(0, 2*np.pi, 20) ]
	psi22 = [ vec2[1]*(np.cos(th2) + 1j*np.sin(th2))*basis(2,1) for th2 in np.linspace(0, 2*np.pi, 20) ]
	
	colors = color*20
	b.vector_color = colors
	 
	for i in range(0, len(psi11)):
		for j in range(0, len(psi12)):
			b.add_states( psi11[i] + psi12[j])  
			b.add_states( psi21[i] + psi22[j]) 
	return b

# plots the quantum bloch sphere associated to a node in a Bayesian Network
def visualise_wigner_function( rho, my_method, X_MAX, path, experiment_title, isPartial, color_map1 = cm.RdBu ):

	xvec = np.linspace(- X_MAX, X_MAX, 300)
	W_density = wigner(rho,xvec,xvec, method = my_method )

	color_map2 = wigner_cmap( W_density )  # Generate Wigner colormap

	nrm = mpl.colors.Normalize(-4, 4)

	fig, axes = plt.subplots(1, 2, figsize=(12, 4), dpi=300)

	plt1 = axes[0].contourf(xvec, xvec, W_density, 100, cmap=color_map1)
	axes[0].set_title(experiment_title );
	cb1 = fig.colorbar(plt1, ax=axes[0])

	plt2 = axes[1].contourf(xvec, xvec, W_density, 100, cmap=color_map2) # Apply Wigner colormap  cm.RdBu
	axes[1].set_title(experiment_title );
	cb2 = fig.colorbar(plt2, ax=axes[1])

	fig.tight_layout()
	if isPartial:
		fig.savefig(path + "wigner/WPARTIAL_TRACE_"+ experiment_title + ".png")
	else:
		fig.savefig(path + "wigner/WDENSITY_"+ experiment_title + ".png")
	plt.show()
	

