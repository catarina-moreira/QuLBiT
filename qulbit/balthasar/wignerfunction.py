#!/usr/bin/python                                                                  
# -*- coding: utf-8 -*-                                                            
#                                                                                  
# wignerfunction.py: A class implementing standard discrete Wigner functions. 
#                                                                                  
# © 2016 Olivia Di Matteo (odimatte@uwaterloo.ca)                                  
#                                                                                  
# This file is part of the project Balthasar.                                      
# Licensed under BSD-3-Clause                                                      
# 

import numpy as np
#from pynitefields import *

class WignerFunction():
    """ Class to store and plot a discrete Wigner function.

        Args:
            mubs (MUBs): A set of MUBs. The displacement operators will
                         be used to construct the Wigner function kernel.
                         If the 'matrices' option in MUBs is False, no kernel
                         will be produced.

        Attributes:
            field (GaloisField): The finite field over which the Wigner 
                                 function is defined
            dim (int): The dimension of the system 
            mubs (MUBs): The MUBs used to construct this Wigner function
            D (dictionary): The dictionary of displacement operators
            kernel (dictionary): Operators at each point in discrete phase 
                                 space ('point operators' according to Wootters)
    """

    def __init__(self, mubs):
        self.field = mubs.field
        self.dim = self.field.dim
        self.mubs = mubs

        if self.mubs.matrices == False:
            print("Warning: Wigner Function constructed without matrices.")
            print("This Wigner Function cannot be used for plotting or \
                any numerical purposes, aside for computing coarse-grained \
                displacement operators.")
            return

        self.D = mubs.D

        # Compute the kernel
        self.kernel = self.compute_kernel()    

         
    def compute_kernel(self):
        """ Compute the 'kernel' of the Wigner function, i.e. the set of 
            what Wootters calls point operators.

            The kernel is a set of operators associated to each point in 
            phase space, :math:`\Delta(\\alpha, \\beta)`. The operator at the
            origin is computed as the sum of all the displacement operators:

            .. math::

                \Delta(0, 0) = \\frac{1}{p^n} \sum_{\\alpha, \\beta} D(\\alpha, \\beta)

            The kernel operator at any other point can be computed by
            translating that of the origin by the appropriate displacement
            operator:

            .. math ::

                \Delta(\\alpha, \\beta) = D(\\alpha, \\beta) \Delta(0, 0) D^\dagger (\\alpha, \\beta)

            Returns:
                A dictionary containing the mapping from points :math:`(\\alpha, \\beta)`
                to the operator :math:`\Delta(\\alpha, \\beta)`.
        """

        kernel = {} 

        if self.mubs.matrices == False:
            print("Error, no matrices, cannot compute Wigner Function kernel.")
            return

        # Computation of the kernel can be accomplished by computing the
        # value at point (0, 0), then translating it using the D operators,
        # i.e.        w(a, b) = D(a, b) w(0, 0) D(a, b)^\dag
          
        kernel_00 = np.zeros((self.dim, self.dim), dtype=np.complex_)
        if self.field.p != 2: # Special case, need numerical values of pthrootofunity
            for key in self.D.keys():
                kernel_00 = kernel_00 + (self.D[key][0].eval() * self.D[key][1])
            kernel_00 = kernel_00 / self.dim
            kernel[(self.field[0], self.field[0])] = kernel_00
        else:
            for key in self.D.keys():
                kernel_00 = kernel_00 + (self.D[key][0] * self.D[key][1])
            kernel_00 = kernel_00 / self.dim
            kernel[(self.field[0], self.field[0])] = kernel_00

        # Compute the rest of the points by translating the first one
        for a in self.field:
            for b in self.field:
                if a == self.field[0] and b == self.field[0]:
                    continue # Don't set the 0 case again
                # Compute the displacement operator with phase included
                if self.field.p == 2: # For p = 2 this is a number
                    dab = self.D[(a, b)][0] * self.D[(a, b)][1]
                else: # For p != 2 we need to evaluate the power of pth root 
                    dab = self.D[(a, b)][0].eval() * self.D[(a, b)][1]

                dab_dag = np.asmatrix(dab).getH()

                kernel[(a, b)] = np.dot(dab, np.dot(kernel_00, dab_dag))

        return kernel

                            
                                
    def compute_wf(self, state):
        """ Compute the Wigner function for a given state.

            For a state :math:`\\rho`, the value of the Wigner function at
            a point :math:`(\\alpha, \\beta)` is given by

            .. math::

                W_\\rho (\\alpha, \\beta) = \\text{Tr}(\\rho \Delta(\\alpha, \\beta))

            Args:
                state (np array/matrix): The state to compute the Wigner Function
                                         of. This is a numpy array and can 
                                         either be a vector, or a full density
                                         matrix.

            Returns:
                A numpy matrix which is the Wigner function of the passed state.
        """

        if self.mubs.matrices == False:
            print("Error, no matrices, cannot compute Wigner Function.")
            return

        W = np.zeros((self.dim, self.dim)) # Holds result

        # Don't discriminate - allow the user to submit either a ket vector
        # or a density operator. If it's a ket, switch to a density operator.
        if state.shape[0] == 1:
            state = np.outer(state, np.conj(state))

        # We sort the elements of the finite field so that Wigner function
        # Gets plotted in the "correct" order, i.e. the computational basis
        # states go 000, 001, 010, etc. and same for the +/- basis. This in 
        # a way imposes an 'order' on the field elements.
        sorted_els = sorted(self.field.elements)

        for a in self.field:
            for b in self.field:
                mat = np.trace(np.dot(state, self.kernel[(a, b)]))
                a_coord = sorted_els.index(a)
                b_coord = sorted_els.index(b)
                W[a_coord][b_coord] = (1.0 / self.dim) * mat

        return W

    
    def plot_mat(self, state):
        """ A simple matrix plot of the Wigner function. 
        
            Args:
                state (np array/matrix): The state to compute the Wigner Function
                                         of. This is a numpy array and can 
                                         either be a vector, or a full density
                                         matrix.
        """
        import matplotlib.pyplot as plt
        W = self.compute_wf(state)
        plt.matshow(W)
        plt.show()


    
    def plot(self, state, filename=""):
        """ Compute and plot the Wigner function of a given state. 

            The plot is a 3D representation of the phase space with the axes
            labeled by kets in the computational and +/- basis.

            For the CoarseWignerFunction class, the axes are labeled by groups
            of kets in the same coset.

            Args:
                state (np array/matrix): The state to compute the Wigner Function
                                         of. This is a numpy array and can 
                                         either be a vector, or a full density
                                         matrix.
                filename (string): If a filename is passed, the Wigner function
                                   plot will not be displayed on screen, but
                                   instead saved as an eps file.
        """
    
        if self.mubs.matrices == False:
            print("Error, no matrices, cannot plot Wigner function.")
            return

        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt

        W = self.compute_wf(state)

        # This code was written in the summer and frankly I still don't 
        # fully understand, or at this point remember, why it works >.<
        fig = plt.figure()

        ax1 = fig.add_subplot(211, projection='3d')

        xpos = []
        ypos = []
        zpos = []

        for x in range(0, len(W)):
            for y in range(0, len(W)):
                xpos.append(x)
                ypos.append(y)
                zpos.append(0)

        dx = np.ones(len(xpos))
        dy = np.ones(len(xpos))
        dz = W.flatten()

        max_colour = dz.max()
        min_colour = dz.min()
        colours = []

        for point in dz:
            if point < 0:
                colours.append(( 1 - (point / min_colour), 1 - (point / min_colour), 1))
            elif point == 0:
                colours.append((1, 1, 1))
            else:
                colours.append((1, 1 - (point / max_colour), 1 - (point / max_colour)))

        ax1.bar3d(xpos, ypos, zpos, dx, dy, dz, color = colours)

        plt.gca().invert_xaxis()

        # For standard Wigner functions, set the ticks to be basis states.
        # When the Wigner function gets created it is ordered in this way.
        if type(self) == WignerFunction:
            fmt_str = "#0" + str(self.field.n + 2) + "b"
            comp_basis = [format(x, fmt_str)[2:] for x in range(self.field.dim)]
            pm_basis = [x.replace('0', '+').replace('1', '\u2013') for x in comp_basis]

            plt.xticks(range(0, len(W)))
            plt.yticks(range(0, len(W)))
            ax1.set_xticklabels(pm_basis)
            ax1.set_yticklabels(comp_basis)
        else:
            # Make labels corresponding to the cosets
            comp_basis = []
            pm_basis = []
            for coset in self.cosets:
                c_labs = ["".join([str(x) for x in el.sdb_coefs]) for el in coset]
                pm = [b.replace('0', '+').replace('1', '\u2013') for b in c_labs]

                
                next_string = ""
                for c in range(len(c_labs)):
                    next_string += (c_labs[c] + " ")
                    if (c + 1) % 2 == 0 and (c > 0):
                        next_string += "\n"
                comp_basis.append(next_string)

                next_string = ""
                for c in range(len(pm)):
                    next_string += (pm[c] + " ")
                    if (c + 1) % 2 == 0 and (c > 0):
                        next_string += "\n" 
                pm_basis.append(next_string)
                       

            plt.xticks(range(0, len(W)))
            plt.yticks(range(0, len(W)))
            ax1.set_xticklabels(pm_basis)
            ax1.set_yticklabels(comp_basis)
#        plt.gca().set_xticks([])
#        plt.gca().set_yticks([])
        
        if dz.min() < 0:
            plt.gca().set_zlim([dz.min() - 0.01, dz.max() + 0.01])
        else:
            plt.gca().set_zlim([0, dz.max() + 0.01])

        plt.tick_params(axis='x', labelsize=7)
        plt.tick_params(axis='y', labelsize=7)
        plt.tick_params(axis='z', labelsize=16)
        ax1.locator_params(nbins=6, axis='z')

        plt.tight_layout()

        if filename == "":
            plt.show()
        else:
            plt.savefig(filename, dpi=1200, bbox_inches='tight')
