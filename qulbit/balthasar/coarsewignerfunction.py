#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# coarsewignerfunction.py: A class for Coarse-grained Wigner functions
#
# © 2016 Olivia Di Matteo (odimatte@uwaterloo.ca)
#
# This file is part of the project Balthasar.
# Licensed under BSD-3-Clause
#

from itertools import product 

import numpy as np
#from pynitefields import *

from balthasar.wignerfunction import WignerFunction

class CoarseWignerFunction(WignerFunction):
    """ Class to store and plot a coarse-grained discrete Wigner function.

        CoarseWignerFunction inherits many properties from WignerFunction, 
        such as the MUBs, field, dimensions, and the original kernel.

        Args:
            wf (WignerFunction): The Wigner function to coarse-grain.
            coarse_field (GaloisField): A field over which to perform 
                                        coarse-graining. 

        Keyword Args:
            mode (string): There are two modes for coarse-graining, "general",
                           or "subfield". The default mode is "general"; for
                           square dimensions it is also possible to coset using
                           the subfield.
            basis (list of FieldElement): A user-specified basis of the field
                                          in which to perform coarse-graining.
                                          If unspecified, default is the 
                                          polynomial basis. Only valid for the
                                          "general" mode.
            cosets (list of lists): A user-specified partitioning of the 
                                    finite field to use as the cosets for
                                    coarse-graining. Use with caution, as no
                                    perfect checking mechanisms are implemeneted.
            
        Attributes:
            subfield (list): The effective subfield within the large field
            subfield_map (list): Map between elements of the large field with those
                                 in the coarse field
            cosets (list of lists): The set of cosets of the big finite field
            coarse_D (dictionary): Surviving displacement operators, indexed by 
                                   points in the coarse-grained space.
            coarse_kernel (dictionary): The kernel operators for the coarse-
                                        grained Wigner function.
    """

    def __init__(self, wf, coarse_field, **kwargs):
        # Initialize all the standard parameters
        WignerFunction.__init__(self, wf.mubs)

        self.coarse_field = coarse_field

        # Loop through the optional arguments 
        self.basis = []
        self.cosets = []
        self.mode = "general"

        # We will probably need to do some better validity checks here but this
        # should be okay for now (hopefully...)
        if kwargs is not None:
            # Choose the basis with which to do coarse-graining (usual this 
            # will matter only in the general case.
            if "basis" in kwargs:
                self.basis = kwargs["basis"]
            else:
                self.basis = self.compute_polynomial_basis()

            if "mode" in kwargs:
                self.mode = kwargs["mode"]
                if self.mode != "general" and self.mode != "subfield":
                    print("Error, the mode must be 'general' or 'subfield'.")
                    return None
                # Dimension needs to be square for subfield mode. In other
                # words, the number of elements in each coset must be the same
                # as the size of each coset.
                if self.mode == "subfield":
                    if self.coarse_field.dim ** 2 != self.field.dim:
                        print("Error, dim must be square for subfield mode.")
                        return None

            if "cosets" in kwargs:
                # Make sure the cosets are valid
                field_as_list = [self.field[i] for i in range(self.dim)]
                cosets_flattened = [el for c in kwargs["cosets"] for el in c]
                if sorted(field_as_list) != sorted(cosets_flattened):
                    print("Error, are not a proper partitioning of the field.")
                    return None
                # If they're valid, choose representatives and carry on
                self.cosets = kwargs["cosets"]

        # Compute the subfield
        self.subfield, self.subfield_map = self.compute_subfield()
        if self.subfield == None:
            return None

        # If the cosets were not provided, compute them now. 
        if len(self.cosets) == 0:
            self.cosets = self.compute_cosets()

        # Compute the coarse displacement operators.
        self.coarse_D = self.compute_coarse_D()

        if self.mubs.matrices == False:
            print("Warning: coarse Wigner function constructed without matrices.")
            return

        # Using the coarse displacement operators, compute the coarse kernel
        # The points here are indexed by their coset representatives.
        self.coarse_kernel = self.compute_coarse_kernel()
        

    def compute_subfield(self):
        """ Compute the subfield from the big field. e.g. if our original
            system is dimension 16, and we are coarse-graining down to 
            dimension 4, we want to find the copy of F4 in F16.
        """
        subfield = []
        subfield_map = {}

        # For now let's concern ourselves only with the square case.
        # More specifically, the square case where we are considering an
        # extension over a non-prime field.
        # Let the original field have dimension d. Then the subfield has
        # is generated by powers of the element f[sqrt(d) + 1]
        if self.coarse_field.dim ** 2 == self.field.dim and self.coarse_field.n > 1:
            prim_element = self.field[self.coarse_field.dim + 1]
            subfield = [prim_element ** x for x in range(self.coarse_field.dim)]    

            # Make the subfield map
            for i in range(len(subfield)):
                subfield_map[subfield[i]] =  self.coarse_field[i]
        else:
            # If the dimension is not square, or field we are coarse-graining
            # to is prime, we can only pull out a copy of
            # what Luis affectionately calls the "mother field", i.e. the
            # base prime field. These are the elements whos exp_coefs are all
            # 0 except for the first one. Find these!

            # We need the field in the polynomial basis to do this properly,
            # because in the self-dual basis the elements are not right form.
            field_poly_basis = GaloisField(self.field.p, self.field.n, \
                self.field.coefs)

            for el in field_poly_basis: 
                if all(a == 0 for a in el.exp_coefs[1:]):
                    el_in_sdb = self.field[el.prim_power] # "Real" element
                    subfield.append(el_in_sdb)
                    subfield_map[el_in_sdb] = self.coarse_field[el.exp_coefs[0]]
                    #subfield_map[el_in_sdb] = self.coarse_field[el.prim_power]

        return subfield, subfield_map


    def compute_polynomial_basis(self):
        """ Return the polynomial basis for coarse graining. 

            When the dimension isn't square, this is just the polynomial basis 
            in the 'big' field. When it is square, we choose the polynomial 
            basis of the big field with respect to the small field. 
            For example, in dimension 16 coarse-graining to dimension 4, we 
            choose the basis :math:`\{1, \\sigma\}` rather than the full poly 
            basis because dim 16 is a 2-dim vector space over the dim 4 case.
        """
        if self.coarse_field.dim ** 2 == self.field.dim: # Square dimensions
            return [self.field[-1], self.field[1]]
        else: # Standard polynomial basis
            return [self.field[-1]] + [self.field[x] for x in range(1, self.field.n)] 


    def compute_cosets(self):
        """ Coset the large finite field over the coarse-grained one. 
            
            Two cases to consider here, one where there is a basis provided 
            to make the cosets, the other where we use the subfield.

            Details can be found in our manuscript.
        """
        cosets = []

        if self.mode == "general":
            # In the general mode, for a system p^mn, we can really only take
            # out one copy of p in the most general case by using this method.
            # However, we can also apply the general mode to square dimensions
            # By choosing a 'small' basis and use the subfield elements
            # as coset representatives.
            first_coset = []

            # First generate all combinations of the coefficients with the
            # subfield except the first element (which is 1 in most cases).
            coefficient_lists = product(self.subfield, repeat = len(self.basis) - 1) 
            for c in coefficient_lists:
                l_comb = [c[i-1] * self.basis[i] for i in range(1, len(c) + 1)]
                s = self.field[0]
                for el in l_comb:
                    s += el
                first_coset.append(s)

            for i in range(len(self.subfield)):
                next_coset = [self.subfield[i] * self.basis[0]
                    + first_coset[t] for t in range(len(first_coset))]
                cosets.append(next_coset)

        else:
            l_om = self.field[1] # Primitive element of original field
            b_om = self.subfield[1] # Primitive element of subfield

            for c_idx in range(len(self.subfield)):
                if c_idx == 0:
                    cosets.append(self.subfield)
                elif c_idx == 1: # Special case because of how we defined pow
                    next_coset = [x + l_om for x in self.subfield]
                    cosets.append(next_coset)
                else:
                    next_coset = [x + l_om * (b_om ** (c_idx - 1)) \
                        for x in self.subfield]
                    cosets.append(next_coset)

        return cosets


    def compute_coarse_D(self):
        """ Compute the coarse-grained displacement operators. 

            This will be a subset of the fine-grained displacement operators
            which 'survive' the sum of eq. x in our paper.
        """
        survivors = []

        for alpha in self.field:
            # Using the gchar here allows us to consider qubits AND qudits
            l = []
            if self.field.p == 2:
                l = [gchar(self.cosets[0][i] * alpha) for i in range(len(self.cosets[0]))]
            else:
                l = [gchar(self.cosets[0][i] * alpha).eval() for i in range(len(self.cosets[0]))]
            if not np.isclose([sum(l)], [0])[0]:
                survivors.append(alpha.prim_power)

        # Collect the surviving operators into a table
        # Note that the MUB table does not contain identity operators
        # at the beginning of each row - thus, we will have to take only 
        # the non-zero survivors, and change the index by -1.
        surviving_ops = {}
        for slope in self.subfield:
            table_row = self.mubs.table[slope.prim_power]
            surviving_ops[slope.prim_power] = \
                [table_row[x-1][0] for x in survivors[1:]]

        # Infinite slope case
        infty_row = self.mubs.table[-1]
        surviving_ops["inf"] = [infty_row[x-1][0] for x in survivors[1:]]

        return surviving_ops 


    def compute_coarse_kernel(self):
        """ Compute the kernel of the coarse-grained Wigner function. 

            This is done by 'globbing' together the fine-grained kernel coset
            by coset: 

            .. math::

              \Delta_c(C_i, C_j) = \sum_{\\alpha \in C_i} \sum_{\\beta \in C_j} \Delta(\\alpha, \\beta)

            The final kernel will be indexed by points in the subfield here, 
            though really it should be indexed by coset representatives.

        """

        ckernel = {}

        # Coarse kernel should be indexed by cosets / subfield elements
        for alpha in range(len(self.subfield)):
            for beta in range(len(self.subfield)):
                coarse_point = (self.coarse_field[alpha], \
                    self.coarse_field[beta])

                mat = np.zeros((self.field.dim, self.field.dim), \
                    dtype = np.complex_)
                
                # Sum up the kernel points from the appropriate cosets
                for x in self.cosets[alpha]:
                    for y in self.cosets[beta]:
                        fine_point = (x, y)
                        mat = mat + self.kernel[fine_point]

                ckernel[coarse_point] = mat
                
        return ckernel
                


    def compute_wf(self, state):                                                
        """ Compute the coarse Wigner function for a given state.

            The idea there is the same as for the normal Wigner function, 
            except the matrix will have the dimensions of the coarse field.

            Args:                                                               
                state (np array/matrix): The state to compute the Wigner Function
                                         of. This is a numpy array and can      
                                         either be a vector, or a full density  
                                         matrix.                                

            Returns:
                A numpy matrix which is the coarse Wigner function of the state.
        """                                                                     

        # For the coarse Wigner function the dimension is that of the 
        # underlying affine plane.
        W = np.zeros((self.coarse_field.dim, self.coarse_field.dim)) 

        # Turn kets into density operators if need be.
        if state.shape[0] == 1:                                                 
            state = np.outer(state, np.conj(state))                             
          
        # A sorted copy of the subfield for plotting
        sorted_els = sorted(self.coarse_field.elements)

        # The coarse Wigner function is indexed by the subfield, so use this.
        for alpha in self.subfield: 
            for beta in self.subfield: 
                coarse_point = (self.subfield_map[alpha], self.subfield_map[beta])
                mat = np.trace(np.dot(state, self.coarse_kernel[coarse_point]))
                
                # Determine where in the Wigner matrix to put this value
                a = sorted_els.index(coarse_point[0])
                b = sorted_els.index(coarse_point[1])

                W[a][b] = (1.0 / self.field.dim) * mat 

        return W   
