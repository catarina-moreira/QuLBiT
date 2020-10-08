#!/usr/bin/python                                                                  
# -*- coding: utf-8 -*-                                                            
#                                                                                  
# mubs.py: A class implementing mutually unbiased bases. 
#                                                                                  
# © 2016 Olivia Di Matteo (odimatte@uwaterloo.ca)                                  
#                                                                                  
# This file is part of the project Balthasar.                                      
# Licensed under BSD-3-Clause                                                      
# 

from functools import reduce
from operator import mul

import numpy as np
#from pynitefields import *

from curve import Curve 
from striations import Striations

class MUBs():
    """ Class to hold a complete set of mutually unbiased bases (MUBs).

        Args:
            f (GaloisField): The finite field over which the MUBs will be 
                             constructed. Must be expressed in self-dual basis.

        Class MUBs also takes the following keyword arguments.

        Keyword Args:
            curves (list): Advanced functionality. Pass a full set of curves 
                           with which to create the MUB table. By default the
                           set of linear curves is used.
            matrix (bool): Tells Balthasar whether to construct the full
                           numerical matrices for the MUB operators. True
                           by default.

        Class MUBs contains the following public attributes.

        Attributes:
            field (GaloisField): The finite field of choice 
            p (int): A prime number, the dimension of a single particle
            n (int): The number of particles
            dim (int): The dimension of the system dim = :math:`p^n`
            w (pthRootOfUnity): The :math:`p^{th}` root of unity in the field.
            curves (list of Curves): The set of curves used to construct the table
            table: The table of operators, in form (name, phase, matrix)
            D (dictionary): Table of displacement operators, mapping of points
                            in phase space to the associated (phase, matrix)
            matrices (bool): Default true. Constructs the matrices associated
                            with the operators. If set to false, only the 
                            string representations of the operators is computed.

    """

    def __init__(self, f, **kwargs):
        if f.n > 1 and f.is_sdb is False:
            print("Warning - you have passed a finite field whose \
                    expansion coefficients are not represented in \
                    the self-dual basis.")
            return None

        # Set the attributes that depend on the field
        self.field = f # Keep a copy of the finite field to do math!
        self.p = f.p 
        self.n = f.n
        self.dim = f.dim
        self.w = f.w

        # Initialize MUB specific attributes
        self.curves = [] # Which curves to use
        self.table = [] # Operator table in operator form
        self.D = {} # Dictionary of displacement operators
        self.matrices = True # Compute matrix reps of operators 

        # The following code is useful only for odd-prime powers.
        # Find and store the inverse of two only once, to later compute phases.
        self.twoinv = None
        if self.p != 2: # Easy to find for prime dimensions
            if self.n == 1:
                self.twoinv = self.field[2].inv()
            else: # Use polynomial basis expansion to find 2
                for el in self.field: 
                    if el.exp_coefs == ([2] + ([0] * (self.n -1))): 
                        self.twoinv = self.field[el.inv().prim_power]

        # Deal with the keyword arguments
        # Set the curves, default to Desarguesian bundle if nothing passed in 
        if "curves" in kwargs: 
            if self.verify_curves(curves) == True:
                self.curves = curves 
        else: # Desarguesian curves
            self.curves = Striations.generate_rays(f)

        # Determine whether to compute the matrix or not.
        if "matrix" in kwargs:
            if kwargs["matrix"] == False:
                self.matrices = False   

        # Build the operator table in matrix form at the same time
        self.table, self.D = self.build_operator_table()



    def phi(self, a, b):
        """ Phase function for the displacement operators.

            Arguments:
                a (FieldElement): Horizontal coordinate in phase-space
                b (FieldElement): Vertical coordinate in phase-space

            Returns:
                The value of :math:`\Phi(a, b)`.

            Currently, the phase functions are expressed as follows:

            .. math::

              \Phi(\\alpha, \\beta) = (-1)^{f(\\alpha \\beta)}i^{\\text{tr}(\\alpha \\beta)}

            where :math:`f(x)` is the companion function f_m which determines
            the sign of the displacement operators (so that we never end up 
            with the negative identity).
        """

        if self.p == 2: # Qubits
            if self.n == 1: # Single qubit case, +/- i^ab
                return 1j ** (a * b).prim_power
            else: 
                poly_res = self.f_m(self.n, a * b)

                if poly_res == None:
                    print("Error evaluating phase polynomial.")
                    print("Result was not 0 or 1.")
                    return

                return ((-1) ** poly_res) * (1j ** tr(a * b))
        else: # Qudits
            if self.n == 1: # Single qudit case, w ^ (2^-1 ab)
                prefactor = (-1) * (self.twoinv * a * b).prim_power
                return pow(self.w, prefactor) 
            else: # Multiple qudits
                what_to_trace = (self.field[0] - self.field[-1]) * self.twoinv * a * b
                return gchar(what_to_trace)


    def f_m(self, m, x):
        """ A phase factor which ensures our displacement operators
            sum to proper projectors for qubit systems.

            This portion of the phase is a polynomial over field elements
            We let :math:`f_0(x) = f_1(x) = 0`. Then we define

            .. math::

                f_m (x) = \sum_{i=1}^{m-1} \sum_{j=0}^{i-1} x^{2^i + 2^j} 

            Args:
                m (int): The number of qubits.
                x (FieldElement): The argument of the function.

            Returns:
                The value of :math:`f_m(x)` as defined above. This will be a
                field element and is either 0 or 1. As the output of the 
                polynomial should be an integer rather than a FieldElement, the
                phase function phi contains some code to make this conversion.
        """

        # Base case f_0 (x) = f_1(x) = 0
        if m == 0 or m == 1:
            return 0
        else:
            res = self.field[0]

            # Compute the value using the closed form of the polynomial
            for i in range(1, m):
                for j in range(i):
                    res += pow(x, pow(2, i) + pow(2, j))

            if res == self.field[0]:
                return 0
            elif res == self.field[-1]:
                return 1
            else:
                return None



    def build_operator_table(self):
        """ Construct the table of MUB operators.

        MUBs are commonly represented in two forms: collections of mutually
        unbiased basis vectors, or a table of disjoint, commuting operators 
        whose sets of mutual eigenvectors are mutually unbiased bases.
        We choose the latter representation, and represent MUBs as a table
        of tuples of the form (name, phase, matrix); the organization of this
        table will be specified after what follows.
        
        An operator in the MUB table is represented as a displacement 
        operator of the form 

        .. math::

         D(\\alpha, \\beta) = \Phi(\\alpha, \\beta) Z_\\alpha X_\\beta

        where :math:`\\alpha, \\beta` are field elements, and 
        :math:`\Phi(\\alpha, \\beta)` is a phase factor.
        The two operators Z and X are defined as having the action
        
        .. math::
 
              Z_\\alpha |\ell \\rangle = \omega(\\alpha \ell) |\ell \\rangle, \quad
              X_\\beta |\ell \\rangle = |\ell + \\beta \\rangle.
              
        For a single qubit these operators are exactly equal to the Pauli
        operators; for multiple qubits, tensor products thereof. Similarly 
        for qudits, they are the generalized Paulis and tensor products thereof.
        These operators satisfy the Heisenberg-Weyl commutation relations, i.e.

        .. math::
              Z_\\alpha X_\\beta = \chi(\\alpha \\beta) X_\\beta Z_\\alpha

        where 

        .. math::
              
             \chi(\\alpha) = \omega ^ {\\text{tr}(\\alpha)}

        is the character of the field, and :math:`\omega` is the pth root of 
        unity, :math:`\exp \\left(\\frac{2 \\pi i}{p} \\right)`.
        
        The trace of a field element is 

        .. math::
 
            \\text{tr}(\\alpha) = \\alpha + \\alpha^p + \cdots + \\alpha^{p^{n - 1}},

        and always yields an element of the mother (base prime) field.

        We can expand :math:`\\alpha` and :math:`\\beta` in terms of a self-
        dual basis :math:`\{\\theta_1, \ldots, \\theta_n\}`:

        .. math::

            \\alpha = \sum_{i=1}^n (a_i \\theta_i), \quad  
            \\beta = \sum_{i=1}^n (b_i \\theta_i),

        where :math:`a_i = \\text{tr}(\\alpha \\theta_i)` and similarly for the 
        :math:`b_i`. With this, we can express the :math:`Z_\\alpha X_\\beta` 
        in tensor product form:

        .. math::

            Z_\\alpha X_\\beta = Z^{a_1} X^{b_1} \otimes \cdots \otimes Z^{a_n} X^{b_n}


        In a way, we can consider this as assigning each self-dual basis 
        element to a single particle.

        Displacement operators will be stored in a dictionary of tuples,
        indexed by the coordinates in phase space :math:`(\\alpha, \\beta)`.
        The values have the form

         (name, :math:`\Phi(\\alpha, \\beta)`, :math:`Z_\\alpha X_\\beta`)
       
        Here name will be something like "Z ZX X", or a string representing
        the tensor product structure above. The phase factor is included 
        separately from the matrix because they're not always needed.

        
        By default, we will produce MUBs of the Desarguesian form, 
        where :math:`\\beta = \lambda \\alpha`.

        """

        table = []
        D = {}

        # Hold the generalized Paulis and the identity  
        Z = np.zeros((self.p, self.p), dtype=np.complex_)
        X = np.zeros((self.p, self.p), dtype=np.complex_)
        I = np.identity(self.p)
        
        if self.p == 2: # Simple case of qubit Paulis
            X = np.array([[0, 1], [1, 0]]) # X
            Z = np.array([[1, 0], [0, -1]]) # Z
        else:
            # X, thanks to 
            # SO questions/10936767/rearranging-matrix-elements-with-numpy
            perm_order = [self.p - 1] + [x for x in range(self.p - 1)] 
            X = I[perm_order, :]

            # Diagonal generalized Z
            powers_of_w = [pow(self.w, i).eval() for i in range(self.p)]
            np.fill_diagonal(Z, powers_of_w)

        # Now it's time to actually build the tuples of the operator table
        for curve in self.curves:
            row = [] # Each curve produces a different row of the table
            for point in curve: # (a, b)
                op = []
                op_mats = []

                a, b = point[0], point[1]

                if a == self.field[0] and b == self.field[0]:
                    continue # We ignore the identity

                phase = self.phi(a, b)

                # Get the expansion coefficients in the sdb if power of prime.
                if self.n > 1:
                    z = a.sdb_coefs 
                    x = b.sdb_coefs
                else:
                    z = a.exp_coefs
                    x = b.exp_coefs

                for idx in range(len(x)):
                    # The Z portion must be treated separately because of almost sdb
                    z_exp = z[idx]

                    # Get the inverses of the SDB norms -> Need the prime field
                    prime_f = GaloisField(self.p)

                    if self.n > 1:
                        z_exp = (z_exp * self.field.sdb_norms[idx]) % self.p

                    if z_exp == 0 and x[idx] == 0: # Both coefs 0
                        op.append("I") # Tensor factor is identity
                    elif z_exp == 0 and x[idx] != 0:
                        op.append("X" + ("" if x[idx] == 1 else str(x[idx])))
                    elif z_exp != 0 and x[idx] == 0:
                        op.append("Z" + ("" if z_exp == 1 else str(z_exp)))
                    else:
                        op.append("Z" + ("" if z_exp == 1 else str(z_exp)) + \
                            "X" + ("" if x[idx] == 1 else str(x[idx])))

                    if self.matrices:
                        # Matrix for this chunk of the tensor product
                        # Recall that for qudit systems we need to take into account the
                        # almost self-dual extra coefficient on the first Z. 
                        Z_part = np.linalg.matrix_power(Z, z_exp) 
                        X_part = np.linalg.matrix_power(X, x[idx])
                        op_mats.append(np.dot(Z_part, X_part))
                          
                if self.matrices:
                    matrix_op = reduce(np.kron, op_mats)
                    row.append( (op, phase, matrix_op) )
                
                    # Add the displacement operator to the matrix
                    D[(a, b)] = (phase, matrix_op)
                else:
                    row.append( (op, phase, None) )

            table.append(row) # Add rows to the table

            # Finally, add the identity operator to the D table if user wants
            if self.matrices:
                id_phase = self.phi(self.field[0], self.field[0])
                D[(self.field[0], self.field[0])] = (id_phase, np.identity(self.dim)) 

        return table, D


    def verify_curves(self, curves):
        """ Check that the properties of user-provided curves are valid.
          
            Curves must be additive and commutative.

            TODO: everything. Currently the only thing that is checked 
            is whether there are enough curves.
        """
        if len(curves) != self.dim + 1:
            print("Error, not enough curves provided.")
            return False
        return True


    def print(self, matrix_form = False):
        """ Print a nice(-ish) formatted MUB table.

            Args:
                matrix_form (bool): If set to True, will print out the matrices
                                    as well as the operator names.
        """

        np.set_printoptions(precision=4, threshold=np.nan, suppress=True)
        for row in self.table:
            for operator in row:
                print(str(operator[1]) + " ", end = "") # Phase
                print(" ".join(operator[0]) + "\t\t", end = "")
                if matrix_form: # Print as matrices
                    print()
                    print(operator[2]) # Matrix
            print("\n")
