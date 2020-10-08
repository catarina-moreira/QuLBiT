#!/usr/bin/python                                                                  
# -*- coding: utf-8 -*-                                                            
#                                                                                  
# galoisfield.py: Implementation of a full finite field. 
#                                                                                  
# © 2016 Olivia Di Matteo (odimatte@uwaterloo.ca)                                  
#                                                                                  
# This file is part of the project PyniteFields.                                      
# Licensed under BSD-3-Clause                                                      
# 

import sys
import math

from fieldelement import FieldElement
from pthrootofunity import pthRootOfUnity

class GaloisField():
    """ A finite field, or Galois field.

        Args:
            p (int): A prime number, the base of the field.
            n (int): An exponent representing the degree of a field extension.
                     By default n = 1 (prime field).
            coefs (list): A list of integers representing the coefficients of
                          an irreducible primitive polynomial of degree n over
                          GF(p). Default is the empty list for prime fields.

        Attributes:
            p (int): The prime dimension of the field
            n (int): The degree of the field extension 
            dim (int): The full order of the field, :math:`p^n`.
            w (pthRootOfUnity): The :math:`p^{\\text{th}}` root of unity.

            coefs (list): The coefficients of the irreducible polynomial
            elements (list): A list of all FieldElements in this finite field.
            is_sdb (bool): A boolean which tells us whether the elements'
                           expansion coefficients are in the self-dual
                           basis (True) or the polynomial basis (False). 
                           The default is False.
    """
    def __init__(self, p, n = 1, coefs = []):
        # TODO implement check for prime number
        self.p = p

        # Field extension parameter. 
        # Base prime field if n = 1, otherwise the field is GF(p^n)
        if n > 1:
            self.n = n
        elif n != 1:
            print("Error, invalid exponent for field extension.")
            print("Please enter an exponent n of 2 or greater.\n")
            sys.exit()
        else:
            self.n = 1
            self.coefs = []
        
        # Set separate parameter for the field dimension
        self.dim = int(math.pow(p, n))

        # Initialize the pth root of unity
        self.w = pthRootOfUnity(p)

        # Initialize the coefficients for the irreducible polynomial 
        # to do the field extension with. 
        if len(coefs) > 0:
            # We should have n + 1 coefficients for GF(p^n)
            if len(coefs) == n + 1:
                self.coefs = coefs
            else:
                print("Error, invalid number of coefficients in the irreducible polynomial.")
                print("Field of size " + str(self.p) + "^" + str(self.n) + " should have ", end = "")
                print(str(n + 1) + " coefficients in its irreducible polynomial.")
                sys.exit()


        # Generate the actual field elements
        if self.n == 1:
            # Prime case is easy. No field basis, just the numbers from 0 to p,
            # stored as FieldElements.
            self.elements = []
            for i in range(0, p):
                self.elements.append(FieldElement(self.p, self.n, [i]))
        else:
            # Use the irreducible polynomial to generate the field elements
            # They'll be stored in order as a list of coefficients in the polynomial basis
            # e.g. in dim 4, x^2 + x + 1 is the polynomial, use the basis (1, x) and store
            # the elements as:
            # 0 -> [0, 0], 1 -> [1, 0], x -> [1, 0], x^2 = [1, 1]

            # Hold all the coefficients for each element
            # For simplicity, rather than a list of list, represent each field element as a 
            # string of coefficients, i.e. [0, 1, 1] -> "011"  
            field_list = []

            # The polynomial basis contains n elements
            # The first element is always 0
            self.elements = []
            self.elements.append(FieldElement(self.p, self.n, [0]*self.n))
            field_list.append("0," * (self.n - 1) + "0")

            # The next few elements are initial terms in the poly basis (i.e. x, x^2 ...)
            for i in range(1, self.n):
                next_coefs = [0]*(i) + [1] + [0]*(self.n - i - 1) 
                self.elements.append(FieldElement(self.p, self.n, next_coefs))
                field_list.append(",".join([str(x) for x in next_coefs]))

            # For the n^th power of x, we need to use the irreducible polynomial
            nth_coefs = [((-1) * self.coefs[i]) % self.p for i in range(0, self.n)]
            self.elements.append(FieldElement(self.p, self.n, nth_coefs))
            field_list.append(",".join([str(x) for x in nth_coefs]))

            # For the remaining powers, multiply the previous element by primitive element
            for el in range(self.n + 1, self.dim):
                # Shift all coefficients ahead by 1 power of x and take the sum because
                # we know all the previous elements, and will never get anything 
                # with such a high exponent we don't know it's basis coefficients
                next_coefs = [0] + self.elements[el - 1].exp_coefs
                
                # Get a list of the powers whose coefficients aren't 0
                which_to_sum = [self.elements[i] * co for i, co in enumerate(next_coefs) if co != 0]
                sum = self.elements[0]

                for sum_el in which_to_sum:
                    sum = sum + sum_el

                # TODO Make sure that this element is not already in the list - if it is, then
                # we did not use a true primitive polynomial.
                str_rep = ",".join([str(x) for x in sum.exp_coefs])
                if str_rep not in field_list:
                    self.elements.append(sum)
                    field_list.append(str_rep)
                else:
                    raise ValueError("Repeated field element detected; please make sure your irreducible polynomial is primitive.")
                 
            # This is really dumb, but make sure each element holds a copy of the whole
            # list of the field elements. This makes field multiplication way easier.
            for i in range(len(self.elements)):
                (self.elements[i]).field_list = field_list 
                (self.elements[i]).prim_power = i

        # SDB information
        self.is_sdb = False # Have we indicated an sdb?
        self.sdb = [] # The indices of the elements that make up the sdb
        self.sdb_norms = [] # The trace of the sdb squared - usually 1, but
                            # if the sdb is almost sd, then one is not 1.


    def __getitem__(self, idx):
        """ Access specific elements in the finite field.

            Args:
                idx (int): The index of the element to retrieve. For primes 
                    this is the same as the number itself; for power-of-primes
                    it represents the power of the primitive element.
          
            Returns:
                The element at the specified index in the field. 

              None if idx is out of bounds.
        """
        if idx < self.dim and idx >= (-1 * self.dim):
            return self.elements[idx]
        else:
            print("Error, element out of bounds.")


    def __iter__(self):
        """ Make the finite field iterable. 

            Returns:
                An iterator to the field elements.
        """
        return iter(self.elements)


    def to_sdb(self, sdb_element_indices):
        """ Transform the expansions coefficients to the self-dual basis.

            Currently valid only for fields whose orders are powers of 2.

            Args:
                sdb_element_indices (list): The indices of the FieldElements 
                    (as powers of the primitive element) that represent the
                    self-dual basis. e.g. if the self-dual basis is 
                    :math:`\{ \sigma^3, \sigma^5, \sigma^6 \}`, this list
                    would be [3, 5, 6].

            TODO:
                Implement automatic construction of some self-dual basis.
        """

        if self.n == 1:
            print("Cannot take self-dual basis of a prime field.")
            return

        # Make sure that the provided sdb is valid. In qudit cases, we may
        # also be shuffling the elements, so make sure to get the shuffled copy.
        valid_sdb, valid_element_indices, valid_sdb_norms = self.verify_sdb(sdb_element_indices)

        if not valid_sdb:
            print("Invalid self-dual basis provided.")
            return

        if valid_element_indices != sdb_element_indices:
            print("The order of your self-dual basis elements has changed.")
            print("This is due to the presence of a non-1 normalization coefficient.")
            print("New ordering is " + str(valid_element_indices) + ".")

        # Set the sdb 
        self.is_sdb = True
        self.sdb = valid_element_indices
        self.sdb_norms = valid_sdb_norms

        first_norm_inverse = 1  

        if self.p % 2 == 1 and self.n % 2 == 0: # If p odd, n even
            for i in range(self.p):
                if (self.sdb_norms[0] * i) % self.p == 1:
                    first_norm_inverse = i

        # If all goes well, we can start computing the coefficients
        # in terms of the new elements by using the trace and multiplication
        # functions.
        sdb_els = [self.elements[self.sdb[i]] for i in range(0, self.n)]
        sdb_field_list = []
        for element in self.elements:
            sdb_coefs = [] # Expansion coefficients in the sdb

            for i in range(len(sdb_els)):
                if i == 0:
                    sdb_coefs.append((first_norm_inverse * tr(element * sdb_els[i])) % self.p)
                else:
                    sdb_coefs.append(tr(element * sdb_els[i]))


            sdb_field_list.append(",".join([str(x) for x in sdb_coefs]))

            element.is_sdb = True
            element.sdb_coefs = sdb_coefs

        # Finally, give every element a copy of the sdb coefficients
        for element in self.elements:
            element.sdb_field_list = sdb_field_list
    


    def verify_sdb(self, sdb_element_indices):
        """ Verify if a set of elements form a proper self-dual normal basis.

            For qubit systems, a proper self-dual basis always exists. There are
            two properties to check for this:
              * The trace of each basis element squared is 1.
              * The trace of each basis element multiplied by every other is 
                0 (orthogonality). 

            In qudit cases, the first condition needs to be modified; there
            will always be some basis element where the trace of the square is *not*
            1, but rather some other element in the prime field. To this end, we will
            keep a list of the 'normalization' coefficients for the almost sdb. We will
            also reorder the almost sdb in this case so that the non-1 element is first.
 
            Returns:
                A triple containing the following values:
                - True if above conditions are satisfied, false if not. 
                - The sdb element indices, the ordering of which may change if 
                  the basis is not perfectly self-dual. None if cond 1 is false.
                - The normalizations of the sdb elements. A list of 1s if the
                  basis is perfectly self-dual, or a positive coefficient plus
                  the rest of the list 1s if almost self-dual. None if cond 1
                  is not satisfied.
        """

        if len(sdb_element_indices) != self.n:
            print("Error, incorrect number of elements in proposed basis.")
            return False, None, None

        sdb_norms = []

        if self.p == 2: # Qubit case
            for i in range(0, self.n):
                for j in range(i, self.n): # Don't double compute things
                    trace_result = tr(self.elements[sdb_element_indices[i]] * self.elements[sdb_element_indices[j]])

                    if i == j: # Same element, should have trace 1
                        if trace_result != 1:
                            return False, None, None
                    else: # Different elements should be orthogonal and have trace 0
                        if trace_result != 0:
                            return False, None, None

            # If successful, set the orthogonality coefficients to 1
            sdb_norms = [1] * self.n

        else: # Qudit case
            for i in range(0, self.n):
                for j in range(0, self.n):
                    trace_result = tr(self.elements[sdb_element_indices[i]] * self.elements[sdb_element_indices[j]])
                    
                    if i == j: # Square the element and trace it
                        # Just needs to be in the prime_field
                        if trace_result < 0 or trace_result >= self.p:
                            return False, None, None
                           
                        # Only one element can have a non-1 normalization 
                        if trace_result == 1: 
                            sdb_norms.append(trace_result)
                        else:
                            non1 = [x for x in sdb_norms if x != 1]
                            if len(non1) >= 1:
                                print("Error, more than one element has a non-one normalization coefficient.")
                                print("Self-dual basis is invalid.")
                                return False, None, None
                            else:
                                sdb_norms.append(trace_result)
                    else: # Different elements must be trace-orthogonal
                        if trace_result != 0:
                            return False, None, None

            # For power of primes, the self-dual basis **might** be real (e.g. dim 27).
            # It's possible that all normalizations are 1, so check this, and carry on if true.
            if sdb_norms.count(1) != len(sdb_norms):
                # If the sdb so far has been okay, let's reshuffle it so the element
                # with coefficient > 1 is at the beginning. I'm honestly not sure why we 
                # do this, but this is what Andrei said to do in our LS paper.
                # Get the index of the non-1 element. Thanks to 
                # http://stackoverflow.com/questions/4111412/how-do-i-get-a-list-of-indices-of-non-zero-elements-in-a-list
                non1 = [i for i, e in enumerate(sdb_norms) if e != 1][0] 

                shuffled_sdb = [sdb_element_indices[non1]] + sdb_element_indices[:non1] + \
                                    sdb_element_indices[non1 + 1:]
                shuffled_norms = [sdb_norms[non1]] + sdb_norms[:non1] + sdb_norms[non1 + 1:]

                sdb_element_indices = shuffled_sdb
                sdb_norms = shuffled_norms
                              
        # If we made it this far, we're golden.
        return True, sdb_element_indices, sdb_norms


    def compute_sdb(self):
        """ Compute a self-dual basis for this field.
       
            .. warning::
                DO NOT USE, still under development.
        """

        # Compute a short list who's trace of their square is equal to 1
        first_round = []   
        for element in self.elements:
            if tr(element * element) == 1:
                first_round.append(element)
    
        for element in first_round:
            print(element)

        second_round = []

        # Of the remaining possible elements, compute traces and see 
        # if we can find n of them which equal 0
        for i in range(0, len(first_round)):
            traces = [tr(first_round[i] * first_round[j]) for j in range(0, len(first_round))] 
            if traces.count(0) == self.n:
                second_round.append(first_round[i])
                print(traces)

        print(second_round)

        return


    def to_poly(self):
        """ Switch back to representation in the polynomial basis. 
        """
        for el in self.elements:
            el.is_sdb = False
            el.sdb_coefs = []
        self.is_sdb = False


    def evaluate(self, coefs, argument):
        """ Evaluate a function, or curve on a finite field element.

            We consider here functions of the form
            
            .. math::
              
              f(\\alpha) = c_0 + c_1 \\alpha + \cdots + c_n \\alpha^n

            This function is primarily meant for use with the Curve class in
            my Balthasar package.

            Args:
                coefs (list): A set of coefficients for the curve, i.e.
                              :math:`[c_0, c_1, \ldots, c_n]`. These should
                              be a mix of integers and FieldElements.
                argument (FieldElement): The argument to the function, i.e.
                      the :math:`\\alpha` in :math:`f(\\alpha)`.
              
            Returns:
                The value of the function of the argument, taken over the
                finite field.
        """
        result = coefs[0] * self.elements[-1] 
        for coef_idx in range(1, len(coefs)):
            result += coefs[coef_idx] * pow(argument, coef_idx)
        return result


    def print(self):
        """ Print out all the useful information about a field."""
        
        print("--- Galois field information ---")
        print("p = " + str(self.p))

        if self.n > 1:
            print("n = " + str(self.n))

            print("\nIrreducible polynomial: ", end = "")
            if self.coefs[0] != 0:
                print(str(self.coefs[0]) + " + ", end ="")

            for i in range(1, len(self.coefs)):
                if self.coefs[i] != 0:
                    # Print coefficient if it's not 1
                    if self.coefs[i] != 1:
                        print(str(self.coefs[i]), end = "") 

                    # Print exponent value
                    print("x", end = "")
                    if i != 1:
                        print("^" + str(i), end = "")

                    if i != len(self.coefs) - 1: 
                        print(" + ", end = "")

        print("\nField elements:")
        for element in self.elements:
            element.print()


def tr(x):
    """ Wrapper trace function so the user can do tr(x) or x.trace()."""
    # Make sure x is a field element
    if type(x) is not FieldElement:
        print("Error, invalid argument to function 'tr'.")
        return None
    else:
        return x.tr()


def gchar(x):
    """ Wrapper so the user can do x.gchar() or gchar(x). """
    if type(x) is not FieldElement:
        print("Error, invalid argument to function 'gchar'.")
        return None
    else:
        return x.gchar()


def inv(x):
    """ Wrapper so the user can do x.inv() or inv(x) interchangeably."""
    # Make sure x is a field element
    if type(x) is not FieldElement:
        print("Error, invalid argument to function 'inv'.")
        return None
    else:
        return x.inv()


