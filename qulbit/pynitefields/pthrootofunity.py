#!/usr/bin/python                                                                  
# -*- coding: utf-8 -*-                                                            
#                                                                                  
# pthrootofunity.py: Implementation of pth roots of unity. 
#                                                                                  
# © 2016 Olivia Di Matteo (odimatte@uwaterloo.ca)                                  
#                                                                                  
# This file is part of the project PyniteFields.                                      
# Licensed under BSD-3-Clause                                                      
# 

import math

import numpy as np

class pthRootOfUnity():
    """ Class to hold :math:`p^{\\text{th}}` roots of unity symbolically over finite fields. 
        
        These guys have the form

        .. math::

            \omega_p = \exp \left( \\frac{2 \\pi i}{p} \\right)

        where :math:`p` is the characteristic of the field.
        They can be evaluated both symbolically and also exactly by 
        explicitly computing the value above using numpy. 
        Here we'll implement only what we need: exponents and multiplication.

        Args:
            p (int): The characteristic of some finite field.
            e (int): An exponent, if desired (default is 1).

        An object of type pthRootOfUnity has the following attributes.

        Attributes:
            p (int): Which root of unity we are dealing with.
            e (int): An additional exponent representing a power, i.e. 
                     :math:`\omega_p ^ e`.

    """

    def __init__(self, p, e = 1):
        # TODO some sort of primality check
        if p < 2:
            print ("Error, please enter a prime to make roots of unity.")
        else:
            self.p = p
            self.e = e


    def __mul__(self, op):
        """ Multiply two roots of unity. Roots of unity are cyclic, i.e.
            :math:`\omega_p ^ p = 1`.

            Args:
                op (pthRootOfUnity): What to multiply by, say, 
                                     :math:`\omega_p^{e^\prime}`.

            Returns: 
                :math:`\omega_p^{e} \cdot \omega_p^{e^\prime}`.
        """
        if self.p != op.p:
            print("Error, cannot multiply roots of unity from different primes.")
            return 

        new_exp = (self.e + op.e) % self.p
        return pthRootOfUnity(self.p, new_exp)  


    def __imul__(self, op):
        """ Multiplication with assignment. """
        return self * op


    def __truediv__(self, op):
        """ Division. 

            Args:
                op (pthRootOfUnity): What to divide by, say, 
                                     :math:`\omega_p^{e^\prime}`.

            Returns: 
                :math:`\omega_p^{e} /  \omega_p^{e^\prime}`.
        """
        if self.p != op.p:
            print("Error, cannot divide roots of unity from different primes.")
            return 

        new_exp = (self.e - op.e) % self.p
        return pthRootOfUnity(self.p, new_exp)  


    def __itruediv__(self, op):
        """ Division with assignment. """
        return self / op


    def __pow__(self, ex):
        """ Exponentiation. 
        
            Args:
                ex (int): Some integer exponent.

            Returns:
                :math:`\omega_p^{e \cdot ex}`, where :math:`ex` is
                the exponent passed in as an argument.
        """
        new_exp = (self.e * ex) % self.p
        return pthRootOfUnity(self.p, new_exp)  

            
    def __eq__(self, op):
        """ Equality of two roots of unity.

            We consider two roots of unity equal if both their characteristics
            :math:`p` and exponents :math:`e` are equal.
        
            Args:
                op (pthRootOfUnity): A pth root to check equality with.

            Returns:
                True if the primes and exponents are the same; false otherwise.
                None if there is a type error.
        """
        if type(op) != pthRootOfUnity:
            print("Error, type cannot be compared with pthRootOfUnity.")
            return None

        if self.p != op.p:
            return False
        if self.e != op.e:
            return False
        return True


    def __repr__(self):
        """ Print pth roots of unity in the command line.

            We represent the :math:`\omega` using a simple 'w'.
        """ 
        return ("w^" + str(self.e))


    def eval(self):
        """ Use numpty to evaluate the actual number. Gross. 
        
            Returns: 
                The numerical value :math:`\exp \\left(\\frac{2 \\pi i \cdot e}{p} \\right)`.
        """
        return pow(np.exp(2 * np.pi * 1j / self.p), self.e)


    def print(self):
        """ Prints the pth root of unity. 
        """
        print("w^" + str(self.e))
