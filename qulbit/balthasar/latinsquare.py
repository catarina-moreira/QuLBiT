#!/usr/bin/python                                                                  
# -*- coding: utf-8 -*-                                                            
#                                                                                  
# latinsquare.py: A class for making Latin squares. 
#                                                                                  
# © 2016 Olivia Di Matteo (odimatte@uwaterloo.ca)                                  
#                                                                                  
# This file is part of the project Balthasar.                                      
# Licensed under BSD-3-Clause                                                      
# 

import numpy as np

#from pynitefields import *

class LatinSquare():
    """ Class to hold a Latin square.

        A Latin square L is computed from a curve :math:`c` according to:

        .. math::
                    L(\\alpha, \\beta) = \\alpha + c(\\beta)
      
        Args:
            curve (Curve): The curve with which to create the Latin square.

        Attributes:
            p (int): Prime dimension
            n (int): Power of prime 
            dim (int): The dimension of the system :math:`p^n`
            curve (Curve): The curve with which this square was made 
            square (array): The actual square (stored as a numpy array)
    """

    def __init__(self, curve):
        # Set some obvious parameters
        self.dim = self.field.dim
        self.curve = curve

        self.square = np.zeros((self.dim, self.dim), dtype=np.int)

        for row in range(self.dim):
            for column in range(self.dim):
                row_el = curve[row][1]
                col_el = curve[column][0]
                self.square[row][column] = (row_el + col_el).prim_power


    def print(self):
        """ Prints the Latin square.
        """
        print(self.square)
