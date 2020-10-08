#!/usr/bin/python                                                                  
# -*- coding: utf-8 -*-                                                            
#                                                                                  
# curve.py: A class for curves over finite fields in discrete phase space.
#                                                                                  
# © 2016 Olivia Di Matteo (odimatte@uwaterloo.ca)                                  
#                                                                                  
# This file is part of the project Balthasar.                                      
# Licensed under BSD-3-Clause                                                      
# 

#from pynitefields import *

class Curve():
    """ Class to hold all points in a curve. 

        Curves are sets of points of the form :math:`(\\alpha, c(\\alpha)` for
        all :math:`\\alpha` in a specified GaloisField, where

        .. math::
          c(\\alpha) = c_0 + c_1 \\alpha + c_2 \\alpha^2 + ...

        They are constructed by passing in a list of coefficients in the form
        :math:`[c_0, c_1, c_2, \ldots]`.

        Args:
            coefs (list): The coefficients of the curve. 
            field (GaloisField): The field over which the curve is defined.

        Attributes:
            field (GaloisField): The finite field in which is curve is defined
            coefs (list): The coefficients :math:`[c_0, c_1, \ldots]`.
            form (string): "beta" or "alpha", tells whether to do the curve 
                           as :math:`beta = f(alpha) or alpha = f(beta).
                           By default, we use the beta form, beta = f(alpha).
            is_ray (bool): A Boolean which indicates whether the curve passes 
                           through the point (0, 0) or not
            points (list): A list of tuples of field elements which are the 
                           points of this curve over the field.
    """

    def __init__(self, coefs, field, form="beta"):
        self.field = field
        self.coefs = coefs
        self.form = form # Default curve in form beta (alpha)
        self.is_ray = False # Does it pass through 0, 0
        self.points = [] # List of points as tuples

        # Determine if the curve is a ray by checking the constant term 
        if type(coefs[0]) is int:
            if coefs[0] == 0:
                self.is_ray = True
        else:
            if coefs[0] == self.field[0]:
                self.is_ray = True

        # Compute all the points on the curve
        for el in field:
            if self.form == "beta":
                self.points.append( (el, field.evaluate(coefs, el)) )
            elif self.form == "alpha":
                self.points.append( (field.evaluate(coefs, el), el) )


    def __getitem__(self, index):
        """ Get a point on the curve.

            Args:
                index (int): The index at which we would like to
                             find the value on the curve. Expressed as a power
                             of the primitive element of the field.

            Returns:
                The tuple (index, value of the curve at index).
        """

        if index < len(self.points) and index >= 0:
            return self.points[index]
        else:
            print("Error, element out of bounds.")
        

    def __iter__(self):
        """ Allow the user to iterate over the curve point by point """
        return iter(self.points)


    def print(self, as_points = False):
        """ Print the curve.

        Args:
            as_points (bool): If True is passed, will print the list of 
                              points on the curve as tuples. By default, prints
                              the curves as a polynomial.

        """
        if as_points == True: # Print the curve as a list of points
            for point in self.points:
                print(str(point[0].prim_power) + ", " + str(point[1].prim_power), end = "\t")
        else: # Print the curve as a polynomial
            print(self.form + "(x) = ", end = "") 
            for i in range(len(self.coefs)):
                if type(self.coefs[i]) == int: # Integers
                    if self.coefs[i] == 0: # Ignore 0 terms unless the curve is 0
                        continue
                    if self.coefs[i] == 1 and i == 0: # Only print 1 if it's the constant
                        print(str(self.coefs[i]), end="")
                        continue
                    print(str(self.coefs[i]), end="")
                else: # Field elements
                    if self.coefs[i] == self.field[0]: # Ignore 0 terms unless curve is 0
                        continue
                    if (self.coefs[i].prim_power == 1):
                        if self.field.n > 1:
                            print("s", end="")
                    else:
                        if self.field.n > 1:
                            print("s" + str(self.coefs[i].prim_power), end="")
                        else:
                            print(str(self.coefs[i].prim_power), end="")

                if i > 0:
                    if i == 1:
                        print(" x", end="")
                    else:
                        print(" x" + str(i), end="")
                if i != len(self.coefs) - 1:
                    print(" + ", end = "")
        print("")

        
        
