#!/usr/bin/python                                                                  
# -*- coding: utf-8 -*-                                                            
#                                                                                  
# striations.py: A class to implement striations in the affine plane. 
#                                                                                  
# © 2016 Olivia Di Matteo (odimatte@uwaterloo.ca)                                  
#                                                                                  
# This file is part of the project Balthasar.                                      
# Licensed under BSD-3-Clause                                                      
# 

#from pynitefields import *
from curve import Curve

class Striations():
    """ Striations, sets of parallel lines in the affine plane.

        The rays in discrete phase space are the set of non-intersecting
        lines, save for the point at the origin. Together they cover all of
        phase space. Given a single one of these rays, and translating it 
        by every element in the field, also fully covers the phase space.
        This set of parallel lines is a striation, and is the basis for
        things like Latin squares, as well as the quantum nets defined by
        Wootters when constructing his Wigner functions.
    
        The Striations class is used to construct rays and their translates.
        Mostly it is here just to keep this business separate from the more
        physical aspects of Balthasar.

        Striations are essentially stored as a set of Curves over the provided
        finite field.

        Args: 
            field (GaloisField): The field / discrete phase space to work in.

        Attributes:
            field (GaloisField): The finite field we work in.
            dim (int): The dimension of the space.
            rays (list): The set of curves in this phase space that pass 
                         through the origin. The last ray in the set is the
                         vertical one; the rest can all be accessed by the
                         index which is the power of the primitive element
                         of the slope (e.g. rays of :math:`\\sigma^3` is rays[3]).
            striations (list of lists): A list containing the full set of 
                                        translates for each ray. 
    """

    @staticmethod
    def generate_rays(field):
        """ Generates the set of rays in phase space.

            Args:
                field (GaloisField): The field over which to generate the rays.

            Returns:
                The list of rays, as objects of type Curve.
        """
        rays = []

        # All normal curves
        for slope in field:
            rays.append(Curve([field[0], slope], field))

        # The vertical curve, in alpha = g(beta) form.
        rays.append(Curve([field[0], field[0]], field, "alpha"))

        return rays


    def __init__(self, field):
        self.field = field
        self.dim = field.dim
        self.striations = []

        # Striations of the form beta = lambda alpha + gamma
        # Do these first so that the slope lambda corresponds to the list element 
        for slope in field:
            new_striation = []
            for intercept in field:
                    curve = Curve([intercept, slope], self.field)
                    new_striation.append(curve)
            self.striations.append(new_striation)

        # Handle the vertical striation (infinite slope) last
        # This way it is easily accessible as striation -1
        horizontal_striation = []
        for intercept in field:
            curve = Curve([intercept, self.field[0]], self.field, True)
            horizontal_striation.append(curve)
        self.striations.append(horizontal_striation)

        # Store the rays as a set for safekeeping
        self.rays = [s[0] for s in self.striations]
        

    def __iter__(self):
        """ Allow the user to iterate through all the striations. """
        return iter(self.striations)


    def __getitem__(self, index):
        """ Grab a striation. 
        
            Args:
                index (int): The index of the striation to collect. Striation 
                            -1 corresponds to the vertical striation,
                            all other striations are indicated by their slope 
                            as a power of the primitive field element.

            Return:
                The set of striations with slope :math:`\\sigma^\\text{index}`.
        """
        if (index >= -1) and (index <= self.dim + 1):
            return self.striations[index]
        else:
            print("Error, element out of bounds.")


    def plot(self, str_idx = 0, colours = []):
        """ Plot a set of striations in discrete phase space.
        
            Args:
                str_idx (int): The index of the striation to plot. By default
                               this function will plot the rays.
                colours (list): A set of colours to use to plot. There are
                                16 basic ones implemented; if you are plotting
                                striations of a larger phase space, you will
                                need to specify more.
        """
        
        import matplotlib.pyplot as plt

        # A subset of 16 colours to use to start
        if colours == []:
            colours = ['black', 'red', 'plum', 'yellow', 'blue', 'lightgrey', 
                       'cyan', 'darkgrey', 'lime', 'darkred', 'darkgreen', 
                       'orange', 'darkblue', 'purple', 'darkorchid', 'deeppink', 
                       'chartreuse']

        for i in range(0, len(self.striations[str_idx])):
            line = self.striations[str_idx][i]
            line_as_ints = [(pt[0].prim_power, pt[1].prim_power) for pt in line]
            plt.scatter(*zip(*line_as_ints), color = colours[i], marker = "s", s =70)
        plt.show()


    def print(self, as_points = False):
        """ Print out all the striations.
       
            Args:
                as_points (bool): If set to true, will print the curves
                                  as a set of points. Otherwise, default
                                  behaviour is to print the curves as sets of
                                  polynomials over the field elements.
        """

        print("============================")
        for s in self.striations:
            print("Ray: ")
            s[0].print()

            for curve in s:
                curve.print(as_points)
                print("")
            print("============================")
