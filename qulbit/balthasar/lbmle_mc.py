#!/usr/bin/python                                                                  
# -*- coding: utf-8 -*-                                                            
#                                                                                  
# lbmle_mc.py: A helper class implementing Monte Carlo
# simulations to generate frequency data for least-biased 
# MLE tomography, done according to Rehacek et. al, 
# PRA 92, 052303 (2015)
#                                                                                  
# © 2017 Olivia Di Matteo (odimatte@uwaterloo.ca)                                  
#                                                                                  
# This file is part of the project Balthasar.                                      
# Licensed under BSD-3-Clause                                                      
# 

from random import random

import numpy as np

class LBMLE_MC():
    """ A wrapper class to hold Monte Carlo routines for 
        testing of least-biased MLE tomography. 

        Args:
            bases (list of lists): The complete set of MUBs
                                   for a given dimension.
                         
        Attributes:
            bases (list): Same as above.
    """

    def __init__(self, basis_vectors):
        """ Initialize a least-bias MLE container.
        """
        # Take the eigenvectors that are passed, and turn
        # them into projectors.
        self.projectors = self.generate_projectors(basis_vectors)
    

    def generate_projectors(self, vectors):
        """ Converts a set of vectors into projector form.
            Keeps order/format the same.
        """
        projectors = []
    
        for basis in vectors:
            projectors_this_basis = []
            for vector in basis:
                projectors_this_basis.append(np.outer(vector, np.conj(vector)))
            projectors.append(projectors_this_basis)
        return projectors

        
    def simulate(self, bases, state, N = 10000):
        """ Run a Monte Carlo simulation to generate measurement statistic 
            data for the provided basis indices.

            Args:
                bases (list): The indices of which bases we want to generate
                              frequencies for.
                state (matrix): An arbitrary density matrix to generate the 
                                statistics for.
                N (int): The number of 'particles'/trials in the simulation. 
                         Default value is 10000.
    
            Return:
                A list of frequencies which are the results of simulating 
                a set of measurements done in each basis listed in bases.
        """
        
        freqs = [] # Output data. A list of lists, one per provided basis.

        # Loop through bases one at a time.
        for basis_idx in bases:
            #print("Simulating measurement data in basis " + str(basis_idx))

            # Check that the user specified valid bases. Note that -1 is a valid option
            # as it corresponds to the "last" basis, i.e. the vertical slope ray.
            if basis_idx not in range(len(self.projectors)) and basis_idx is not -1:
                print("Error, invalid basis index '" + str(basis_idx) + "' provided.") 
                return
    
            # Gather the projectors
            projs = self.projectors[basis_idx]

            d = len(projs) # Number of possible measurement outcomes
            
            # Compute the probabilities of each measurement using the Born rule,
            # p_i = Tr(state * M_i)
            vn_probs = [float(np.trace(np.dot(state, projs[x])).real) for x in range(d)]
            
            # Create the MC simulation "range" array. This is the continuous
            # line of the probabilities we use to determine where the particle
            # lands in the simulation whenever a random value is drawn.
            mc_range = [0] * d 
            mc_range[0] = vn_probs[0]

            for i in range(1, d):
                mc_range[i] = mc_range[i-1] + vn_probs[i]

            # Vector to store the outcome of the experiments
            # Each slot represents a measurement outcome; every time    
            # we get that outcome (as per mc_range), we'll add a count to the slot.
            counts = [0] * d 

            # Go through N particles, generate a measurement outcome for each, and 
            # update the counts.
            for n in range(N):
                r = random()
        
                for outcome in range(d):
                    if outcome == 0:
                        if r >= 0 and r < mc_range[outcome]:
                            counts[outcome] += 1
                    else:
                        if r >= mc_range[outcome - 1] and r < mc_range[outcome]:
                            counts[outcome] += 1
    
            # Compute the frequencies via the MC simulation data and append
            these_freqs = [counts[x] / N for x in range(d)]
            freqs.append(these_freqs)

        return freqs












 
        

