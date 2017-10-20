#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Python implementation of the Oslo Ricepile model.
"""


import numpy as np

class Oslo:
    """ Docstring """

    def __init__(self, L):
        if type(L) != int:
            raise ValueError("Grid size, L, must be integer type.")
        self.__L = L
        self.__t = 0
        self.__z = np.zeros(L,dtype='int')
        self.__z_c = np.random.randint(1,3,L)
        self.s = []
        self.d = []
        self.index_choice, self.topple_check, self.topple_dependencies = (
                self.index_dep_gen(L))

    def index_dep_gen(self, L):
        """Internal method for generating list of indices for possible toppling
            locations and toppling dependencies of each site."""
        index_choice = []
        topple_check = []
        dependencies = []
        for i in range(L):
            index_temp = np.arange(i%2,i+1,2)
            topple_check_temp = np.zeros_like(index_temp,dtype='bool')
            index_choice.append(index_temp)
            topple_check.append(topple_check_temp)
            if i == 0:
                dependencies.append([1])
            elif i == L-1:
                dependencies.append([L-2])
            else:
                dependencies.append([i-1,i+1])

        return index_choice, topple_check, dependencies

    def info(self,single = False):
        """Returns key information about current state of the ricepile.

            Returns tuple (L,t,z,z_c) if single == False.
            If single == i for any i in {0,1,2,3}, returns (L,t,z,z_c)[i].

            L:      int
                    System size of ricepile. Equal to number of grid spaces.

            t:      int
                    Macroscopic time of system. Increases by 1 each time a grain
                    is added to the pile.

            z:      numpy.ndarray, shape(L)
                    Current slope at each grid location. Grains propagate from
                    right to left. Grains are dissipated at right boundary.

            z_c:    numpy.ndarray, shape(L)
                    Current critical slope at each grid location. Possible
                    values at each site {1,2}.
        """
        if single not in [False,1,2,3,4]:
            raise ValueError("single must take value in [False,1,2,3,4]")
        data = (self.__L, self.__t, self.__z, self.__z_c)
        if single == False:
            return data
        else:
            return data[single]
