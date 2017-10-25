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
        self.__z_c = np.random.randint(1,2,L)
        self.s = []
        self.d = []
        self.point = ([],[])

    def micro_run(self):
        """Docstring"""
        tm = 0
        sm = 0
        dm = 0
        for i in self.point[tm%2]:
            if self.__z[i] - self.__z_c[i] > 0:
                self.__z[i] -= 2
                self.__z_c[i] = np.random.randint(1,2)
                sm += 1
                if i == 0:
                    self.point[(tm+1)%2].append(i+1)
                    self.__z[i+1] += 1
                elif i == self.__L - 1:
                    self.point[(tm+1)%2].append(i-1)
                    self.point[(tm+1)%2].append(i)
                    self.__z[i-1] += 1
                    self.__z[i] += 2
                    dm += 1
                else:
                    self.point[(tm+1)%2].append(i+1)
                    self.__z[i+1] += 1
                    self.point[(tm+1)%2].append(i-1)
                    self.__z[i-1] += 1
            self.point[tm%2] = []
            tm += 1
            self.s.append(sm)
            self.d.append(dm)

    def run(self,N):
        """Docstring"""
        t = 0
        index = np.arange(self.__L)
        z_t = ((self.__z - self.__z_c) > 0)
        self.point[0] = index[z_t]
        for i in range(N):
            t += 1
            self.micro_run()

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
