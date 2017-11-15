#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Python implementation of the Oslo Ricepile model.
"""


import numpy as np

class Oslo:
    """ Docstring """

    def __init__(self, L,mode = 'n'):
        if type(L) != int:
            raise ValueError("Grid size, L, must be integer type.")
        self.__L = L
        self.__t = 0
        self.__z = np.zeros(L,dtype='int')
        self.__z_c = np.random.randint(1,3,L)
        self.s = []
        self.d = []
        self.cor = []
        self.r = []
        self.point = [[],[]]
        if mode == 'r':
            self.__z += 2
            Oslo.run(self,1)
            self.__t = 0
            self.s = []
            self.d = []
            self.cor = []
            self.r = []
            self.point = [[],[]]
            self.d_offset = Oslo.height(self)
        else:
            self.d_offset = 0

    def height(self):
        j = 0
        h = 0
        for i in self.__z[::-1]:
            j += i
            h += j
        return h

    def custom_z(self,X):
        """Input custom pile configuration.

            Parameters:

                X:      numpy.ndarray, shape(L)
                        Numpy array with entries in range {0,1,2,3}.
        """
        X = X.astype('int')
        if np.shape(X) != (self.__L,):
            raise ValueError('Input array is not of shape (L,)')
        if np.all(np.in1d(X,[0,1,2,3])):
            self.__z = X
        else:
            raise ValueError('Custom array contains values other\
                                than [0,1,2,3]')



    def newslope(self):
        if np.random.random() > 0.5:
            return 1
        else:
            return 2

    def micro_run(self):
        """Docstring"""
        tm = 0
        sm = 0
        dm = 0
        cor = False
        r = 0
        # print('break')
        while len(self.point[tm%2]) != 0:
            # print(self.point[tm%2],self.point,tm)
            for i in self.point[tm%2]:
                if i == self.__L - 1:
                    cor = True
                if i > r:
                    r = i
                if self.__z[i] > self.__z_c[i]:
                    self.__z[i] -= 2
                    self.__z_c[i] = self.newslope()
                    sm += 1
                    if i == 0:
                        self.point[(tm+1)%2].append(i+1)
                        self.__z[i+1] += 1
                    elif i == self.__L - 1:
                        self.point[(tm+1)%2].append(i-1)
                        self.point[(tm+1)%2].append(i)
                        self.__z[i-1] += 1
                        self.__z[i] += 1
                        dm += 1
                    else:
                        self.point[(tm+1)%2].append(i+1)
                        self.__z[i+1] += 1
                        self.point[(tm+1)%2].append(i-1)
                        self.__z[i-1] += 1
            self.point[tm%2] = []
            tm += 1
        self.cor.append(cor)
        self.r.append(r)
        self.s.append(sm)
        self.d.append(dm)

    def run(self,N):
        """Docstring"""
        index = np.arange(self.__L)
        self.__z[0] += 1
        z_t = ((self.__z - self.__z_c) > 0)
        self.__z[0] -= 1
        self.point[0] = list(index[z_t])
        checks = 0
        for j in range(N):
            if 0 not in self.point[0]:
                self.point[0].append(0)
            self.__z[0] += 1
            self.__t += 1
            self.micro_run()
            # print(self.__z)

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

    def micro_run(self):
        """Internal function.
            Executes all topplings of the ricepile at the macroscopic time
            the function is called. Records total number of toppings (avalance
            size) and the total number of grains leaving the system at the open
            boundary (drop size)."""
        tm = 0 #microscopic counter (not time!)
        sm = 0 #avalance index
        dm = 0 #drop index
        boundary_topple = True #has the boundary site toppled at last microtime

        while boundary_topple:
            boundary_topple = False
            continue_topple = True #any sites left to topple
            while continue_topple:
                print(tm)
                #Sites with possible toppings
                print(self.index_choice[tm])
                print(self.topple_check[tm])
                dep_tm = self.index_choice[tm][self.topple_check[tm]]
                z_dif = (self.__z[dep_tm] -
                            self.__z_c[dep_tm])
                #Sites that do topple
                z_top = (z_dif + np.absolute(z_dif)).astype('bool')
                top_sum = np.sum(z_top)
                sm += top_sum
                continue_topple = bool(top_sum)
                if continue_topple: #Reset critical slopes at toppled sites
                    self.__z_c[dep_tm[z_top]] = np.random.randint(1,2,top_sum)
                    self.__z[dep_tm[z_top]] -= 2
                    if tm == self.__L - 1:
                        self.topple_check[tm-1] *= 0
                        if (tm+1)%2:
                            self.topple_check[tm-1] += z_top
                            self.topple_check[tm-1][:-1] += z_top[1:]
                            self.__z[dep_tm[z_top[1:]] - 1] += 1
                        else:
                            self.topple_check[tm-1] += z_top[1:]
                            self.topple_check[tm-1] += z_top[:-1]
                            self.__z[dep_tm[z_top] - 1] += 1
                        if z_top[-1]:
                            self.__z[-1] += 1
                        tm -= 1
                        boundary_topple = z_top[-1]
                        dm += boundary_topple
                    elif tm == self.__L - 2 and boundary_topple:
                        if (tm+1)%2:
                            self.topple_check[tm+1] += z_top
                            self.topple_check[tm+1][:-1] += z_top[1:]
                            self.topple_check[tm+1][-1] = True
                            self.__z[dep_tm[z_top[1:]] - 1] += 1
                            self.__z[dep_tm[z_top] + 1] += 1
                        else:
                            self.topple_check[tm+1][:-1] += z_top
                            self.topple_check[tm+1][1:] += z_top
                            self.topple_check[tm+1][-1] = True
                            self.__z[dep_tm[z_top] - 1] += 1
                            self.__z[dep_tm[z_top] + 1] += 1
                        tm += 1
                        continue_topple = True
                        #Ensures boundary is checked again for double toppling
                    else:
                        if (tm+1)%2:
                            self.topple_check[tm+1] += z_top
                            self.topple_check[tm+1][:-1] += z_top[1:]
                            self.__z[(dep_tm[z_top] - 1)[1:]] += 1
                            self.__z[dep_tm[z_top] + 1] += 1
                        else:
                            self.topple_check[tm+1][:-1] += z_top
                            self.topple_check[tm+1][1:] += z_top
                            self.__z[dep_tm[z_top] - 1] += 1
                            self.__z[dep_tm[z_top] + 1] += 1
                        tm += 1

        self.s.append(sm)
        self.d.append(dm)

a = Oslo(10)
b = np.ones(10,dtype='int')
b[0] = 2
a.custom_z(b)
a.micro_run()
