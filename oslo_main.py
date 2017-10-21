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
        self.__L = L #grid size
        self.__t = 0 #macroscopic time
        self.__z = np.zeros(L,dtype='int') #grid of ricepile slopes
        self.__z_c = np.random.randint(1,2,L) #critical slopes
        self.s = [] #avalance list
        self.d = [] #drop list
        self.index_choice, self.topple_check, self.topple_dependencies = (
                self.index_dep_gen(L))

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
                dependencies.append([0])
            elif i == L-1:
                dependencies.append([(L-2)//2])
            else:
                dependencies.append([(i-1)//2,(i+1)//2])
        topple_check[0] = np.array([1],dtype='bool')

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
