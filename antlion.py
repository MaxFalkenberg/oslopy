#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Python implementation of the Oslo Ricepile model.
"""


import numpy as np
import pickle
import os
import binascii
import copy

class Oslo:
    """ Docstring """

    def __init__(self, L,p = 0.5, load = False):
        self.p = p
        if type(L) != int:
            raise ValueError("Grid size, L, must be integer type.")
        if load != False:
            with open (load + '/meta.pickle', 'rb') as fp:
                self.__L,self.__t,self.point = pickle.load(fp)
            self.__z = np.load(load + '/z.npy')
            self.__z_c = np.load(load + '/z_c.npy')
            self.s = list(np.load(load + '/s.npy'))
            self.d = list(np.load(load + '/d.npy'))
        else:
            self.__L = L
            self.__t = 0
            self.__z = np.zeros(L,dtype='int')
            self.__z_c = [[3,2] for i in range(L)]
            self.s = []
            self.d = []
            self.point = [[],[]]
            self.datadump = []


    def newslope(self):
        if np.random.random() > self.p:
            return 1
        else:
            return 2

    # @profile
    def micro_run(self,mode= 'build'):
        """Docstring"""
        tm = 0
        sm = 0
        dm = 0
        topple = True
        # print('break')
        while topple:
            topple = False
            if self.__z[0] > (2 * self.__z_c[0][-2] - self.__z_c[0][-1] + 1):
                self.__z_c[1].append(self.__z_c[0][-1])
                del self.__z_c[0][-1]
                self.__z[0] -= 2
                self.__z[1] += 1
                sm+=1
                topple = True
            for i in range(1,self.__L -1):
                # print(i)
                if self.__z[i] > (2 * self.__z_c[i][-2] - self.__z_c[i][-1] + 1):
                    self.__z_c[i+1].append(self.__z_c[i][-1])
                    del self.__z_c[i][-1]
                    self.__z[i] -= 2
                    self.__z[i-1] += 1
                    self.__z[i+1] += 1
                    sm+=1
                    topple = True
            if mode == 'build':
                if self.__z[self.__L -1] > (2 * self.__z_c[self.__L -1][-2] - self.__z_c[self.__L -1][-1] + 1):
                    #self.__z_c[self.__L].append(self.__z_c[self.__L -1][-1])
                    del self.__z_c[self.__L -1][-1]
                    self.__z[self.__L -1] -= 1
                    self.__z[self.__L -2] += 1
                    sm+=1
                    dm+=1
                    topple = True

            tm += 1
        self.s.append(sm)
        self.d.append(dm)

    # @profile
    def run(self,N):
        """Docstring"""
        for i in range(N):
            self.__t += 1
            self.__z[0] += 1
            self.__z_c[0].append(self.newslope())
            self.micro_run()
            self.datadump.append(copy.deepcopy(self.__z_c))

    def dig(self):
        self.__z_c = []
        self.remove = []
        for i in range(self.__L):
            self.__z_c.append([self.newslope() for j in range(self.__L)])
        # for i in range(self.__L/2):
        #     self.__z_c.append([3,3])
        # self.__z[(self.__L/2)-1] = 3 * self.__L
        self.micro_run(mode = 'dig')
        while len(self.__z_c[-1]) > 50:
        # for j in range(10):
            # print self.__z_c[-5][-5:] + self.__z_c[-4][-5:] + self.__z_c[-3][-5:] + self.__z_c[-2][-5:] + self.__z_c[-1][-5:]
            # break
            self.remove.append(np.array(self.__z_c[-5][-5:] + self.__z_c[-4][-5:] + self.__z_c[-3][-5:] + self.__z_c[-2][-5:] + self.__z_c[-1][-5:]))
            self.remove[-1] = self.remove[-1].reshape((5,5)).T[::-1]
            del self.__z_c[-1][-5:]
            del self.__z_c[-2][-5:]
            del self.__z_c[-3][-5:]
            del self.__z_c[-4][-5:]
            del self.__z_c[-5][-5:]
            self.__z[-6] += 5
            self.micro_run(mode = 'dig')
            self.datadump.append(copy.deepcopy(self.__z_c))
        #return self.__z_c


    def plot_pile(self,z_c, mode = 'dig'):
        if mode == 'dig':
            grid =  np.zeros((self.__L,self.__L),dtype ='int')
        else:
            grid = np.zeros((3*self.__L +5,self.__L),dtype ='int')
        for i in range(self.__L):
            grid[:,i][:len(z_c[i])] = z_c[i]
        return grid[::-1]

    def animate(self, mode = 'dig'):
        if mode == 'dig':
            grid = np.zeros((len(self.datadump),self.__L,self.__L),dtype ='int')
        else:
            grid = np.zeros((len(self.datadump),3*self.__L +5,self.__L),dtype ='int')
        for i in range(len(self.datadump)):
            grid[i] = self.plot_pile(self.datadump[i])
        return grid




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

# a = Oslo(32)
# a.run(100000)
