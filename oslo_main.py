#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Python implementation of the Oslo Ricepile model.
"""


import numpy as np
import pickle
import os
import binascii
import oslo_information as oi

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

    def save(self, foldername = None):

        files = (self.__L,self.__t,self.point)
        if foldername == None:
            folder = str('filedump_L' + str(self.__L) + '_t' + str(self.__t) +
                                        '_' + binascii.b2a_hex(os.urandom(6)))
        else:
            folder = foldername
        os.makedirs(folder)
        with open(folder + '/meta.pickle', 'wb') as f:
            pickle.dump(files, f)
        np.save(folder + '/z',self.__z)
        np.save(folder + '/z_c',self.__z_c)
        np.save(folder + '/s',self.s)
        np.save(folder + '/d',self.d)

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

    def def_run(self,N):
        self.d_offset = []
        for i in range(N):
            if i % 5 == 0:
                print i
            z_stack,zc_stack,offset = oi.slope_gen(self.__L)
            self.d_offset.append(offset)
            for j in range(self.__L+1):
                self.__z = z_stack[j]
                self.__z_c = zc_stack[j]
                self.point = [[0],[]]
                self.__z[0] += 1
                self.__t += 1
                self.micro_run()
        self.d_offset = np.block(self.d_offset)
        self.d_offset_norm = 2. - 2.*self.d_offset/(self.__L * (self.__L + 1.))

    def def_save(self,foldername=None):
        files = {'L':self.__L,'t':self.__t,'N':float(self.__t)/(self.__L + 1)}
        if foldername == None:
            folder = str('defdata_dump_L' + str(self.__L) + '_t' + str(self.__t) +
                                        '_' + binascii.b2a_hex(os.urandom(6)))
        else:
            folder = foldername
        os.makedirs(folder)
        with open(folder + '/meta.pickle', 'wb') as f:
            pickle.dump(files, f)
        np.save(folder + '/d_offset',self.d_offset)
        np.save(folder + '/d_offset_norm',self.d_offset_norm)
        np.save(folder + '/s',self.s)
        np.save(folder + '/d',self.d)
        np.save(folder + '/r',self.r)
        np.save(folder + '/cor',self.cor)



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
