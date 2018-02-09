import numpy as np
import matplotlib.pyplot as plt
import copy
from random import shuffle

def density(L = 150.):
    x = np.linspace(0.1,L,1000).astype('float')
    d2 = (x**2)/(2*(L - x))
    d3 = (x ** 3)/(12*(((L/2)**2) - ((x/2)**2)))
    plt.plot(x,d2,label = '2D Density')
    plt.plot(x,d3,label = '3D Density')
    plt.plot([50,50],[0.1,100],label = 'Width = 50',ls = '--')
    plt.xlabel('Pit Width [Grain Units]', fontsize = 18)
    plt.ylabel('Grain Density', fontsize = 18)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc = 'best',fontsize=18)
    plt.savefig('graindensity.png')
    plt.show()

def lining_count(p = 0.5):
    N = 25
    x = [2,4,6,8,10,15,20,30,50,100,200,300,500,1000]
    dx = [2,2,2,2,2,5,5,10,20,50,100,100,200,500]
    yraw = [[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
    for i in range(N):
        print i
        a = Oslo(351,p_large=p,mode='dig',w=5)
        for j in range(len(dx)):
            a.run(dx[j])
            a.count()
            yraw[j].append(np.mean(a.lining))
    y = np.array([np.mean(i) for i in yraw])
    yerr = np.array([np.std(i,ddof=1) for i in yraw])
    return x,y,yerr



class Oslo:

    def __init__(self,L,p_large=0.5,mode='build',centre=None,w=1,grain_throw = 'bucket'):
        self.L = int(L)
        self.z = np.zeros((2,L),dtype='int8')
        self.p = p_large
        self.rdrop = []
        self.ldrop = []
        self.removed = []
        self.w = w
        self.width_dump = []
        self.depth_dump = []
        self.topview = []
        self.lr = []
        self.t = 0
        self.removal_window = np.zeros((self.w,self.w),dtype = 'float')
        self.grain_throw = grain_throw
        if mode == 'dig':
            # self.grains = [[0]*L for i in range(L)]
            self.grains = np.random.rand(L,L) + 1 + p_large
            self.grains = self.grains.astype('int').tolist()
        elif mode == 'build':
            self.grains = [[] for i in range(L)]
        else:
            raise ValueError('mode variable must be set to \'build\' or \'dig\'.')
        self.mode = mode
        if centre==None:
            self.centre = int(L - 1)/2
        else:
            self.centre = centre
        self.check = np.zeros(self.L,dtype='bool')

    def plot_pile(self):
        grid = copy.deepcopy(self.grains)
        l = 0
        for i in grid:
            if len(i) > l:
                l = len(i)
        for i in grid:
            if len(i) != l:
                i += [0] * (l - len(i))
        plt.clf()
        plt.figure()
        plt.imshow(np.array(grid).T, origin='lower', interpolation = 'none')
        plt.colorbar()
        plt.show()
        return np.array(grid).T

    def count(self):
        l = True
        r = True
        len_t = np.array([len(i) for i in self.grains])
        rlim = self.centre
        llim = self.centre
        lining = [self.grains[self.centre][-1]]


        for i in range(self.centre+1,self.L):
            if r:
                # if len_t[i] > len_t[i-1]:
                #      lining += self.grains[i][len_t[i-1]:]
                # else:
                #     lining.append(self.grains[i][-1])
                lining.append(self.grains[i][-1])
                if len_t[i] >= self.L:
                    rlim = i
                    r = False
                # else:
                #     r = False

        for i in range(self.centre-1,-1,-1):
            if l:
                # if len_t[i] > len_t[i+1]:
                #     lining += self.grains[i][len_t[i+1]:]
                # else:
                #     lining.append(self.grains[i][-1])
                lining.append(self.grains[i][-1])
                if len_t[i] >= self.L:
                    llim = i
                    l = False
                # else:
                #     l = False

        w = rlim - llim
        self.width = w
        self.width_dump.append(w)
        self.lr.append([llim,rlim])
        self.depth_dump.append(len_t[self.centre])
        self.topview.append([j[-1] for j in self.grains])
        self.lining = lining
        # print w, llim, rlim, lining

    def run(self,N=1):
        if self.mode=='build':
            self.check_sites = [self.centre]
            self.check[self.centre] = True
            for i in range(N):
                self.t += 1
                if np.random.rand() > self.p:
                    self.grains[self.centre].append(1)
                else:
                    self.grains[self.centre].append(2)
                self.z[:,self.centre] += 1
                if self.centre != self.L - 1:
                    self.z[:,self.centre+1][0] -= 1
                if self.centre != 0:
                    self.z[:,self.centre-1][1] -=1
                self.micro_run()
        else:
            # if N>=self.L:
            #     raise ValueError('N cannot exceed L for mode=\'dig\'')
            self.check[:] = True
            self.check_sites = range(self.L)
            for i in range(N):
                self.t+=1
                if len(self.grains[self.centre]) <= self.w:
                    print('Digging ended because lower limit has been reached.')
                    print('Digging ended after ' + str(self.t) + ' iterations.')
                    break
                self.count()
                self.remove()
                self.micro_run()

    def remove(self):

        r = np.arange(self.centre,self.centre + self.w)
        r -= int(self.w - 1)/2
        if np.min(r) < 0 or np.max(r) >= self.L:
            raise ValueError('Removal window size and location overlaps with system boundaries.')
        else:
            for i in range(self.w):
                for j in r:
                    g = self.grains[j].pop()
                    self.removal_window[i][j-r[0]] += (g - 1)
                    z = len(self.grains[j])
                    if self.grain_throw != 'bucket':
                        p = self.redistribute(g,i,j,z)
            self.z[:,r[0]][0] -= self.w
            self.z[:,r[-1]][1] -= self.w
            if r[0] != 0:
                self.z[:,r[0]-1][1] += self.w
            if r[-1] != self.L - 1:
                self.z[:,r[-1]+1][0] += self.w

    def micro_run(self):
        cont = 1
        while cont:
            sm = 0
            shuffle(self.check_sites)
            for i in self.check_sites:
                if len(self.grains[i]) < 2:
                    pass
                else:
                    zc = (2*self.grains[i][-2]) - self.grains[i][-1] + np.random.randint(0,2)
                    # if i == 0:
                    #     pass
                    # elif i == self.L -1:
                    #     pass
                    # else:
                    lr = self.z[:,i]>zc
                    s = np.sum(lr)
                    sm += s
                    if s == 0:
                        pass
                    else:
                        if s == 2:
                            if np.random.rand() > 0.5:
                                lr[0] = True
                                lr[1] = False
                            else:
                                lr[1] = True
                                lr[0] = False
                        if lr[0]:
                            if i == 0:
                                self.ldrop.append(self.grains[i].pop())
                                self.z[:,i] -= 1
                                self.z[:,i+1][0] += 1
                            elif i == self.L - 1:
                                self.grains[i-1].append(self.grains[i].pop())
                                self.z[:,i][0] -= 2
                                self.z[:,i][1] -= 1
                                self.z[:,i-1][0] += 1
                                self.z[:,i-1][1] += 2
                                self.z[:,i-2][1] -= 1
                                if self.check[i-1] == False:
                                    self.check_sites.append(i-1)
                                    self.check[i-1] = True
                            else:
                                self.grains[i-1].append(self.grains[i].pop())
                                self.z[:,i][0] -= 2
                                self.z[:,i][1] -= 1
                                self.z[:,i-1][0] += 1
                                self.z[:,i-1][1] += 2
                                self.z[:,i+1][0] += 1
                                if i != 1:
                                    self.z[:,i-2][1] -=1
                                if self.check[i-1] == False:
                                    self.check_sites.append(i-1)
                                    self.check[i-1] = True
                        elif lr[1]:
                            if i == 0:
                                self.grains[i+1].append(self.grains[i].pop())
                                self.z[:,i][0] -= 1
                                self.z[:,i][1] -= 2
                                self.z[:,i+1][0] += 2
                                self.z[:,i+1][1] += 1
                                self.z[:,i+2][0] -= 1
                                if self.check[i+1] == False:
                                    self.check_sites.append(i+1)
                                    self.check[i+1] = True
                            elif i == self.L - 1:
                                self.rdrop.append(self.grains[i].pop())
                                self.z[:,i] -= 1
                                self.z[:,i-1][1] += 1
                            else:
                                self.grains[i+1].append(self.grains[i].pop())
                                self.z[:,i][0] -= 1
                                self.z[:,i][1] -= 2
                                self.z[:,i+1][0] += 2
                                self.z[:,i+1][1] += 1
                                self.z[:,i-1][1] += 1
                                if i != self.L - 2:
                                    self.z[:,i+2][0] -=1
                                if self.check[i+1] == False:
                                    self.check_sites.append(i+1)
                                    self.check[i+1] = True
            if sm == 0:
                cont = 0

    def redistribute(self,gsize,layer,x0,z0):
        v0 = [70.,70.,70.,70.,70.]
        dt = (np.random.rand()*20)-10
        dv = (np.random.rand()*60)-30
        x,z = self.traj(gsize, v0[layer] + dv, 50. + dt)
        if np.random.rand() > 0.5:
            x = -x
        x += x0
        z += z0
        x = x.astype('int')
        z = z.astype('int')
        xr = np.arange(self.L)
        yr = [len(self.grains[i]) for i in xr]
        for k in xr:
            if k == 0:
                pass
            else:
                if x[k] > self.L - 2 or x[k] < 1:
                    self.removed.append(gsize)
                    return 0
                elif z[k] < yr[x[k]]:
                    # print x[k], v0[layer]
                    # self.removed.append(gsize)
                    self.grains[x[k]].append(gsize)
                    self.z[:,x[k]] += 1
                    self.z[:,x[k]-1][1] -= 1
                    self.z[:,x[k]+1][0] -= 1
                    return 0

    def traj(self,size = 1,v0=100.,theta=50.):
        n = 101
        g = 981
        theta *= (np.pi/180.)
        vx = v0 * np.cos(theta)
        vz = v0 * np.sin(theta)
        t = np.linspace(0,0.2,101)
        if size == 1:
            vt = 150.
        else:
            vt = 450.
        if self.grain_throw == 'drag':
            x = ((v0 * vt * np.cos(theta))/g)*(1 - np.exp(-((g * t)/vt)))
            z = (vt/g)*(v0 * np.sin(theta) + vt)*(1 - np.exp(-((g * t)/vt))) - (vt * t)
        else:
            x = v0 * np.cos(theta) * t
            z = (v0 * np.sin(theta) * t) - ((g/2)*(t**2))
        x *= 10
        z *= 35
        return x,z
