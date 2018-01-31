import numpy as np

class dig:

    def __init__(self, p_large = 0.5, L = 100, w = 5, mean_throw = 35, scale = 20, ratio = 1.):
        self.grid = np.random.rand(L,L) + 1 + p_large
        self.grid = self.grid.astype('int').tolist()
        self.redrop = []
        self.d = 0
        self.L = L
        #self.z = np.zeros(L,dtype='int')
        self.w = w
        self.p = p_large
        self.mean_throw = mean_throw
        self.scale = scale
        self.ratio = ratio
        self.dump = []

    def run(self,depth = None):
        if depth == None:
            depth = int(0.9 * self.L)
        while self.d < depth:
            self.remove()
            self.micro_run()
        x = []
        for i in self.grid:
            x.append(i[-1])
            self.surface = np.copy(x)

    def remove(self):
        for i in range(self.w):
            for j in range(self.w):
                x = self.grid[self.L - (j + 1)].pop()
                r = -1
                while r <= 0:
                    if x == 2:
                        r = int(np.random.normal(self.mean_throw*self.ratio,self.scale*self.ratio)) - (self.d/2)
                    else:
                        r = int(np.random.normal(self.mean_throw,self.scale)) - (self.d/2)
                if r > self.L:
                    self.dump.append(x)
                else:
                    self.grid[self.L - r].append(x)
                    #self.z[self.L - r] += 1
        #self.z[-(self.w + 1)] += self.w

    def micro_run(self):
        sm = 1
        while sm != 0:
            sm = 0
            for i in range(self.L - 1):
                z = len(self.grid[i]) - len(self.grid[i+1])
                if z <= 0:
                    pass
                else:
                    zc = ((2 * self.grid[i][-2]) - self.grid[i][-1]) + 1
                    if z > zc:
                        sm += 1
                        self.grid[i+1].append(self.grid[i].pop())
        self.d = self.L - len(self.grid[-1])

    def plot_pile(self):
        l = 0
        for i in self.grid:
            if len(i) > l:
                l = len(i)
        for i in self.grid:
            if len(i) != l:
                i += [0] * (l - len(i))
        return np.array(self.grid).T
