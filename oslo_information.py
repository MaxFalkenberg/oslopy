import numpy as np
import matplotlib.pyplot as plt
import sympy

def pile_excess(data, mode = 'outflux', L=None):
    """ Calculate deviation between net outflux/avalanches and expected
    outflux/avalanches at each time step.

        Parameters:

            data:       numpy.ndarray, shape(t), dtype('int')
                        Numpy array containing total number of
                        outflux/avalanche grains at each time step
                        corresponding to index in 1d array. Total time t.

            mode:       'outflux' or 'avalanches'
                        String defining mode. Either calculates for outflux
                        or avalanche data. If mode == 'avalanches', must
                        define size of pile from which excess is calculated.

            L:          integer
                        Integer size of Oslo Ricepile to calculate deviation
                        of avalanche behaviour.
        Returns:

            excess:     numpy.ndarray, shape(t), dtype('int')
                        Numpy array containing excess data at each timestep.

    """

    if mode == 'outflux':
        d = np.array(data) - 1 #Expected outflux is 1 grain each timestep.
    if mode == 'avalanches':
        if type(L) != int:
            raise ValueError('Must input pile size using variable L in integer\
                                type to calculate excess data for avalanches')
        else:
            d = np.array(data) - L #Expected s is L grains each timestep.
    excess = []
    j = 0
    for i in d:
        j += i
        excess.append(j)
    return np.array(excess).astype('float') - np.mean(excess)

def excess_compress(excess):
    """Convert excess data (from function 'pile_excess') into compressed form
        where an entry of 1 corresponds to more grains in pile than expected
        and 0 corresponds to fewer grains in pile than expected.

        Parameters:

            excess:     numpy.ndarray, shape(t), dtype('int')
                        Excess data produced from 'pile_excess' function. Total
                        length corresponds to t timesteps.

        Returns:

            excess_compressed:  numpy.ndarray, shape(t), dtype('bool')
                                Compressed excess data of length t.

    """

    mean = np.mean(excess)
    return (excess < mean)

def excess_chain(compress):
    """Converts information about excess grain data into integer array where
        each value corresponds to the number of consecutive timesteps where
        pile has been in excess state.

        Parameters:

            excess_compressed:  numpy.ndarray, shape(t), dtype('bool')
                                Compressed excess data of length t. Generated
                                from function 'excess_compress'.

        Returns:

            excess_chain:       numpy.ndarray, shape(t), dtype('int')
                                Array containing number of consecutive
                                timesteps where pile is in excess state.
    """
    chain = []
    j = 0
    for i in compress:
        j *= i
        j += i
        chain.append(j)
    return np.array(chain)

def chain_count(chain):
    """Calculates the probability that given the pile has been in the excess
        state for a given number of timesteps, the pile will remain in the excess
        state at the next time step.

        Parameters:

            excess_chain:       numpy.ndarray, shape(t), dtype('int')
                                Array containing number of consecutive
                                timesteps where pile is in excess state.
                                generated from function 'excess_chain'.

        Returns:

            chain_count:        numpy.ndarray, shape(np.max(excess_chain)),
                                dtype('float')
    """
    chain_max = np.max(chain[:-1])
    prob = np.zeros(chain_max+1, dtype = 'float')
    count = np.zeros(chain_max+1, dtype = 'float')

    for i in range(chain_max+1):
        i_arg = np.argwhere(chain[:-1] == i)
        t = float(len(np.argwhere(chain[i_arg + 1] > 0)))/len(i_arg)
        prob[i] = t
        count[i] = len(i_arg)

    return prob,count

def prob_chain(chain_count):
    """Total probability that the pile stays in an excess state for a given
        number of timesteps.

        Parameters:

            chain_count:        numpy.ndarray, shape(np.max(excess_chain)),
                                dtype('float')

        Returns:

            chain_prob:         numpy.ndarray, shape(np.max(excess_chain)),
                                dtype('float')

    """
    prob_chain = []
    j = 1
    for i in chain_count:
        j *= i
        prob_chain.append(j)
    return np.array(prob_chain)


def block_binary(compress,blocksize):
    """For a given blocksize, convert binary string of compressed data (from
        function 'excess_compressed') for consecutive chunks of length
        'blocksize' into a decimal number in range 0 to (2**blocksize) - 1.

        Parameters:

            compress:   numpy.ndarray, shape(t), dtype('bool')
                        Compressed excess data. Total length corresponds to t
                        timesteps.

            blocksize:  int
                        Length of desired blocksize in integers. blocksize must
                        be smaller than length of compressed data.

        Returns:

            excess_dec: numpy.ndarray, shape(t + 1 - blocksize), dtype('int')
                        Numpy array of decimal integers converted to binary
                        strings of length 'blocksize'.
    """

    dec = []
    powers = 2 ** np.arange(blocksize)
    for i in range(len(compress) + 1 - blocksize):
        dec.append(np.sum(compress[i:i+blocksize]*powers))
    return np.array(dec)

def excess_difference(excess,T, single = False):
    """Docstring
    """
    count_bin = []
    range_bin = []

    if single:
        dif1 = excess[T:-T] - excess[:-2*T]
        dif2 = excess[2*T:] - excess[:-2*T]
        count = np.zeros(np.ptp(dif1) + 1, dtype='float')
        d_range = np.arange(np.min(dif1),np.max(dif1) + 1)
        for j in range(len(d_range)):
            ind = np.argwhere(dif1 == d_range[j])
            if len(ind) != 0:
                count[j] = (float(len(np.argwhere(dif2[ind] > dif1[ind])))/
                                    len(ind))
            else:
                count[j] = None
        count_bin.append(count)
        range_bin.append(d_range)
    else:
        for i in range(1,T+1):
            dif1 = excess[i:-i] - excess[:-2*i]
            dif2 = excess[2*i:] - excess[:-2*i]
            count = np.zeros(np.ptp(dif1) + 1, dtype='float')
            d_range = np.arange(np.min(dif1),np.max(dif1) + 1)
            for j in range(len(d_range)):
                ind = np.argwhere(dif1 == d_range[j])
                if len(ind) != 0:
                    count[j] = (float(len(np.argwhere(dif2[ind] > dif1[ind])))/
                                        len(ind))
                else:
                    count[j] = None
            count_bin.append(count)
            range_bin.append(d_range)
    return range_bin,count_bin


def deficit(outflux,L,offset = 0):
    """Docstring
    """
    outflux = np.array(outflux)
    outflux -= 1
    outflux *= -1
    N_g = []
    j = offset
    for i in outflux:
        j += i
        N_g.append(j)
    N_g = np.array(N_g,dtype = 'float')
    deficit = 2. - 2.*N_g/(L * (L + 1.))
    return deficit

def cr_def(deficit,correlate):
    """Input pile.cor or pile.r as correlate.
    """
    correlate = np.array(correlate)
    deficit = deficit[:]
    correlate = correlate[:] #This deficit calculation only valid for correlated calculations!
    u = np.unique(deficit)

    dump = []
    l = []

    for i in u:
        ind = np.argwhere(deficit == i).flatten()
        dump.append(np.mean(correlate[ind]))
        l.append(len(ind))

    l = np.array(l)
    dump = np.array(dump)
    return u,dump,l

def Ng(slopes):
    """Docstring
    """
    i = np.arange(1,len(slopes) + 1)
    return np.sum(slopes * i)

def state_gen(L):
    """Docstring
    """
    N = []
    x = np.ones(L)
    for i in x:
        N.append(Ng(x))
        i = np.argwhere(x == 1).flatten()
        x[i[np.random.randint(0,len(i))]] = 2
    N.append(Ng(x))
    return np.array(N)

def state_dist(N,L):
    """Docstring
    """
    dump = []
    for i in range(N):
        dump.append(state_gen(L))
    return np.block(dump)

def slope_gen(L):
    """Generate spectrum of slope and critical slope configurations for an
        approximately uniform distribution of deficit values for a pile of
        size L.
    """
    x = np.ones((L+1,L))
    zc = np.ones((L+1,L))
    zc[0] = np.random.randint(1,3,L)
    N = [Ng(x[0])]
    for i in range(1,len(x)):
        x[i] = x[i-1]
        j = np.argwhere(x[i] == 1).flatten()
        x[i][j[np.random.randint(0,len(j))]] = 2
        zc[i] = x[i]
        k = np.argwhere(zc[i] == 1).flatten()
        zc[i][k] = np.random.randint(1,3,len(k))
        N.append(Ng(x[i]))
    return x,zc,np.array(N)


def crit_gen(z):
    """Docstring
    """
    zc = np.copy(z)
    zc[z==1] = np.random.randint(1,3,size = len(np.argwhere(z==1)))
    return zc

def lin_bin(x,y):
    b = np.linspace(0.02,1.,50)
    xb = []
    yb = []
    for i in b:
        t1 = x<i
        t2 = x>= i - 0.02
        t = t1 * t2
        xb.append(i-0.01)
        yb.append(np.mean(y[t]))
    return np.array(xb),np.array(yb)

def open_data(foldername):
    x = np.load(foldername + '/d_offset_norm.npy')
    y = np.load(foldername + '/r.npy')
    z = np.load(foldername + '/cor.npy')
    return x,y,z

def slope_block_init(L):
    #This method preferentially result in long blocks dominating
    #generated slope. Weight according to length.
    primes = list(sympy.primerange(0,L+1))
    blocks = [np.array([1])]
    for i in primes:
        b = np.ones(i,dtype = 'int8')
        b[-1] = 2
        b[np.random.randint(0,i-1)] = 0
        blocks.append(b)
    l = [1] + primes
    w = 1. / np.array(l)
    slope = []
    while L != 0:
        w /= np.sum(w[:sympy.primepi(L)+1])
        r = np.random.choice(sympy.primepi(L)+1,1,p = w[:sympy.primepi(L)+1])
        slope = [blocks[int(r)]] + slope
        L -= l[int(r)]
    slope = np.block(slope)
    slope_c = np.copy(slope)
    m = np.argwhere(slope_c != 2).flatten()
    slope_c[m] = np.random.randint(1,3,len(m))
    return slope, slope_c

def slopes(z_i,zc_i):
    z = [z_i]
    zc = [zc_i]
    N = [Ng(z_i)]
    l = len(np.argwhere(z_i != 2).flatten())
    while l != 0:
        z.append(np.copy(z[-1]))
        zc.append(np.copy(zc[-1]))
        i = np.argwhere(z[-1] != 2).flatten()
        zc[-1][i] = np.random.randint(1,3,len(i))
        r = np.random.randint(0,len(i))
        if z[-1][i[r]] == 1:
            z[-1][i[r]] = 2
            zc[-1][i[r]] = 2
        else:
            z[-1][i[r]] = 1
            zc[-1][i[r]] = np.random.randint(1,3)
        l = len(i) - 1
        N.append(Ng(z[-1]))
    z = np.vstack(z)
    zc = np.vstack(zc)
    N = np.array(N)
    return z,zc,N
