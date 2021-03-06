import numpy as np
import matplotlib.pyplot as plt

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
    deficit = deficit[:-1]
    correlate = correlate[1:]
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
    i = np.arange(1,len(slopes) + 1)
    return np.sum(slopes * i)

def state_gen(L):
    N = []
    x = np.ones(L)
    for i in x:
        N.append(Ng(x))
        i = np.argwhere(x == 1).flatten()
        x[i[np.random.randint(0,len(i))]] = 2
    N.append(Ng(x))
    return np.array(N)

def state_dist(N,L):
    dump = []
    for i in range(N):
        dump.append(state_gen(L))
    return np.block(dump)
