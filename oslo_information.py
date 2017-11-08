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
    return np.array(excess)

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
