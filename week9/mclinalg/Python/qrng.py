import sfaur
import ssobol
# import flinalg

import numpy as np
import scipy.sparse as sp

class RNG:
    def __init__(self, dimen, seed=0):
        self.s = dimen
        np.random.seed(seed)
        
    def __call__(self):
        return np.random.uniform(low=0.0, high=1.0, size=self.s)
    
    def rng_type(self):
        res = "RNG" 
        return res
    
    def restart(self, seed=0):
        np.random.seed()
        

class Sobol:
    """       User Define: 
        dimen : dimension 
        seq_length : sequence length
        max_scrambled_digits : Maximum Digits of Scrambling Of Owen type Scrambling
        IFLAG: User Choice of Sequences
        IFLAG = 0 : No Scrambling
        IFLAG = 1 : Owen type Scrambling
        IFLAG = 2 : Faure-Tezuka type Scrambling
        IFLAG = 3 : Owen + Faure-Tezuka type Scrambling
    """
    def __init__(self, dimen, seq_length, max_scrambled_digits, iflag):
        self.s = dimen
        self.atmost = seq_length
        self.max_scrambled_digits = max_scrambled_digits
        self.iflag = iflag
        res = ssobol.inssobl(dimen, seq_length, max_scrambled_digits, iflag)
        if not all(res[0]):
            raise Exception("Check parameters of Sobol")
    
    def __call__(self):
        return ssobol.gossobl()[1:self.s]
       
    def restart(self):
        res = ssobol.inssobl(self.s, self.atmost, self.max_scrambled_digits, self.iflag)
        if not all(res[0]):
            raise Exception("Check parameters of Sobol")

    
    def rng_type(self):
        res = "Sobol" 
        if(self.iflag == 0): res += "; No Scrambling"
        elif(self.iflag == 1): res += "; Owen type Scrambling"
        elif(self.iflag == 2): res += "; Faure-Tezuka type Scrambling"
        elif(self.iflag == 3): res += "; Owen + Faure-Tezuka type Scrambling"
        return res
    
    
class Faur:
    """       User Define: 
        dimen : dimension 
        seq_length : sequence length
        max_scrambled_digits : Maximum Digits of Scrambling Of Owen type Scrambling
        IFLAG: User Choice of Sequences
        IFLAG = 0 : No Scrambling
        IFLAG = 1 : Owen type Scrambling
        IFLAG = 2 : Faure-Tezuka type Scrambling
        IFLAG = 3 : Owen + Faure-Tezuka type Scrambling
    """
    def __init__(self, dimen, seq_length, max_scrambled_digits, iflag):
        self.s = dimen
        self.atmost = seq_length
        self.max_scrambled_digits = max_scrambled_digits
        self.iflag = iflag
        res = sfaur.insfaur(dimen, seq_length, max_scrambled_digits, iflag)
        if not all(res[0]):
            raise Exception("Check parameters of Faur")
    
    def __call__(self):
        # return sfaur.gosfaur(self.max_scrambled_digits)[0:self.s]
        return sfaur.gosfaur(self.max_scrambled_digits)[1:(self.s)]
    
    def restart(self):
        res = sfaur.insfaur(self.s, self.atmost, self.max_scrambled_digits, self.iflag)
        if not all(res[0]):
            raise Exception("Check parameters of Faur")
    
    def rng_type(self):
        res = "Faur" 
        if(self.iflag == 0): res += "; No Scrambling"
        elif(self.iflag == 1): res += "; Owen type Scrambling"
        elif(self.iflag == 2): res += "; Faure-Tezuka type Scrambling"
        elif(self.iflag == 3): res += "; Owen + Faure-Tezuka type Scrambling"
        return res

def normalize(A):
    (data, indices, indptr) = (A.data.copy(), A.indices, A.indptr)
    
    for i in range(len(indptr) - 1):
        row = np.abs(data[indptr[i]:indptr[i + 1]])
        rsum = np.sum(row)
        if(rsum > 0.0):
            data[indptr[i]:indptr[i + 1]] = row / rsum
        else:
            raise Exception("zero line")    
    return sp.csr_matrix((data, indices, indptr), shape=A.shape)
    
def getWeighted(vals, weights, rand):
    if (np.abs(np.sum(weights) - 1) > 1e-6):
        raise Exception("weigths should sum up to 1!")
    r = rand()
    s = 0.0
    i = 0
    while(s <= r):
        while(weights[i] <= 0):
            i += 1
        s += weights[i]
        i += 1
    return vals[i - 1]
#     return flinalg.getweightedint(vals, weights,r)
        
def eigmcPowerMethod(A, m, N, numseq):
    (n, nn) = A.shape
    if(n != nn):
        raise Exception("non square matrix")
    h = np.repeat(1.0, n)
    p0 = h / sum(h)
    P = normalize(A)
    sumW = np.float96(0.0)
    sumW_1 = np.float96(0.0)
    for i in range(N):
        k = getWeighted(range(n), p0, numseq)
        W = h[k] / p0[k]
        W_1 = W
        for j in range(m):
            k_1 = k
            vals = P.indices[P.indptr[k_1]:P.indptr[k_1 + 1]]
            weights = P.data[P.indptr[k_1]:P.indptr[k_1 + 1]]
            k = getWeighted(vals, weights, numseq)
            W_1 = W
            W *= A[k_1, k] / P[k_1, k] 
        sumW += W
        sumW_1 += W_1
    return sumW / sumW_1


