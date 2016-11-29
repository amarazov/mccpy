from qrng import *

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy as sc
import scipy.sparse.linalg as sla
from time import gmtime, strftime, clock
from mpi4py import MPI
import math
import numpy.linalg
import pickle

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def readMatrix(file_name):
    s = open(file_name,'r')
    X = []
    for line in s:
        X.append(eval(line))
    s.close()
    return np.array(X)

def parseMatrixMarketFile(filename):
    f = open(filename, 'r')
    header = f.readline()
    symmetric = header.find("symmetric")>0
    line = f.readline()
    while(line.find(r"%")>=0):
        line = f.readline()
    args = line.split()    
    row = []
    col = []
    nz = []
    for line in f:
        vals = line.split()
        row.append(int(vals[0])-1)
        col.append(int(vals[1])-1)
        nz.append(float(vals[2]))
        if(symmetric and vals[1] != vals[0]):
            row.append(int(vals[1])-1)
            col.append(int(vals[0])-1)
            nz.append(float(vals[2]))
    res = sp.csr_matrix((sc.array(nz),(sc.array(row),sc.array(col))),shape=(int(args[0]),int(args[1])))
    return res
        
def enlarge(A, nq):
    q =  sc.array([1.0/(nq-j+1) if i >= j else 0 for i in range(nq) for j in range(nq)])
    Q = np.reshape(q, (nq,nq), order='F')
    
    return sp.csr_matrix(sp.kron(A, Q))
    
def main():
    # IO only in the root
    A = []
    X = []
    if rank==0:
       A = pickle.load( open( "matA.p", "rb" ) )
       X = pickle.load( open( "matX.p", "rb" ) )
    
    A = comm.bcast(A, root = 0)
    X = comm.bcast(X, root = 0)
    
    
    seq_length=2**16
    max_scrambled_digits=30
    # define RNGs
    numseq = [Faur(2,seq_length,max_scrambled_digits,iflag) for iflag in (0,1,2,3)]
    #for iflag in (0,): numseq.append(Sobol(2,seq_length,max_scrambled_digits,iflag)) 
    numseq.insert(0, RNG(1))


    # iteration parameters
    mm = [5]
    NN = [1000, 10000, 100000]
    repeat = 5
    
    # matricies
    matrix_eigenvalue = (('A', A, 1.0),('X', X, 1.0))
    results=[]
    i=0
    experiment_start = 0.0
    if rank == 0:
        experiment_start = clock()
        
    for seq in numseq:
        i+=1
        if rank == i % size:
            print "Start tests with", seq.rng_type(), " on core ", rank
            for name, X,eigv in matrix_eigenvalue:
                print "matrix ", name
                for m in mm:
                    print "m=",m
                    for N in NN:
                        print "N=",N
                        for k in range(repeat):
                            seq.restart()
                            start = clock()
                            res = eigmcPowerMethod(X, m, N, seq)
                            end = clock()
                            results.append({'generator':seq.rng_type(), 'matrix':name, 'm':m, 'N':N, 
                                            'result':res, 'eigval':round( eigv,4), 'error':round(res - eigv, 4), 
                                            'time':round((end - start), 4)})
    total = comm.gather(results, root=0)
    
    if rank==0:
        print "Experiment finished in ", round((clock() - experiment_start), 4)
        total = [item for sublist in total for item in sublist]               
        df = pd.DataFrame(total)
        df.to_csv('../results/' + 'run' + strftime("%Y-%m-%d-%H-%M-%S", gmtime())+ '.csv')
        print df

if __name__ == "__main__":
    main()
    
