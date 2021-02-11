import numpy as np
from macarons import macarons
from util import computeSKAT

#   Load the data 
npzfile = np.load(r'../../data/data.npz')
#   X : genotype  - n x m  matrix - n : Samples, m : Features (SNPs)
X = npzfile['X']
#   Y : phenotype - n x 1 column vector
Y = npzfile['Y']
#   snp: positions and chromosomes of the SNPs - n x 2 matrix
#   snp[:, 1] = chromosome indices
#   snp[:, 2] = Position on the chromosome
snp = npzfile['snp']

nPCs = 1    #   Number of Principal Components used
            #   to correct population stratification
k = 200     #   Number of Features to be selected
D = 2e4     #   Intra-chromosomal distance in base pairs
            #   to limit the search space of Macarons

#   C : scores - m x 1 column vector
C = computeSKAT(X, Y, nPCs)

#   S : selected SNPs - k x 1 vector containing the selected indices
S = macarons(C, X, snp, k, D)

#   Save the result to a file
#np.save('result.npy', S)