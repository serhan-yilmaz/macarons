import numpy as np
from util import corr

# Implementation of the Macarons (Improved)
# C - n x 1 score vector
# X - m x n genotype matrix
# snp - n x 2 SNP location matrix (row is chromosome, location)
# k - the number of SNPs to select
# D - threshold parameter
# Returns: a k x 1 vector containing the indices of the chosen SNPs
def macarons(C, X, snp, k, D, i_size = 1000, g_factor = 2):
    sel = []
    inds = np.argsort(-C, kind='stable')
    scores = C[inds]
    snp = snp[inds,:]
    N = i_size
    active_range = np.zeros(X.shape[1], dtype=bool)
    in_range = np.zeros(X.shape[1], dtype=bool)
    active_snp = snp[:N,:]
    active_range[:N] = True
    while(len(sel) < k):
        chosen = np.argmax(scores[:N+1]) # best SNP *will* be either in active, or the next SNP
        if (chosen < N): 
            # if the chosen SNP is in the active set, select it and
            # penalize the SNPs that are close (within D)
            sel.append(chosen)
            chr = snp[chosen,0]
            loc = snp[chosen,1]
            in_range[:N] = (active_snp[:,0] == chr) & (np.abs(active_snp[:,1] - loc) <= D)
            sqcorr = corr(X[:, inds[in_range]], X[:, inds[chosen]])**2
            t = scores[:N]
            t[in_range[:N]] = t[in_range[:N]]*(1-sqcorr)
            scores[:N] = t
        else: 
            # if the chosen SNP is not in the active set, expand the active set
            # and update the scores of the newly added SNPs
            N_next = min(N * g_factor, X.shape[1])
            active_range = False* active_range
            active_snp = snp[N:N_next,:]
            in_range = False * in_range
            for i in sel:
                chr = snp[i,0]
                loc = snp[i,1]
                in_range[N:N_next] = (active_snp[:,0] == chr) & (np.abs(active_snp[:,1] - loc) <= D)
                sqcorr = corr(X[:, inds[in_range]], X[:, inds[i]])**2
                t = scores[N:N_next]
                t[in_range[N:N_next]] = t[in_range[N:N_next]]*(1-sqcorr)
                scores[N:N_next] = t
            N = N_next
            active_range[:N] = True
            active_snp = snp[:N,:]
    # translate sel to original indices
    sel = np.array(inds[sel])
    return sel.astype(int)

# Original Macarons implementation
# C - n x 1 score vector
# X - m x n genotype matrix
# snp - n x 2 SNP location matrix (row is chromosome, location)
# k - the number of SNPs to select
# D - threshold parameter
# Returns: a k x 1 vector containing the indices of the chosen SNPs
def macarons_original(C, X, snp, k, D):
    sel = np.zeros(k)
    scores = C.copy()
    for i in range(k):
        chosen = np.argmax(scores)
        sel[i] = chosen
        chr = snp[chosen,0]
        loc = snp[chosen,1]
        in_range = (snp[:,0] == chr) & (np.abs(snp[:,1] - loc) <= D)
        sqcorr = corr(X[:, in_range], X[:, chosen])**2
        scores[in_range] = scores[in_range]*(1-sqcorr)
    return sel.astype(int)