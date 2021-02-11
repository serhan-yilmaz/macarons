import numpy as np

# Returns the top n principal componenents of the matrix X
def pca(X, n):
    Xn = (X - X.mean(0)) / X.std(0)
    w,v = np.linalg.eig(np.dot(Xn, Xn.T))
    inds = w.argsort()[::-1][:n]
    PC = v[:,inds]
    return PC

# Computes the SKAT scores for the features of X w.r.t Y 
# while correcting the population statification using the 
# k-first principal components
def computeSKAT(X, Y, k):
    # calculate the pc components
    n = X.shape[0]
    P = pca(X, k)


    cov = np.ones((n, k+1))
    cov[:,1::] = P

    est = np.dot(np.linalg.pinv(cov), Y)
    yp = np.dot(cov, est)

    r = Y - yp
    C = ((np.dot(X.T, r)) / n)**2
    return C

# Calculates the correlation between a vector vec
# and the columns of m. Returns a vector with
# dimension equal to the number of columns of m
def corr(vec, m):
    assert(vec.shape[0] == m.shape[0])

    if(m.ndim == 1 and vec.ndim > 1): 
        t = m
        m = vec
        vec = t

    vec = vec - vec.mean(0)
    m = m - m.mean(0)

    ssv = (vec**2).sum(0)
    ssm = (m**2).sum(0)

    return np.dot(vec, m) / (np.sqrt(ssv*ssm))