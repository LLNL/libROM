import numpy as np
from scipy.special import binom
from scipy.integrate import odeint
from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve

def library_size(n, poly_order, use_sine=False, use_cosine=False, include_constant=True):
    l = 0
    for k in range(poly_order+1):
        l += int(binom(n+k-1,k))
    if use_sine:
        l += n
    if use_cosine:
        l += n
    if not include_constant:
        l -= 1
    return l


def sindy_library(X, poly_order, include_sine=False, include_cosine=False):
    m,n = X.shape
    l = library_size(n, poly_order, include_sine, include_cosine, True)
    library = np.ones((m,l))
    index = 1

    if poly_order > 0: 
        for i in range(n):
            library[:,index] = X[:,i]
            index += 1
        
    if poly_order > 1:
        for i in range(n):
            for j in range(i,n):
                library[:,index] = X[:,i]*X[:,j]
                index += 1

    if poly_order > 2:
        for i in range(n):
            for j in range(i,n):
                for k in range(j,n):
                    library[:,index] = X[:,i]*X[:,j]*X[:,k]
                    index += 1

    if poly_order > 3:
        for i in range(n):
            for j in range(i,n):
                for k in range(j,n):
                    for q in range(k,n):
                        library[:,index] = X[:,i]*X[:,j]*X[:,k]*X[:,q]
                        index += 1
                    
    if poly_order > 4:
        for i in range(n):
            for j in range(i,n):
                for k in range(j,n):
                    for q in range(k,n):
                        for r in range(q,n):
                            library[:,index] = X[:,i]*X[:,j]*X[:,k]*X[:,q]*X[:,r]
                            index += 1

    if include_sine:
        for i in range(n):
            library[:,index] = np.sin(X[:,i])
            index += 1
            
    if include_cosine:
        for i in range(n):
            library[:,index] = np.cos(X[:,i])
            index += 1
            
    return library


def sindy_fit(RHS, LHS, coefficient_threshold):
    m,n = LHS.shape
    Xi = np.linalg.lstsq(RHS,LHS, rcond=None)[0]
    
    for k in range(10):
        small_inds = (np.abs(Xi) < coefficient_threshold)
        Xi[small_inds] = 0
        for i in range(n):
            big_inds = ~small_inds[:,i]
            if np.where(big_inds)[0].size == 0:
                continue
            Xi[big_inds,i] = np.linalg.lstsq(RHS[:,big_inds], LHS[:,i], rcond=None)[0]
    return Xi


def sindy_simulate(x0, t, Xi, poly_order, include_sine, include_cosine=False):
    m = t.size
    n = x0.size
    f = lambda x,t : np.dot(sindy_library(np.array(x).reshape((1,n)), poly_order, include_sine, include_cosine), Xi).reshape((n,))

    x = odeint(f, x0, t)
    return x


### Compact Finite Difference Method ###
def D_Lele(N,h):
    d=[-1, 0, 1];
    B1=3/8*np.ones(N)
    B2=3/8*np.ones(N)
    B1[-2]=3
    B1[-3]=1/4
    B1[-4]=1/3
    B1[0]=1/4
    B1[1]=1/3
    B2[1]=3
    B2[2]=1/4
    B2[3]=1/3
    B2[-1]=1/4
    B2[-2]=1/3
    A=spdiags([B1, np.ones(N), B2],d, N,N)
    
    alf=25/32/h
    bet=1/20/h
    gam=-1/480/h
    d=np.arange(-3,4)
    
    # d=-3
    B1=-gam*np.ones(N)
    B1[-4]=1/6/h
    B1[-5]=0
    B1[-6]=0
    # d=-2
    B2=-bet*np.ones(N)
    B2[0]=-1/36/h
    B2[-3]=-3/2/h
    B2[-4]=0
    B2[-5]=-1/36/h
    # d = -1
    B3=-alf*np.ones(N)
    B3[0]=-3/4/h
    B3[1]=-7/9/h
    B3[-2]=-3/2/h
    B3[-3]=-3/4/h
    B3[-4]=-7/9/h
    # d = 0
    B4=np.zeros(N)
    B4[0]=-17/6/h
    B4[-1]=17/6/h
    # d = 1
    B5=alf*np.ones(N)
    B5[1]=3/2/h
    B5[2]=3/4/h
    B5[3]=7/9/h
    B5[-1]=3/4/h
    B5[-2]=7/9/h
    # d = 2
    B6=bet*np.ones(N)
    B6[2]=3/2/h
    B6[3]=0
    B6[4]=1/36/h
    B6[-1]=1/36/h
    # d = 3
    B7=gam*np.ones(N)
    B7[3]=-1/6/h
    B7[4]=0
    B7[5]=0
    B=spdiags([B1, B2, B3, B4, B5, B6, B7],d,N,N)
    return spsolve(A.tocsc(),B.tocsc())


def derivative(x, tstop):
    dxdt = np.empty(x.shape)
    D = D_Lele(x.shape[0], tstop/x.shape[0])   
    for i in range(x.shape[1]):
        dxdt[:,i] = np.dot(D.toarray(), x[:,i])
    del D
    return dxdt