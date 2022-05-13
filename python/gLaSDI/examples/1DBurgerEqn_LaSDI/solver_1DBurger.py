import numpy as np
from scipy import sparse
from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve
from time import time

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
    
    
def kth_diag_indices(a, k):
    rows, cols = np.diag_indices_from(a)
    if k < 0:
        return rows[-k:], cols[:k]
    elif k > 0:
        return rows[:-k], cols[k:]
    else:
        return rows, cols
    
def sine_wave(amp, width):
    
    u0 = np.zeros(nx)
    u0[1:int(width/dx+1)] =amp/2*(np.sin(2*np.pi/(x[int(width/dx+1)]-x[1])*x[1:int(width/dx+1)]-np.pi/2)+1)
    u0[-1] = u0[0]
    
    return u0


def gaussian(amp, width, x):
    
    u0 = amp*np.exp(-(x-0.0)**2/(2*width**2))
    u0[-1] = u0[0]
    
    return u0

def residual(un, uw, c, idxn1):
    
    # r = -u^{n} + u^{n+1} -dt*f(u^{n+1})
    
    f = c*(uw**2 - uw*uw[idxn1]) 
    
    r = -un + uw + f
    
    return r

def jacobian(u, c, idxn1, nx):
    diag_comp = 1.0 + c*(2*u - u[idxn1])
    subdiag_comp = np.ones(nx-1)
    subdiag_comp[:-1] = -c*u[1:]
        
    data = np.array([diag_comp, subdiag_comp])
    J = spdiags(data,[0,-1],nx-1,nx-1,format='csr')
    J[0,-1] = -c*u[0]
    
    return J

def solve(u0, maxk, convergence_threshold, nt, nx, idxn1, c):
    u = np.zeros((nt+1,nx))
    u_inter=np.array([])
    u[0] = u0
    u_inter=np.append(u_inter,u0[:-1])
    I = sparse.eye(nx,format='csr')
    for n in range(nt): 
        uw = u[n,:-1].copy()
        r = residual(u[n,:-1], uw, c, idxn1)
        
        for k in range(maxk):
            J = jacobian(uw, c, idxn1, nx)
            duw = spsolve(J, -r)
            uw = uw + duw
            r = residual(u[n,:-1], uw, c, idxn1)
            u_inter = np.append(u_inter, uw)

            rel_residual = np.linalg.norm(r)/np.linalg.norm(u[n,:-1])
            if rel_residual < convergence_threshold:
                u[n+1,:-1] = uw.copy()
                u[n+1,-1] = u[n+1,0]
                break
    
    return u,u_inter.reshape((-1,nx-1))

 
def gen_data(amp_arr, width_arr, nx=1001, nt=1000):
    maxk = 10
    convergence_threshold = 1.0e-8
    x = np.linspace(-3, 3, nx)
    dx = 6 / (nx - 1)
    tstop = 2
    t = np.linspace(0, tstop, nt)
    dt = tstop / nt 
    c = dt/dx
    degree = 1 
    thres = 1

    idxn1 = np.zeros(nx-1,dtype='int')
    idxn1[1:] = np.arange(nx-2)
    idxn1[0] = nx-2
    timer = []
    
    # compute x
    timer.append(time())
    num_amp = amp_arr.shape[0]
    num_width = width_arr.shape[0]
    soln = []
    for i in range(num_amp):
        for j in range(num_width):
            u0 = gaussian(amp_arr[i],width_arr[j],x)
#             u0=sine_wave(amp_arr[i],width_arr[j])
            u,_ = solve(u0, maxk, convergence_threshold, nt, nx, idxn1, c)
            soln.append(u)
    soln = np.vstack(soln)    
    
    # compute dx/dt
    timer.append(time())
    data_x = soln.reshape(-1,nt+1,nx)
    dxdt_full = []
    t_full = []
    for num in range(len(data_x)):
        X = data_x[num] # exclude the last time step
        dxdt = np.empty(X.shape)
        D = D_Lele(X.shape[0],tstop/X.shape[0])   
        for i in range(X.shape[1]):
            dxdt[:,i] = np.dot(D.toarray(), X[:,i])
        dxdt_full.extend(dxdt)
        t_full.append(np.linspace(0,tstop,nt+1).reshape(-1,1))
    del D
    dxdt = np.array(dxdt_full)
    t = np.vstack(t_full)
    timer.append(time())
    
    time_fom = [timer[2]-timer[1], timer[1]-timer[0]]
    data = {}
    data['t'] = t
    data['x'] = soln
    data['dx'] = dxdt
    data['time_fom'] = time_fom
    return data


def get_data_1DBurger(amp, width):
    data = []
    param = []
    for i in amp:
        for j in width:
            tmp = gen_data(np.array([i]), np.array([j]))
            data.append(tmp)
            param.append(np.array([i, j]))
            
    data_all = {}
    data_all['data'] = data
    data_all['param'] = param

    return data_all