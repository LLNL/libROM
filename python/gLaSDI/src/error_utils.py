import numpy as np
from scipy import sparse as sp
import subprocess
from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve


def err_indicator(u_sim, params, data=None, err_type=1):
    """
    This function computes errors using a speciffied error indicator.
    inputs:
        data: dict, data of the evalution case
        err_type: int, types of error indicator
                1: max relative error (if test data is available)
                2: residual norm (mean), 1D Burger's eqn
                3: residual norm (mean), 2D Burger's eqn
                4: MFEM example 16: Time dependent heat conduction
                5: MFEM example 9: DG advection
    outputs:
        err: float, error
    """
    if err_type == 1:
        err = (np.linalg.norm(data - u_sim, axis=1) / np.linalg.norm(data, axis=1)*100).max()
    elif err_type == 2:
        res = []
        for k in range(u_sim.shape[0]-1):
            res.append(residual_1Dburger(u_sim[k,:], u_sim[k+1,:], params))
        err = np.stack(res).mean()
    elif err_type == 3:
        res = []
        for k in range(u_sim.shape[0]-1):
            res.append(residual_2Dburger(u_sim[k,:], u_sim[k+1,:], params))
        err = np.stack(res).mean()
    elif err_type == 4:
        u_file = params['pde']['u_file']
        np.savetxt(u_file, np.array([[u_sim.shape[0], u_sim.shape[1]]]), fmt='%d')
        with open(u_file, 'ab') as f:
            np.savetxt(f, u_sim, fmt='%.18f')
        err = residual_MFEMex16(params) / int(u_sim.shape[0]*params['pde']['res_ns']-1)
    elif err_type == 5:
        u_file = params['pde']['u_file']
        np.savetxt(u_file, np.array([[u_sim.shape[0], u_sim.shape[1]]]), fmt='%d')
        with open(u_file, 'ab') as f:
            np.savetxt(f, u_sim, fmt='%.18f')
        err = residual_MFEMex9(params) / int(u_sim.shape[0]*params['pde']['res_ns']-1)
    return err


def residual_1Dburger(u0, u1, params):
    """
    r = -u^{n} + u^{n+1} -dt*f(u^{n+1})
    """
    nx = params['pde']['nx']
    nt = params['pde']['nt']
    tstop = params['pde']['tstop']
    dx = 6 / (nx - 1)
    dt = tstop / nt 
    c = dt / dx

    idxn1 = np.zeros(nx,dtype='int')
    idxn1[1:] = np.arange(nx-1)
    idxn1[0] = nx-1
    
    f = c*(u1**2 - u1*u1[idxn1]) 
    r = -u0 + u1 + f
    return np.linalg.norm(r)


def residual_2Dburger(x_prev, x, params):
    Re = params['pde']['Re']
    nx = params['pde']['nx']
    ny = nx
    nt = params['pde']['nt']
    tstop = params['pde']['tstop']
    ic = params['pde']['ic'] # initial condition, 1: Sine, 2: Gaussian
    u_prev = x_prev[:nx*ny]
    u = x[:nx*ny]
    v_prev = x_prev[nx*ny:]
    v = x[nx*ny:]
    
    dt=tstop/nt
    t = np.linspace(0, tstop, nt+1)
    nxy = (nx-2)*(ny-2)
    dx = 1/(nx-1)
    dy = 1/(ny-1)

    if ic == 1: # sine
        xmin = 0
        xmax = 1
        ymin = 0
        ymax = 1
    elif ic == 2: # Gaussian
        xmin = -3
        xmax = 3
        ymin = -3
        ymax = 3
        x0 = 0 # Gaussian center
        y0 = 0 # Gaussian center
    else: 
        print('wrong values for IC!')
    I=sp.eye(nxy,format='csr')

    # full indices, free indices, fixed indices
    [xv,yv] = np.meshgrid(np.linspace(xmin,xmax,nx),np.linspace(ymin,ymax,ny),indexing='xy')
    x = xv.flatten()
    y = yv.flatten()

    multi_index_i,multi_index_j = np.meshgrid(np.arange(nx),np.arange(ny),indexing='xy')
    full_multi_index = (multi_index_j.flatten(),multi_index_i.flatten())
    free_multi_index = (multi_index_j[1:-1,1:-1].flatten(),multi_index_i[1:-1,1:-1].flatten())
    x0_multi_index = (multi_index_j[1:-1,0].flatten(),multi_index_i[1:-1,0].flatten())
    x1_multi_index = (multi_index_j[1:-1,-1].flatten(),multi_index_i[1:-1,-1].flatten())
    y0_multi_index = (multi_index_j[0,1:-1].flatten(),multi_index_i[0,1:-1].flatten())
    y1_multi_index = (multi_index_j[-1,1:-1].flatten(),multi_index_i[-1,1:-1].flatten())

    dims=(ny,nx)
    full_raveled_indices = np.ravel_multi_index(full_multi_index,dims)
    free_raveled_indices = np.ravel_multi_index(free_multi_index,dims)
    x0_raveled_indices = np.ravel_multi_index(x0_multi_index,dims)
    x1_raveled_indices = np.ravel_multi_index(x1_multi_index,dims)
    x01_raveled_indices = np.concatenate((x0_raveled_indices,x1_raveled_indices))
    y0_raveled_indices = np.ravel_multi_index(y0_multi_index,dims)
    y1_raveled_indices = np.ravel_multi_index(y1_multi_index,dims)
    y01_raveled_indices = np.concatenate((y0_raveled_indices,y1_raveled_indices))
    fixed_raveled_indices = np.setdiff1d(full_raveled_indices,free_raveled_indices)

    # boundary one-hot vector
    x0_one_hot = np.eye(nx-2)[0]
    y0_one_hot = np.eye(ny-2)[0]
    x1_one_hot = np.eye(nx-2)[-1]
    y1_one_hot = np.eye(ny-2)[-1]

    # inner grid
    inner_multi_index_i,inner_multi_index_j=np.meshgrid(np.arange(nx-2),np.arange(ny-2),indexing='xy')
    inner_x_multi_index = (np.concatenate((inner_multi_index_j[:,0].flatten(),inner_multi_index_j[:,-1].flatten())),
                         np.concatenate((inner_multi_index_i[:,0].flatten(),inner_multi_index_i[:,-1].flatten())))
    inner_y_multi_index = (np.concatenate((inner_multi_index_j[0,:].flatten(),inner_multi_index_j[-1,:].flatten())),
                         np.concatenate((inner_multi_index_i[0,:].flatten(),inner_multi_index_i[-1,:].flatten())))

    inner_dims = (ny-2,nx-2)
    inner_x_raveled_indices = np.ravel_multi_index(inner_x_multi_index,inner_dims)
    inner_y_raveled_indices = np.ravel_multi_index(inner_y_multi_index,inner_dims)


    # first order derivative
    # central
    Mcb = sp.diags([np.zeros(nx-2),-np.ones(nx-2),np.ones(nx-2)],[0,-1,1],(nx-2,nx-2))
    Mc = sp.kron(sp.eye(ny-2),Mcb,format="csr")

    Ib = sp.eye(nx-2)
    Nc = sp.kron(sp.diags([np.zeros(ny-2),-np.ones(ny-2),np.ones(ny-2)],[0,-1,1],(ny-2,ny-2)),Ib,format="csr")

    # forward
    Mfb = sp.diags([-np.ones(nx-2),np.ones(nx-2)],[0,1],(nx-2,nx-2))
    Mf = sp.kron(sp.eye(ny-2),Mfb,format="csr")

    Ib = sp.eye(nx-2)
    Nf = sp.kron(sp.diags([-np.ones(ny-2),np.ones(ny-2)],[0,1],(ny-2,ny-2)),Ib,format="csr")

    # backward
    Mbb = sp.diags([np.ones(nx-2),-np.ones(nx-2)],[0,-1],(nx-2,nx-2))
    Mb = sp.kron(sp.eye(ny-2),Mbb,format="csr")

    Ib = sp.eye(nx-2)
    Nb = sp.kron(sp.diags([np.ones(ny-2),-np.ones(ny-2)],[0,-1],(ny-2,ny-2)),Ib,format="csr")

    # laplacian operator
    Dxb = sp.diags([-2*np.ones(nx-2),np.ones(nx-2),np.ones(nx-2)],[0,-1,1],(nx-2,nx-2))
    Dx = sp.kron(sp.eye(ny-2),Dxb,format="csr")

    Ib = sp.eye(nx-2)
    Dy = sp.kron(sp.diags([-2*np.ones(ny-2),np.ones(ny-2),np.ones(ny-2)],[0,-1,1],(ny-2,ny-2)),Ib,format="csr")

    # Initial condition
    amp = params['pde']['param'][0]
    width = params['pde']['param'][1]
    if ic == 1: # IC: sine
        zv = amp*np.sin(2*np.pi*xv)*np.sin(2*np.pi*yv)
        zv[np.nonzero(xv>0.5)] = 0.0
        zv[np.nonzero(yv>0.5)] = 0.0
    elif ic == 2: # IC: Gaussian
        zv = amp * np.exp(-((xv-x0)**2 + (yv-y0)**2) / width)
        z = zv.flatten()
    u0 = z.copy()
    v0 = z.copy()
    
    
    # boundary for first order derivative term
    Bdudx0_cur=np.kron(u0[x0_raveled_indices],x0_one_hot)
    Bdudy0_cur=np.kron(y0_one_hot,u0[y0_raveled_indices])
    Bdvdx0_cur=np.kron(v0[x0_raveled_indices],x0_one_hot)
    Bdvdy0_cur=np.kron(y0_one_hot,v0[y0_raveled_indices])
    Bdudx1_cur=np.kron(u0[x1_raveled_indices],x1_one_hot)
    Bdudy1_cur=np.kron(y1_one_hot,u0[y1_raveled_indices])
    Bdvdx1_cur=np.kron(v0[x1_raveled_indices],x1_one_hot)
    Bdvdy1_cur=np.kron(y1_one_hot,v0[y1_raveled_indices])

    # boundary for second order derivative term
    bxu_cur=np.zeros(nxy)
    byu_cur=np.zeros(nxy)
    bxv_cur=np.zeros(nxy)
    byv_cur=np.zeros(nxy)

    bxu_cur[inner_x_raveled_indices]=u0[x01_raveled_indices]
    byu_cur[inner_y_raveled_indices]=u0[y01_raveled_indices]
    bxv_cur[inner_x_raveled_indices]=v0[x01_raveled_indices]
    byv_cur[inner_y_raveled_indices]=v0[y01_raveled_indices]
    
    u_free_prev=np.copy(u_prev[free_raveled_indices])
    v_free_prev=np.copy(v_prev[free_raveled_indices])

    u_free=np.copy(u[free_raveled_indices])
    v_free=np.copy(v[free_raveled_indices])

    Mu_free=Mb.dot(u_free)
    Mv_free=Mb.dot(v_free)
    Nu_free=Nb.dot(u_free)
    Nv_free=Nb.dot(v_free)


    f_u=(-1/dx*(u_free*(Mu_free - Bdudx0_cur))
    -1/dy*(v_free*(Nu_free - Bdudy0_cur))
    +1/(Re*dx**2)*(Dx.dot(u_free) + bxu_cur)
    +1/(Re*dy**2)*(Dy.dot(u_free) + byu_cur))

    f_v=(-1/dx*(u_free*(Mv_free - Bdvdx0_cur))
    -1/dy*(v_free*(Nv_free - Bdvdy0_cur))
    +1/(Re*dx**2)*(Dx.dot(v_free) + bxv_cur)
    +1/(Re*dy**2)*(Dy.dot(v_free) + byv_cur))

    r_u = u_free - u_free_prev-dt*f_u
    r_v = v_free - v_free_prev-dt*f_v
        
    return np.linalg.norm(r_u)+np.linalg.norm(r_v)


def residual_MFEMex16(params):  
    """
    This function calculates the residual error of the 
    heat conduction problem (MFEM example 16).
    The initial condition is parameterized by w0 and h0.
    Their values are stored in `params['pde']['param']`
    """
    subprocess.call([params['pde']['exe_file'],
                     '-m', params['pde']['m_file'],
                     '-uf', params['pde']['u_file'],
                     '-rf', params['pde']['res_file'],
                     '-r', str(params['pde']['rl']),
                     '-o', str(params['pde']['order']),
                     '-s', str(params['pde']['ODEsolver']),
                     '-tf', str(params['pde']['tstop']),
                     '-dt', str(params['pde']['dt']),
                     '-w0', str(params['pde']['param'][0]),
                     '-h0', str(params['pde']['param'][1]),
                     '-x1', str(params['pde']['x1']),
                     '-x2', str(params['pde']['x2']),
                     '-res_ns', str(params['pde']['res_ns']),
                     '-Tmax_iter', str(params['pde']['Tmax_iter']),
                     '-no-vis',
                     '-no-visit',
                     '-res'])
    res = np.loadtxt(params['pde']['res_file'])
    return res
    
    
def residual_MFEMex9(params):
    """
    This function calculates the residual error of the 
    radial advection problem (MFEM example 9).
    The initial condition is parameterized by w1 and w2.
    Their values are stored in `params['pde']['param']`
    """
    subprocess.call([params['pde']['exe_file'],
                     '-m', params['pde']['m_file'],
                     '-uf', params['pde']['u_file'],
                     '-rf', params['pde']['res_file'],
                     '-p', str(params['pde']['prob']),
                     '-r', str(params['pde']['rl']),
                     '-o', str(params['pde']['order']),
                     '-tf', str(params['pde']['tstop']),
                     '-dt', str(params['pde']['dt']),
                     '-w1', str(params['pde']['param'][0]),
                     '-w2', str(params['pde']['param'][1]),
                     '-res_ns', str(params['pde']['res_ns']),
                     '-Mmax_iter', str(params['pde']['Mmax_iter']),
                     '-no-vis',
                     '-no-visit',
                     '-res'])
    res = np.loadtxt(params['pde']['res_file'])
    return res