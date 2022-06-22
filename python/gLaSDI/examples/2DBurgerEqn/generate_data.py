#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys
sys.path.append('../../src')
import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import spsolve
from scipy.sparse import spdiags
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import pickle
from time import time
import os
from sindy_utils import *
get_ipython().run_line_magic('matplotlib', 'inline')
np.set_printoptions(threshold=sys.maxsize)


# In[ ]:


plot_fig = False
save_data = True
Re = 10000
nx = 60
ny = nx
ic = 2  # initial condition, 1: Sine, 2: Gaussian
nt = 200
tstop = 1.0
dt = tstop/nt
t = np.linspace(0, tstop, nt+1)
nxy = (nx-2)*(ny-2)
dx = 1/(nx-1)
dy = 1/(ny-1)
maxitr = 10
tol = 1e-8

# parameters
amp_arr = np.array([0.7])
width_arr = np.array([0.9])

# amp_arr = np.linspace(0.7, 0.9, 5)
# width_arr = np.linspace(0.9, 1.1, 5)


# In[ ]:


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

# full indicies, free indicies, fixed indicies
[xv,yv]=np.meshgrid(np.linspace(xmin,xmax,nx),np.linspace(ymin,ymax,ny),indexing='xy')
x=xv.flatten()
y=yv.flatten()

multi_index_i,multi_index_j=np.meshgrid(np.arange(nx),np.arange(ny),indexing='xy')
full_multi_index=(multi_index_j.flatten(),multi_index_i.flatten())
free_multi_index=(multi_index_j[1:-1,1:-1].flatten(),multi_index_i[1:-1,1:-1].flatten())
x0_multi_index=(multi_index_j[1:-1,0].flatten(),multi_index_i[1:-1,0].flatten())
x1_multi_index=(multi_index_j[1:-1,-1].flatten(),multi_index_i[1:-1,-1].flatten())
y0_multi_index=(multi_index_j[0,1:-1].flatten(),multi_index_i[0,1:-1].flatten())
y1_multi_index=(multi_index_j[-1,1:-1].flatten(),multi_index_i[-1,1:-1].flatten())

dims=(ny,nx)
full_raveled_indicies=np.ravel_multi_index(full_multi_index,dims)
free_raveled_indicies=np.ravel_multi_index(free_multi_index,dims)
x0_raveled_indicies=np.ravel_multi_index(x0_multi_index,dims)
x1_raveled_indicies=np.ravel_multi_index(x1_multi_index,dims)
x01_raveled_indicies=np.concatenate((x0_raveled_indicies,x1_raveled_indicies))
y0_raveled_indicies=np.ravel_multi_index(y0_multi_index,dims)
y1_raveled_indicies=np.ravel_multi_index(y1_multi_index,dims)
y01_raveled_indicies=np.concatenate((y0_raveled_indicies,y1_raveled_indicies))
fixed_raveled_indicies=np.setdiff1d(full_raveled_indicies,free_raveled_indicies)

# boundary one-hot vector
x0_one_hot=np.eye(nx-2)[0]
y0_one_hot=np.eye(ny-2)[0]
x1_one_hot=np.eye(nx-2)[-1]
y1_one_hot=np.eye(ny-2)[-1]

# inner grid
inner_multi_index_i,inner_multi_index_j=np.meshgrid(np.arange(nx-2),np.arange(ny-2),indexing='xy')
inner_x_multi_index=(np.concatenate((inner_multi_index_j[:,0].flatten(),inner_multi_index_j[:,-1].flatten())),
                     np.concatenate((inner_multi_index_i[:,0].flatten(),inner_multi_index_i[:,-1].flatten())))
inner_y_multi_index=(np.concatenate((inner_multi_index_j[0,:].flatten(),inner_multi_index_j[-1,:].flatten())),
                     np.concatenate((inner_multi_index_i[0,:].flatten(),inner_multi_index_i[-1,:].flatten())))

inner_dims=(ny-2,nx-2)
inner_x_raveled_indicies=np.ravel_multi_index(inner_x_multi_index,inner_dims)
inner_y_raveled_indicies=np.ravel_multi_index(inner_y_multi_index,inner_dims)


# first order derivative
# central
Mcb=sp.diags([np.zeros(nx-2),-np.ones(nx-2),np.ones(nx-2)],[0,-1,1],(nx-2,nx-2))
Mc=sp.kron(sp.eye(ny-2),Mcb,format="csr")

Ib=sp.eye(nx-2)
Nc=sp.kron(sp.diags([np.zeros(ny-2),-np.ones(ny-2),np.ones(ny-2)],[0,-1,1],(ny-2,ny-2)),Ib,format="csr")

# forward
Mfb=sp.diags([-np.ones(nx-2),np.ones(nx-2)],[0,1],(nx-2,nx-2))
Mf=sp.kron(sp.eye(ny-2),Mfb,format="csr")

Ib=sp.eye(nx-2)
Nf=sp.kron(sp.diags([-np.ones(ny-2),np.ones(ny-2)],[0,1],(ny-2,ny-2)),Ib,format="csr")

# backward
Mbb=sp.diags([np.ones(nx-2),-np.ones(nx-2)],[0,-1],(nx-2,nx-2))
Mb=sp.kron(sp.eye(ny-2),Mbb,format="csr")

Ib=sp.eye(nx-2)
Nb=sp.kron(sp.diags([np.ones(ny-2),-np.ones(ny-2)],[0,-1],(ny-2,ny-2)),Ib,format="csr")

# laplacian operator
Dxb=sp.diags([-2*np.ones(nx-2),np.ones(nx-2),np.ones(nx-2)],[0,-1,1],(nx-2,nx-2))
Dx=sp.kron(sp.eye(ny-2),Dxb,format="csr")

Ib=sp.eye(nx-2)
Dy=sp.kron(sp.diags([-2*np.ones(ny-2),np.ones(ny-2),np.ones(ny-2)],[0,-1,1],(ny-2,ny-2)),Ib,format="csr")


# In[ ]:


def generate_dataset(ic, amp, width, plot_fig=False):
    timer = []
    timer.append(time())

    # compute u_full and v_full
    if ic == 1: # IC: sine
        zv=amp*np.sin(2*np.pi*xv)*np.sin(2*np.pi*yv)
        zv[np.nonzero(xv>0.5)]=0.0
        zv[np.nonzero(yv>0.5)]=0.0
    elif ic == 2: # IC: Gaussian
        zv = amp * np.exp(-((xv-x0)**2 + (yv-y0)**2) / width)
        z = zv.flatten()
    u0 = z.copy()
    v0 = z.copy()

    # plot IC
    if plot_fig:
        fig = plt.figure(figsize=(15,3))
        ax_u = fig.add_subplot(141)
        p_u=ax_u.pcolor(x.reshape(ny,nx), y.reshape(ny,nx), u0.reshape(ny,nx))
        cb_u=fig.colorbar(p_u, ax=ax_u)
        ax_u.set_xlabel('$x$')
        ax_u.set_ylabel('$y$')
        ax_u.set_title('$u$')
        
        ax_v = fig.add_subplot(143)
        p_v=ax_v.pcolor(x.reshape(ny,nx), y.reshape(ny,nx), v0.reshape(ny,nx))
        cb_v=fig.colorbar(p_v, ax=ax_v)
        ax_v.set_xlabel('$x$')
        ax_v.set_ylabel('$y$')
        ax_v.set_title('$v$')

    # boundary for first order derivative term
    Bdudx0_cur=np.kron(u0[x0_raveled_indicies],x0_one_hot)
    Bdudy0_cur=np.kron(y0_one_hot,u0[y0_raveled_indicies])
    Bdvdx0_cur=np.kron(v0[x0_raveled_indicies],x0_one_hot)
    Bdvdy0_cur=np.kron(y0_one_hot,v0[y0_raveled_indicies])
    Bdudx1_cur=np.kron(u0[x1_raveled_indicies],x1_one_hot)
    Bdudy1_cur=np.kron(y1_one_hot,u0[y1_raveled_indicies])
    Bdvdx1_cur=np.kron(v0[x1_raveled_indicies],x1_one_hot)
    Bdvdy1_cur=np.kron(y1_one_hot,v0[y1_raveled_indicies])

    # boundary for second order derivative term
    bxu_cur=np.zeros(nxy)
    byu_cur=np.zeros(nxy)
    bxv_cur=np.zeros(nxy)
    byv_cur=np.zeros(nxy)

    bxu_cur[inner_x_raveled_indicies]=u0[x01_raveled_indicies]
    byu_cur[inner_y_raveled_indicies]=u0[y01_raveled_indicies]
    bxv_cur[inner_x_raveled_indicies]=v0[x01_raveled_indicies]
    byv_cur[inner_y_raveled_indicies]=v0[y01_raveled_indicies]

    def r(u_free,v_free,u_free_prev,v_free_prev,Mu_free,Mv_free,Nu_free,Nv_free,
          Bdudx0_cur,Bdvdx0_cur,Bdudx1_cur,Bdvdx1_cur,Bdudy0_cur,Bdvdy0_cur,Bdudy1_cur,Bdvdy1_cur,
          bxu_cur,bxv_cur,byu_cur,byv_cur):

        f_u=(-1/dx*(u_free*(Mu_free - Bdudx0_cur))
        -1/dy*(v_free*(Nu_free - Bdudy0_cur))
        +1/(Re*dx**2)*(Dx.dot(u_free) + bxu_cur)
        +1/(Re*dy**2)*(Dy.dot(u_free) + byu_cur))
        
        f_v=(-1/dx*(u_free*(Mv_free - Bdvdx0_cur))
        -1/dy*(v_free*(Nv_free - Bdvdy0_cur))
        +1/(Re*dx**2)*(Dx.dot(v_free) + bxv_cur)
        +1/(Re*dy**2)*(Dy.dot(v_free) + byv_cur))

        r_u=u_free-u_free_prev-dt*f_u
        r_v=v_free-v_free_prev-dt*f_v

        return np.concatenate((r_u,r_v))



    def J(u_free,v_free,Mu_free,Mv_free,Nu_free,Nv_free,
          Bdudx0_cur,Bdvdx0_cur,Bdudx1_cur,Bdvdx1_cur,Bdudy0_cur,Bdvdy0_cur,Bdudy1_cur,Bdvdy1_cur):

        df_udu = (-1/dx*(sp.diags(Mu_free - Bdudx0_cur,0,(nxy,nxy),format="csr") 
                            + sp.diags(u_free,0,(nxy,nxy),format="csr").dot(Mb))
        -1/dy*sp.diags(v_free,0,(nxy,nxy),format="csr").dot(Nb)
        +1/(Re*dx**2)*Dx
        +1/(Re*dy**2)*Dy)
        df_udv = -1/dy*sp.diags(Nu_free - Bdudy0_cur,0,(nxy,nxy),format="csr")
        df_vdu = -1/dx*sp.diags(Mv_free - Bdvdx0_cur,0,(nxy,nxy),format="csr")
        df_vdv = (-1/dx*sp.diags(u_free,0,(nxy,nxy),format="csr").dot(Mb)
        -1/dy*(sp.diags(Nv_free - Bdvdy0_cur,0,(nxy,nxy),format="csr")
                   + sp.diags(v_free,0,(nxy,nxy),format="csr").dot(Nb))
        +1/(Re*dx**2)*Dx
        +1/(Re*dy**2)*Dy)

        return sp.bmat([[I-dt*df_udu,-dt*df_udv],[-dt*df_vdu,I-dt*df_vdv]],format='csr')

    # solution snapshot
    u_full=np.zeros(((nt+1),ny*nx))
    v_full=np.zeros(((nt+1),ny*nx))
    
    # solution + intermediate snapshot
    u_full_inter=np.array([])
    v_full_inter=np.array([])

    # IC
    u_full[0]=np.copy(u0)
    v_full[0]=np.copy(v0)
    u0_free=u0[free_raveled_indicies]
    v0_free=v0[free_raveled_indicies]

    for k in range(nt):
        u_free_prev=np.copy(u_full[k,free_raveled_indicies])
        v_free_prev=np.copy(v_full[k,free_raveled_indicies])

        u_free=np.copy(u_full[k,free_raveled_indicies])
        v_free=np.copy(v_full[k,free_raveled_indicies])

        Mu_free=Mb.dot(u_free)
        Mv_free=Mb.dot(v_free)
        Nu_free=Nb.dot(u_free)
        Nv_free=Nb.dot(v_free)

        residual=r(u_free,v_free,u_free_prev,v_free_prev,Mu_free,Mv_free,Nu_free,Nv_free,
                   Bdudx0_cur,Bdvdx0_cur,Bdudx1_cur,Bdvdx1_cur,Bdudy0_cur,Bdvdy0_cur,Bdudy1_cur,Bdvdy1_cur,
                   bxu_cur,bxv_cur,byu_cur,byv_cur)

        for itr in range(maxitr):
            jacobian=J(u_free,v_free,Mu_free,Mv_free,Nu_free,Nv_free,
                       Bdudx0_cur,Bdvdx0_cur,Bdudx1_cur,Bdvdx1_cur,Bdudy0_cur,Bdvdy0_cur,Bdudy1_cur,Bdvdy1_cur)

            delta_free=spsolve(jacobian, -residual)

            u_free+=delta_free[:nxy]
            v_free+=delta_free[nxy:]

            Mu_free=Mb.dot(u_free)
            Mv_free=Mb.dot(v_free)
            Nu_free=Nb.dot(u_free)
            Nv_free=Nb.dot(v_free)

            residual=r(u_free,v_free,u_free_prev,v_free_prev,Mu_free,Mv_free,Nu_free,Nv_free,
                       Bdudx0_cur,Bdvdx0_cur,Bdudx1_cur,Bdvdx1_cur,Bdudy0_cur,Bdvdy0_cur,Bdudy1_cur,Bdvdy1_cur,
                       bxu_cur,bxv_cur,byu_cur,byv_cur)
       
            # store itermediate values
            R = np.linalg.norm(residual)

            if R<tol:
                u_full[k+1,free_raveled_indicies]=np.copy(u_free)
                v_full[k+1,free_raveled_indicies]=np.copy(v_free)
                # BC from exact solution
                u_full[k+1,fixed_raveled_indicies]=np.copy(u0[fixed_raveled_indicies])
                v_full[k+1,fixed_raveled_indicies]=np.copy(v0[fixed_raveled_indicies])
                break

        if R>=tol:
            print("\n non converged after {}th iteration".format(maxitr))
            break     
    timer.append(time())
    time_fom = timer[1]-timer[0]
    
    # plot solution at final time step
    if plot_fig:
        ax_u = fig.add_subplot(142)
        ax_u = plt.gca()
        p_u=ax_u.pcolor(xv, yv, (u_full[-1]).reshape(ny,nx))
        cb_u=fig.colorbar(p_u, ax=ax_u)
        ax_u.set_xlabel('$x$')
        ax_u.set_ylabel('$y$')
        ax_u.set_title('$u$')

        ax_v = fig.add_subplot(144)
        ax_v = plt.gca()
        p_v=ax_v.pcolor(xv, yv, (v_full[-1]).reshape(ny,nx))
        cb_v=fig.colorbar(p_v, ax=ax_v)
        ax_v.set_xlabel('$x$')
        ax_v.set_ylabel('$y$')
        ax_v.set_title('$v$')
        plt.tight_layout()
        plt.show()  
    return u_full.reshape(-1,ny*nx), v_full.reshape(-1,ny*nx),time_fom


# ### Generate and save data in one file

# In[ ]:


# generate training set
data = []
param = []
time_fom = 0
for i in amp_arr:
    for j in width_arr:
        snapshot = {}
        param.append(np.array([i, j]))

        # compute u, v
        u, v, t_fom = generate_dataset(ic, i, j, plot_fig=plot_fig)
        print(f'amplitude: {i}, width: {j}, shape: {u.shape}, {v.shape}')
        snapshot['u'] = u
        snapshot['v'] = v
        snapshot['t'] = t
        snapshot['time_fom'] = t_fom
        time_fom += t_fom
        
        # compute du/dt, dv/dt
        dudt = derivative(u,tstop)
        dvdt = derivative(v,tstop)
        snapshot['du'] = dudt
        snapshot['dv'] = dvdt

        data.append(snapshot)

data_all = {}
data_all['data'] = data
data_all['param'] = param
print(f"time for computing x: {time_fom/len(data_all['data']):.2f} s")


# In[ ]:


nplot = 8
step_list = np.linspace(0,nt,nplot,dtype=int)
vmin = data[0]['u'].min()
vmax = data[0]['u'].max()

fig = plt.figure(figsize=(18,5))
for i,step in enumerate(step_list):
    ax_u = fig.add_subplot(2,nplot,i+1)
    ax_u = plt.gca()
    p_u=ax_u.pcolor(xv, yv, (data[0]['u'][step]).reshape(ny,nx))
    ax_u.set_title(f'$u$ - t: {step/nt*tstop:.1f} s',fontsize=16)
    
for i,step in enumerate(step_list):
    ax_u = fig.add_subplot(2,nplot,i+1+nplot)
    ax_u = plt.gca()
    p_u=ax_u.pcolor(xv, yv, (data[0]['v'][step]).reshape(ny,nx))
    ax_u.set_title(f'$v$ - t: {step/nt*tstop:.1f} s',fontsize=16)
plt.tight_layout()
# plt.savefig('./2dBurger_u_multisteps.png')


# In[ ]:


# save data
if save_data:
    num_case = len(data)
    data_path = '/usr/workspace/he10/data/2DBurgerEqn'
    if not os.path.exists(data_path):
        os.makedirs(data_path)
    if num_case > 1:
        pickle.dump(data_all, open(data_path+f"/local{num_case}_Re{Re}_tstop{tstop:.1f}_nt{nt}_nx{nx}.p", "wb"))
    else:
        pickle.dump(data_all, open(data_path+f"/local{num_case}_Re{Re}_A{amp_arr[0]:.2f}_W{width_arr[0]:.2f}_tstop{tstop:.1f}_nt{nt}_nx{nx}.p", "wb"))

