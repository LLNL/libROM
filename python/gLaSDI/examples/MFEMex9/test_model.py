#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys
sys.path.append("../../src")
import os
import numpy as np
import pickle
from autoencoder import full_network
from training import create_feed_dictionary, create_feed_dictionary2, eval_model, max_err_heatmap
from sindy_utils import *
from error_utils import *
import tensorflow.compat.v1 as tf
tf.compat.v1.disable_v2_behavior()
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from time import time
from copy import deepcopy
import subprocess as sp
get_ipython().run_line_magic('matplotlib', 'inline')
import random
from sklearn.linear_model import LinearRegression
import copy

def get_cmap(n, name='tab20'):
    return plt.cm.get_cmap(name, n)
cmap = get_cmap(10)


# In[ ]:


def get_gpu_memory():
  _output_to_list = lambda x: x.decode('ascii').split('\n')[:-1]

  ACCEPTABLE_AVAILABLE_MEMORY = 1024
  COMMAND = "nvidia-smi --query-gpu=memory.free --format=csv"
  memory_free_info = _output_to_list(sp.check_output(COMMAND.split()))[1:]
  memory_free_values = [int(x.split()[0]) for i, x in enumerate(memory_free_info)]
  return memory_free_values

device_list = tf.config.list_physical_devices('GPU')
free_mem = get_gpu_memory()
for i,gpu in enumerate(device_list):
    print(f'{gpu}: free memory: {free_mem[i]}')


# In[ ]:


# specify which GPU to use
config = tf.ConfigProto(log_device_placement=False, gpu_options=tf.GPUOptions(allow_growth=True,
                                                                              visible_device_list='1'))


# ## Load the Trained gLaSDI

# In[ ]:


data_path = os.getcwd() + '/fig/test/'
save_name = 'ex9_2022_05_21_09_49_43'
params = pickle.load(open(data_path + save_name + '_params.pkl', 'rb'))
params['save_name'] = data_path + save_name
params['config'] = config


# ## Evaluation by One Parameter Case

# In[ ]:


knn = 1
w1 = 1.5
w2 = 2
tstop = 3
t_test = tstop
test_data = pickle.load(open(f"/usr/workspace/he10/data/MFEMex9/local1_p1{w1:.1f}_p2{w2:.1f}_tstop{tstop:.1f}.p", "rb"))
nt = test_data['data'][0]['x'].shape[0]
t = np.linspace(0,t_test,nt)

# discretization coordinates
vert = np.loadtxt('/usr/workspace/he10/data/MFEMex9/vertex.txt')
triang = tri.Triangulation(vert[:,0], vert[:,1])

test_data_x = test_data['data'][0]['x']
test_data_dx = test_data['data'][0]['dx']


# In[ ]:


u_decoder,du_decoder,u_sim,du_sim,z_encoder,dz_encoder,z_sim,dz_sim,idx,timer_rom = eval_model(test_data['data'][0], params,
                                                                                               test_data['param'][0], knn=knn,
                                                                                               calc_dz=True, calc_du=True)
u_decoder = u_decoder.squeeze()
print(z_sim.shape, u_sim.shape, u_decoder.shape)


# In[ ]:


# max relative error
err_decoder = np.linalg.norm(test_data_x - u_decoder, axis=1) / np.linalg.norm(test_data_x, axis=1)*100
err_sindy = np.linalg.norm(test_data_x - u_sim, axis=1) / np.linalg.norm(test_data_x, axis=1)*100
print(f'max autoencoder error: {err_decoder.max():.2f} %')
print(f'max sindy-decoder error: {err_sindy.max():.2f} %')


# In[ ]:


nstep = 8
step_list = np.linspace(0,nt-1,nstep).astype(int)
vmin_x = test_data_x.min()
vmax_x = test_data_x.max()
vmin_dx = test_data_dx.min()/10
vmax_dx = test_data_dx.max()/10

fig = plt.figure(figsize=(18,10))
for i,step in enumerate(step_list):
    ax = fig.add_subplot(4,nstep,i+1)
    ax.tripcolor(triang, test_data_x[step], shading='gouraud', vmin=vmin_x, vmax=vmax_x)
    ax.set_title(f'u - step: {step}', fontsize=16)
    ax.set_aspect('equal')
    plt.axis('off')
plt.tight_layout()

for i,step in enumerate(step_list):
    ax = fig.add_subplot(4,nstep,i+1+nstep)
    ax.tripcolor(triang, u_sim[step], shading='gouraud', vmin=vmin_x, vmax=vmax_x)
    ax.set_title(r'$u_{pred}$' + f' - step: {step}', fontsize=16)
    ax.set_aspect('equal')
    plt.axis('off')
plt.tight_layout()
# plt.savefig(data_path + f"inter_pred_multisteps.png")


# In[ ]:


plt.rcParams.update({"font.size": 24,
                     "font.family": "sans-serif"}) # fontsize for figures

fig1 = plt.figure(figsize=(12,5))
line_type = ['-','-*','-.','-^','-s']
idx = np.arange(0,t.size,10)
ax = fig1.add_subplot(121)
for i in range(z_encoder.shape[1]):
    ax.plot(t, z_encoder[:,i], '-', lw=2, c=cmap(i))
    ax.plot(t[idx], z_sim[idx,i], '--o', lw=2, markersize=5, c=cmap(i))
ax.set_xlabel('Time')
ax.set_ylabel('z')
ax.set_xticks(np.linspace(0,t.max(),5))
ax.set_ylim(z_sim.min()*1.1,z_sim.max()+1.2)
ax.tick_params(axis='both', labelsize=24)
ax.legend(['Encoder', 'DI'], loc='upper right', frameon=False, fontsize=24)
ax.set_xlim(t.min(),t.max())

ax = fig1.add_subplot(122)
for i in range(z_sim.shape[1]):
    ax.plot(t, dz_encoder[:,i], '-', linewidth=2, c=cmap(i))
    ax.plot(t[idx], dz_sim[idx,i], '--o', linewidth=2, markersize=5, c=cmap(i))
ax.set_xlabel('Time')
ax.set_ylabel('dz/dt')
ax.set_xticks(np.linspace(0,t.max(),5))
ax.set_xlim(0,t.max())
ax.set_ylim(dz_sim.min()*1.1,dz_sim.max()*1.3)
ax.tick_params(axis='both', labelsize=24)
ax.legend(['Encoder', 'DI'], loc='upper right', frameon=False, fontsize=24)

plt.tight_layout()
# plt.savefig(data_path + f"advection_latent_dynamics.png",bbox_inches='tight')


# ## Evaluation by the Prescribed Parameter Space

# In[ ]:


knn = 4
res_name = f'mean'

na = 21
nw = na
p1_test = np.linspace(1.5, 2, na)
p2_test = np.linspace(2, 2.5, nw)
num_case = p1_test.size * p2_test.size
max_err = np.zeros([len(p1_test), len(p2_test)])
res_norm = np.zeros([len(p1_test), len(p2_test)])
sindy_idx = np.zeros([len(p1_test), len(p2_test)])
test_data_all = pickle.load(open(f"/usr/workspace/he10/data/MFEMex9/local{num_case}_tstop{tstop:.1f}.p", "rb"))

speed_up = 0
count = 0
timer_rom = np.zeros(4)
start_time = time()
for i,a in enumerate(p1_test):
    for j,w in enumerate(p2_test):
        print(f"{count}/{num_case}: {test_data_all['param'][count]}")
        test_data = {}
        test_data['data'] = [deepcopy(test_data_all['data'][count])]
        test_data['param'] = [deepcopy(test_data_all['param'][count])]
        
        _,_,u_sim,_,_,_,_,_,idx,t_rom = eval_model(test_data['data'][0], params, 
                                                   test_data['param'][0], knn=knn)
        timer_rom += t_rom
        sindy_idx[i,j] = idx+1
        
        # Max error of all time steps
        max_err[i,j] = (np.linalg.norm(test_data['data'][0]['x'] - u_sim, axis=1) \
                        / np.linalg.norm(test_data['data'][0]['x'], axis=1)*100).max()
        
        # residual norm
        res_norm[i,j] = err_indicator(u_sim, params, err_type=params['err_type'])
        count += 1

end_time = time()
time_rom = timer_rom[1:].sum()/num_case # from Step 2 to 4
time_sim = 3.8 # seconds
speed_up = time_sim / time_rom
print(f'Time taken: {end_time-start_time:.2f} s, {(end_time-start_time)/60:.2f} mins')
print(f'Max relative error: {max_err.max() :.2f} %')
print(f'Average speed up: {speed_up:.2f}')


# In[ ]:


a_grid, w_grid = np.meshgrid(p1_test, p2_test)
param_list = np.hstack([a_grid.flatten().reshape(-1,1), w_grid.flatten().reshape(-1,1)])
a_grid, w_grid = np.meshgrid(np.arange(p1_test.size), np.arange(p2_test.size))
idx_list = np.hstack([a_grid.flatten().reshape(-1,1), w_grid.flatten().reshape(-1,1)])

idx_param = []
for i,ip in enumerate(params['param']):
    idx = np.argmin(np.linalg.norm(param_list-ip, axis=1))
    idx_param.append((idx, np.array([param_list[idx,0], param_list[idx,1]])))


# In[ ]:


max_err_heatmap(max_err, sindy_idx, params, p1_test, p2_test, data_path, idx_list, idx_param,
                xlabel='w2', ylabel='w1', dtype='float')


# In[ ]:


# the residual-based errors are scaled by the "scale" parameter
max_err_heatmap(res_norm, sindy_idx, params, p1_test, p2_test, data_path, idx_list, idx_param,
                xlabel='w2', ylabel='w1', label='Residual Norm', dtype='float', scale=100)

