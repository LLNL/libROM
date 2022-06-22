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
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from time import time
get_ipython().run_line_magic('matplotlib', 'inline')
from copy import deepcopy
import subprocess as sp
from error_utils import residual_2Dburger
from sklearn.linear_model import LinearRegression
from matplotlib.ticker import FormatStrFormatter

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
save_name = 'burger_2022_05_21_09_49_44'
params = pickle.load(open(data_path + save_name + '_params.pkl', 'rb'))
params['save_name'] = data_path + save_name
params['config'] = config


# ## Evaluation by One Parameter Case

# In[ ]:


def process_data(data, vel, nt):
    # select component
    if vel == 1:
        data['data'][0]['x'] = data['data'][0].pop('u')
        data['data'][0]['dx'] = data['data'][0].pop('du')
        data['data'][0].pop('v')
        data['data'][0].pop('dv')
    elif vel == 2:
        data['data'][0]['x'] = data['data'][0].pop('v')
        data['data'][0]['dx'] = data['data'][0].pop('dv')
        data['data'][0].pop('u')
        data['data'][0].pop('du')
    elif vel == 3:
        data['data'][0]['x'] = np.hstack((data['data'][0]['u'], data['data'][0]['v']))
        data['data'][0]['dx'] = np.hstack((data['data'][0]['du'], data['data'][0]['dv']))
        
    # select time steps
    data['data'][0]['x'] = data['data'][0]['x'][:nt+1]
    data['data'][0]['dx'] = data['data'][0]['dx'][:nt+1]
    data_x = np.copy(data['data'][0]['x'])
    data_dx = np.copy(data['data'][0]['dx'])
    return data, data_x, data_dx


# In[ ]:


Re = params['pde']['Re']
nx = params['pde']['nx']
ny = nx
nt = params['pde']['nt']
tstop = params['pde']['tstop']
ic = params['pde']['ic']
t_test = tstop
vel = 3 # 1: u, 2: v, 3: u and v
knn = 1
amp_arr = np.array([0.7])
width_arr = np.array([0.9])
test_data = pickle.load(open(f"/usr/workspace/he10/data/2DBurgerEqn/local1_Re{Re}_A{amp_arr[0]:.2f}_W{width_arr[0]:.2f}_tstop{tstop:.1f}_nt{nt}_nx{nx}.p", "rb"))


# In[ ]:


nt_test = int(t_test/tstop*nt)
t = np.linspace(0,t_test,nt_test+1)
test_data, test_data_x, test_data_dx = process_data(test_data, vel, nt_test)

u_decoder,du_decoder,u_sim,du_sim,z_encoder,dz_encoder,z_sim,dz_sim,idx,timer_rom = eval_model(test_data['data'][0], params,
                                                                                               test_data['param'][0], knn=knn,
                                                                                               calc_dz=True, calc_du=True)
u_decoder = u_decoder.squeeze()
time_rom = timer_rom[1:].sum()
print(z_sim.shape, u_sim.shape, u_decoder.shape)
print(f'time: {time_rom:.2f} s')


# In[ ]:


# max relative error
err_decoder = np.linalg.norm(test_data_x - u_decoder, axis=1) / np.linalg.norm(test_data_x, axis=1)*100
err_sindy = np.linalg.norm(test_data_x - u_sim, axis=1) / np.linalg.norm(test_data_x, axis=1)*100
print(f'max autoencoder error: {err_decoder.max():.2f} %')
print(f'max sindy-decoder error: {err_sindy.max():.2f} %')


# In[ ]:


step_list = np.linspace(0,nt_test,10).astype(int)
fig = plt.figure(figsize=(18,3))
for i,step in enumerate(step_list):
    ax = fig.add_subplot(1,10,i+1)
    ax.imshow(test_data_x[step,:nx*ny].reshape(ny,nx))
    ax.set_title(f'u - step: {step}')
plt.tight_layout()

fig = plt.figure(figsize=(18,3))
for i,step in enumerate(step_list):
    ax = fig.add_subplot(1,10,i+1)
    ax.imshow(u_sim[step,:nx*ny].reshape(ny,nx))
    ax.set_title(f'u_pred - step: {step}')
plt.tight_layout()


# In[ ]:


step_list = np.linspace(0,nt_test,10).astype(int)
fig = plt.figure(figsize=(18,3))
for i,step in enumerate(step_list):
    ax = fig.add_subplot(1,10,i+1)
    ax.imshow(test_data_x[step,nx*ny:].reshape(ny,nx))
    ax.set_title(f'v - step: {step}')
plt.tight_layout()

fig = plt.figure(figsize=(18,3))
for i,step in enumerate(step_list):
    ax = fig.add_subplot(1,10,i+1)
    ax.imshow(u_sim[step,nx*ny:].reshape(ny,nx))
    ax.set_title(f'v_pred - step: {step}')
plt.tight_layout()


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
ax.set_ylim(z_sim.min()*1.1,z_sim.max()*2)
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
ax.set_ylim(dz_sim.min()*1.1,dz_sim.max()*2.2)
ax.tick_params(axis='both', labelsize=24)
ax.legend(['Encoder', 'DI'], loc='upper right', frameon=False, fontsize=24)

plt.tight_layout()
# plt.savefig(data_path + f"2Dburger_latent_dynamics.png",bbox_inches='tight')


# ## Evaluation by the Prescribed Parameter Space

# In[ ]:


knn = 3
res_name = f'mean'
nt_test = int(t_test/tstop*nt)
t = np.linspace(0,t_test,nt_test+1)
vel = 3 # 1: u, 2: v, 3: u and v

amp_test = params['test_param'][:,0]
width_test = params['test_param'][:,1]
amp_size = amp_test.size
width_size = width_test.size
num_case = amp_size * width_size
max_err = np.zeros([amp_size, width_size])
res_norm = np.zeros([amp_size, width_size])
sindy_idx = np.zeros([amp_size, width_size])
test_data_all = pickle.load(open(f"/usr/workspace/he10/data/2DBurgerEqn/local{num_case}_Re{Re}_tstop{tstop:.1f}_nt{nt}_nx{nx}.p", "rb"))

speed_up = 0
count = 0
timer_rom = np.zeros(4)
start_time = time()
for i,a in enumerate(amp_test):
    for j,w in enumerate(width_test):
        print(f"{count+1}/{num_case}: {test_data_all['param'][count]}")
        test_data = {}
        test_data['data'] = [deepcopy(test_data_all['data'][count])]
        test_data['param'] = [deepcopy(test_data_all['param'][count])]
        test_data, test_data_x,_ = process_data(test_data, vel, nt_test)
        _,_,u_sim,_,_,_,_,_,idx,t_rom = eval_model(test_data['data'][0], params, 
                                                   test_data['param'][0], knn=knn)
        timer_rom += t_rom
        sindy_idx[i,j] = idx+1
        
        # max error of all time steps
        max_err[i,j] = (np.linalg.norm(test_data_x - u_sim, axis=1) \
                                        / np.linalg.norm(test_data_x, axis=1)*100).max()
        
        # residual norm
        res_norm[i,j] = err_indicator(u_sim, params, err_type=params['err_type'])
        count += 1
end_time = time()
time_rom = timer_rom[1:].sum()/num_case # from Step 2 to 4
time_sim = 38.4 # seconds
speed_up = time_sim / time_rom
print(f'Time taken: {end_time-start_time:.2f} s, {(end_time-start_time)/60:.2f} mins')
print(f'Average speed up: {speed_up:.2f}')


# In[ ]:


a_grid, w_grid = np.meshgrid(amp_test, width_test)
param_list = np.hstack([a_grid.flatten().reshape(-1,1), w_grid.flatten().reshape(-1,1)])
a_grid, w_grid = np.meshgrid(np.arange(amp_test.size), np.arange(width_test.size))
idx_list = np.hstack([a_grid.flatten().reshape(-1,1), w_grid.flatten().reshape(-1,1)])

idx_param = []
for i,ip in enumerate(params['param']):
    idx = np.argmin(np.linalg.norm(param_list-ip, axis=1))
    idx_param.append((idx, np.array([param_list[idx,0], param_list[idx,1]])))


# In[ ]:


max_err_train = []
res_norm_train = []
for i in idx_param:
    idd = i[0]
    max_err_train.append(max_err[idx_list[idd,0],idx_list[idd,1]])
    res_norm_train.append(res_norm[idx_list[idd,0],idx_list[idd,1]])
max_err_train = np.stack(max_err_train)
res_norm_train = np.stack(res_norm_train)
err_ratio_train = res_norm_train / max_err_train
err_ratio_train_mean = err_ratio_train.mean()
res_norm_tol = err_ratio_train_mean * params['tol2']
print(f'tolerance of residual norm (mean): {res_norm_tol:.5f}')
print(f"tolerance of residual norm (reg_max): {params['tol']:.5f}")


# In[ ]:


fig = plt.figure(figsize=(18,4))
ax = fig.add_subplot(131)
ax.plot(np.arange(1,res_norm.size+1), res_norm.flatten(), 'r.', label='res_norm_'+res_name)
ax.set_ylabel('Errors', fontsize=14)
ax.tick_params(labelsize=14)
ax.grid()
ax.legend(fontsize=14)

ax = fig.add_subplot(132)
ax.plot(np.arange(1,res_norm.size+1), max_err.flatten(), 'b.', label='maxRelErr')
ax.set_ylabel('Errors', fontsize=14)
ax.tick_params(labelsize=14)
ax.grid()
ax.legend(fontsize=14)

ax = fig.add_subplot(133)
x = max_err.flatten().reshape(-1,1)
y = res_norm.flatten().reshape(-1,1)
reg = LinearRegression().fit(x, y)
y_pred = reg.predict(x)
y_diff = y - y_pred
x_test = np.linspace(0,max_err.max(),2).reshape(-1,1)
y_test = reg.predict(x_test)
y_test1 = reg.coef_ * x_test + reg.intercept_+ y_diff.min()
y_test2 = reg.coef_ * x_test + reg.intercept_+ y_diff.max()
y_test3 = reg.coef_ * params['tol2'] + reg.intercept_+ y_diff.max()
ax.plot(max_err.flatten(), res_norm.flatten(), 'g.')
ax.plot(x_test,y_test,'r-')
ax.plot(x_test,y_test1,'r--')
ax.plot(x_test,y_test2,'r--')

x_test2 = np.linspace(0,params['tol2'],2).reshape(-1,1)
y_test2 = err_ratio_train_mean * x_test2
ax.plot([x_test2.max(),x_test2.max()],[0,max(y_test3,res_norm_tol)],'k--')
ax.plot(x_test2.max(),res_norm_tol,'b.',markersize=15, label='mean')
ax.plot(params['tol2'],y_test3,'r.',markersize=15, label='reg')
plt.legend(fontsize=16)

ax.set_xlabel('$Max relative error (%)$', fontsize=16)
ax.set_ylabel('Residual norm', fontsize=16)
ax.tick_params(labelsize=16)
ax.grid()

plt.tight_layout()
plt.savefig(data_path + f'residual_norm_{res_name}_knn{knn}.png')
print(max_err.max())


# In[ ]:


print(f'number of DIs: {len(idx_param)}')
max_err_heatmap(max_err, sindy_idx, params, amp_test, width_test, data_path, idx_list, idx_param,
                xlabel='Width', ylabel='Amplitude', dtype='float')


# In[ ]:


# the residual-based errors are scaled by the "scale" parameter
max_err_heatmap(res_norm, sindy_idx, params, amp_test, width_test, data_path, idx_list, idx_param,
                xlabel='Width', ylabel='Amplitude', label='Residual Norm', dtype='float', scale=100)

