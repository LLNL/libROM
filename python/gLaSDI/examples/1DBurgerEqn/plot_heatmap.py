#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import pickle 
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
import seaborn as sns
from copy import deepcopy
from matplotlib.ticker import FormatStrFormatter


# In[ ]:


data_path = "./fig/nCase441_k1_MRN_ld5_subsize50_upfreq2000/"
save_name = 'burger_2021_11_26_21_50_45'
params = pickle.load(open(data_path + save_name + '_params.pkl', 'rb'))


# In[ ]:


# heat map of max relative errors
sns.set(font_scale=1.3)
def max_err_heatmap(max_err, sindy_idx, idx_list, idx_param, dtype='int', scale=1, mask=None):
    if dtype == 'int':
        max_err = max_err.astype(int)
    rect = []
    rect.append(patches.Rectangle((0, 0), 1, 1, linewidth=2, edgecolor='k', facecolor='none'))
    rect.append(patches.Rectangle((amp_test.size-1, 0), 1, 1, linewidth=2, edgecolor='k', facecolor='none'))
    rect.append(patches.Rectangle((amp_test.size-1, width_test.size-1), 1, 1, linewidth=2, edgecolor='k', facecolor='none'))
    rect.append(patches.Rectangle((0, width_test.size-1), 1, 1, linewidth=2, edgecolor='k', facecolor='none'))
    for i in range(len(idx_param)):
        print(f"idx: {idx_param[i][0]}, param: {idx_param[i][1]}")
        idd = idx_param[i][0]
        rect.append(patches.Rectangle((idx_list[idd,0], idx_list[idd,1]), 1, 1, 
                                      linewidth=2, edgecolor='k', facecolor='none'))
    idx_max_err = np.argmax(max_err)
    rect.append(patches.Rectangle((idx_list[idx_max_err,0], idx_list[idx_max_err,1]), 1, 1, 
                                      linewidth=3, edgecolor='lime', facecolor='none'))
    rect2 = deepcopy(rect)
    
    if max_err.size < 100:
        fig = plt.figure(figsize=(10,5))
    else:
        fig = plt.figure(figsize=(18,9))
    
    # local DI indices
    ax = fig.add_subplot(121)
    sindy_idx = sindy_idx.astype(int)
    sns.heatmap(sindy_idx, ax=ax, square=True, xticklabels=width_test, yticklabels=amp_test, 
                annot=True, fmt='d', annot_kws={'size':14}, 
                cbar=False, cmap='Spectral', robust=True, vmin=1, vmax=15, mask=mask)
    for i in rect:
        ax.add_patch(i)
        
    # format text labels
    fmt = '{:0.2f}'
    xticklabels = []
    for item in ax.get_xticklabels():
        item.set_text(fmt.format(float(item.get_text())))
        xticklabels += [item]
    yticklabels = []
    for item in ax.get_yticklabels():
        item.set_text(fmt.format(float(item.get_text())))
        yticklabels += [item]
    ax.set_xticklabels(xticklabels)
    ax.set_yticklabels(yticklabels)
    ax.set_xlabel('Width', fontsize=16)
    ax.set_ylabel('Amplitude', fontsize=16)
    ax.set_title('Index of Selected Local DI', fontsize=16)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30) 

    # heatmap of residual-based error indicator
    ax = fig.add_subplot(122)
    cbar_ax = fig.add_axes([0.99, 0.19, 0.018, 0.7])
    sns.heatmap(max_err*scale, ax=ax, square=True, xticklabels=width_test, yticklabels=amp_test, 
                annot=True, annot_kws={'size':14}, fmt='.1f', 
                cbar_ax=cbar_ax, cbar=True, cmap='vlag', robust=True, vmin=0, vmax=max_err.max()*scale, mask=mask)
    for i in rect2:
        ax.add_patch(i)
        
    # format text labels
    fmt = '{:0.2f}'
    xticklabels = []
    for item in ax.get_xticklabels():
        item.set_text(fmt.format(float(item.get_text())))
        xticklabels += [item]
    yticklabels = []
    for item in ax.get_yticklabels():
        item.set_text(fmt.format(float(item.get_text())))
        yticklabels += [item]
    ax.set_xticklabels(xticklabels)
    ax.set_yticklabels(yticklabels)
    ax.set_xlabel('Width', fontsize=16)
    ax.set_ylabel('Amplitude', fontsize=16)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30) 
    
    plt.tight_layout()
    plt.savefig(data_path + f'heatmap_{len(idx_param)}.png', bbox_inches='tight')
    plt.show()


# In[ ]:


# plot heatmaps of residual-based errors at all greedy sampling steps
amp_test = params['test_param'][:,0]
width_test = params['test_param'][:,1]
a_grid, w_grid = np.meshgrid(np.arange(amp_test.size), np.arange(width_test.size))
idx_list = np.hstack([a_grid.flatten().reshape(-1,1), w_grid.flatten().reshape(-1,1)])
scale = 10
for i in range(len(params['err_array'])):
    idx_param = params['max_err_idx_param'][i]
    mask = params['sindy_idx'][i]==-1
    max_err_heatmap(params['err_array'][i], params['sindy_idx'][i], idx_list, 
                    idx_param, dtype='float', scale=scale, mask=mask)

