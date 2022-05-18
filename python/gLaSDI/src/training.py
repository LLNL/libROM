import numpy as np
import tensorflow.compat.v1 as tf
tf.compat.v1.disable_v2_behavior()
import pickle
from autoencoder import full_network, define_loss, silu
from sindy_utils import sindy_simulate, derivative
from time import time
import seaborn as sns
import random
from error_utils import *
import copy
from sklearn.linear_model import LinearRegression
import seaborn as sns
import matplotlib.patches as patches
from copy import deepcopy
import matplotlib.pyplot as plt

def train_network(training_data, params):
    latent_dim = params['latent_dim']
    library_dim = params['library_dim']
    test_data = params['test_data']
    if isinstance(test_data, str):
        print('Load test data')
        test_data = pickle.load(open(test_data, "rb"))
    
    if 'save_frequency' not in params.keys():
        params['save_frequency'] = 1e8
    
    if params['retrain']:
        testing_losses = params['testing_losses']
        training_losses = params['training_losses']
        if 'err_array' in params.keys():
            err_array = params['err_array']
        else:
            err_array = []
            
        if 'max_err_idx_param' in params.keys():
            max_err_idx_param = params['max_err_idx_param']
        else:
            max_err_idx_param = [[]]
            
        if 'sindy_idx' in params.keys():
            sindy_idx = params['sindy_idx']
        else:
            sindy_idx = []
            
        if 'epoch_count' in params.keys():
            epoch_count = params['epoch_count']   # epoch counter
            params['max_epochs'] += epoch_count
        else:
            epoch_count = 0
        param_flag = False
        
    else:
        testing_losses = []
        training_losses = []
        err_array = []
        max_err_idx_param = [[]]
        sindy_idx = []
        epoch_count = 0   # epoch counter
    sindy_model_terms = [np.sum(params['coefficient_mask'])]
    
    train_flag = True # flag to indicate whether to perform training
    w = 1 # counter for the level of random subset evalutions
    while epoch_count < params['max_epochs'] and train_flag:
        tf.reset_default_graph() # clear graph
        tf.set_random_seed(params['seed'])
        np.random.seed(params['seed'])
        random.seed(params['seed'])
        
        if epoch_count > 0:
            if param_flag:
                print(f'* Updating graph')
            else:
                print(f'* NOT updating graph and continue Training')
        
        # Create models
        autoencoder_network = full_network(params)
        loss,losses,_ = define_loss(autoencoder_network, params)
        learning_rate = tf.placeholder(tf.float32, name='learning_rate')
        train_op = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(loss)
        saver = tf.train.Saver(var_list=tf.get_collection(tf.GraphKeys.GLOBAL_VARIABLES)) 
        
        param_flag = True # flag to indicate whether to update parameters
        with tf.Session(config = params['config']) as sess:
            sess.run(tf.global_variables_initializer())
            if epoch_count == 0:
                print(f'* Evaluating')
                err_array_tmp, max_err, idx, param_tmp, sindy_idx_tmp = err_map_subset(sess, autoencoder_network, 
                                                                                    params, test_data=test_data,
                                                                                    err_type=params['err_type'])
                testing_losses.append(max_err)
                err_array.append(err_array_tmp)
                max_err_idx_param.append([])
                sindy_idx.append(sindy_idx_tmp)
                
            for i in range(params['update_epoch']):
                for j in range(params['epoch_size']//params['batch_size']):
                    batch_idxs = np.arange(j*params['batch_size'], (j+1)*params['batch_size'])
                    train_dict = create_feed_dictionary(training_data, params, idxs=batch_idxs)
                    sess.run(train_op, feed_dict=train_dict)
                
                if params['print_progress'] and (i % params['print_frequency'] == 0):
                    train_loss = print_progress(sess, epoch_count, loss, losses, train_dict)
                    training_losses.append(train_loss)
                epoch_count += 1
                
                if params['sequential_thresholding'] and (i % params['threshold_frequency'] == 0) and (i > 0):
                    params['coefficient_mask'] = np.abs(sess.run(autoencoder_network['sindy_coefficients'])) > params['coefficient_threshold']
                    print('THRESHOLDING: %d active coefficients' % np.sum(params['coefficient_mask']))
                    sindy_model_terms.append(np.sum(params['coefficient_mask']))
                    
                if i % params['save_frequency'] == 0: # Save current model state 
                    params['training_losses'] = training_losses
                    params['epoch_count'] = epoch_count
                    sindy_coefficients = []
                    for i in range(len(autoencoder_network['sindy_coefficients'])):
                        sindy_coefficients.append(sess.run(autoencoder_network['sindy_coefficients'][i], feed_dict={}))
                    params['model_params'] = []
                    params['model_params'].append(sindy_coefficients)
                    params['model_params'].append(sess.run(autoencoder_network['encoder_weights'], feed_dict={}))
                    params['model_params'].append(sess.run(autoencoder_network['encoder_biases'], feed_dict={}))
                    params['model_params'].append(sess.run(autoencoder_network['decoder_weights'], feed_dict={}))
                    params['model_params'].append(sess.run(autoencoder_network['decoder_biases'], feed_dict={}))
                    params['coeff_exist'] = True
                    pickle.dump(params, open(params['fig_path'] + params['save_name'] + '_params.pkl', 'wb'))
                    saver.save(sess, params['fig_path'] + params['save_name'])
            
            
            # Evaluate current model on a random subset of parameters using a specified error indicator
            print(f'* Evaluating')
            err_array_tmp, max_err, idx, param_tmp, sindy_idx_tmp = err_map_subset(sess, autoencoder_network, 
                                                                                params, test_data=test_data,
                                                                                err_type=params['err_type'])
            testing_losses.append(max_err)
            err_array.append(err_array_tmp)
            max_err_idx_param.append(copy.copy(max_err_idx_param[-1]))
            sindy_idx.append(sindy_idx_tmp)
            
            
            # Update tolerance for error indicator
            tol_old = params['tol']
            params['tol'], max_rel_err = update_tol(sess, autoencoder_network, params, 
                                       training_data, err_type=params['err_type'])
            print(f"  Max rel. err.: {max_rel_err:.1f}%, Update tolerance for error indicator from {tol_old:.5f} to {params['tol']:.5f}")
            
            
            # Save current model state 
            # should save the state before new parameters are added, consistent with the graph
            params['training_losses'] = training_losses
            params['testing_losses'] = testing_losses
            params['err_array'] = err_array
            params['max_err_idx_param'] = max_err_idx_param
            params['sindy_idx'] = sindy_idx
            params['epoch_count'] = epoch_count
            
            sindy_coefficients = []
            for i in range(len(autoencoder_network['sindy_coefficients'])):
                sindy_coefficients.append(sess.run(autoencoder_network['sindy_coefficients'][i], feed_dict={}))
            params['model_params'] = []
            params['model_params'].append(sindy_coefficients)
            params['model_params'].append(sess.run(autoencoder_network['encoder_weights'], feed_dict={}))
            params['model_params'].append(sess.run(autoencoder_network['encoder_biases'], feed_dict={}))
            params['model_params'].append(sess.run(autoencoder_network['decoder_weights'], feed_dict={}))
            params['model_params'].append(sess.run(autoencoder_network['decoder_biases'], feed_dict={}))
            params['coeff_exist'] = True
            
            pickle.dump(params, open(params['fig_path'] + params['save_name'] + '_params.pkl', 'wb'))
            final_losses = sess.run((losses['decoder'], losses['sindy_x'], losses['sindy_z'],
                                     losses['sindy_regularization']), feed_dict=train_dict)
            saver.save(sess, params['fig_path'] + params['save_name'])
            
            
            # Update training data
            for i in training_data['param']:
                if np.linalg.norm(i-test_data['param'][idx]) < 1e-8:
                    print(f"  PARAMETERS EXIST, NOT adding it!")
                    param_flag = False
                    break
            if param_flag:
                print(f'* Update Training set: add case {param_tmp}')
                training_data['data'].append(test_data['data'][idx])
                training_data['param'].append(test_data['param'][idx])
                params['num_sindy'] += 1
                params['param'] = training_data['param']
                params['train_idx'].append(idx)
                max_err_idx_param[-1].append((idx, param_tmp))
                    
            # Update random subset size
            subset_ratio = params['subsize']/params['num_test']*100 # new subset size
            if max_err <= params['tol']:
                w += 1
                if params['subsize']*2 <= params['num_test']:
                    params['subsize'] *= 2 # double the random subset size for evaluation
                else:
                    params['subsize'] = params['num_test']
                subset_ratio = params['subsize']/params['num_test']*100 # new subset size
                print(f"  Max error indicator <= Tol! Current subset ratio {subset_ratio:.1f}%")

            if 'sindy_max' in params.keys() and params['sindy_max'] != None:
                if params['num_sindy'] == params['sindy_max']+1:
                    print(f"  Max # SINDys {params['sindy_max']:d} is reached! Training done!")
                    train_flag = False 
            elif subset_ratio >= params['subsize_max']: 
                print(f"  Current subset ratio {subset_ratio:.1f}% >= Target subset ratio {params['subsize_max']:.1f}%!")
                train_flag = False 
                
    results_dict = {}
    results_dict['num_epochs'] = epoch_count
    results_dict['sindy_coefficients'] = sindy_coefficients
    results_dict['model_params'] = params['model_params']
    results_dict['loss_decoder'] = final_losses[0]
    results_dict['loss_decoder_sindy'] = final_losses[1]
    results_dict['loss_sindy'] = final_losses[2]
    results_dict['loss_sindy_regularization'] = final_losses[3]
    results_dict['testing_losses'] = np.array(testing_losses)
    results_dict['training_losses'] = np.array(training_losses)
    results_dict['sindy_model_terms'] = np.array(sindy_model_terms)
    results_dict['err_array'] = err_array
    results_dict['max_err_idx_param'] = max_err_idx_param
    results_dict['sindy_idx'] = sindy_idx
    return results_dict
    

def NN(x, weights, biases, activation):
    """
    This networks serve as either an encoder or a decoder.
    """
    num_layers = len(weights)
    for i in range(num_layers-1):
        x = np.matmul(x, weights[i]) + biases[i]
        if activation == 'tanh':
            x = np.tanh(x)
        elif activation == 'sigmoid':
            x = 1.0 / (1.0 + np.exp(-x)) # sigmoid
            
    # output layer (linear activation)
    x = np.matmul(x, weights[-1]) + biases[-1]
    return x


def eval_perf(sess, tensorflow_run_tuple, autoencoder_network, params, test_data, test_param, idx=None):
    """
    This function evaluates the model on a given testing parameter case.
    
    inputs:
        sess: the tensorflow session
        autoencoder_network: the autoencoder network
        test_data: data of testing parameter case (could be updated to just provide the 
                    initial condition of the testing parameter case)
        idx: the index of the DI used for evaluation; used when knn=1; if knn=1 and it is None, 
                the DI closest to the testing parameter will be used.
    
    outputs:
        u_sim: prediction of full-order model solutions by DI and decoder
        idx: the index of the DI closest to the testing parameter based on the Euclidean distance
    """
    if 'include_sine' in params.keys():
        include_sine = params['include_sine']
    else:
        include_sine = False
        
    if 'include_cosine' in params.keys():
        include_cosine = params['include_cosine']
    else:
        include_cosine = False
        
    test_dictionary = create_feed_dictionary2(test_data, params, idxs=1)
    tf_results = sess.run(tensorflow_run_tuple, feed_dict=test_dictionary)
    
    test_set_results = {}
    for i,key in enumerate(autoencoder_network.keys()):
        test_set_results[key] = tf_results[i]

    if idx == None:
        train_param = np.stack(params['param'])
        idx = np.argmin(np.linalg.norm(train_param-test_param, axis=1))

    # calculate SINDy coefficients
    if params['convex_knn'] == 1: # nearest SINDy's coefficient
        sindy_coeff = test_set_results['sindy_coefficients'][idx]   
        
    else: # KNN convex interpolation of coefficients
        dist = np.linalg.norm(train_param-test_param, axis=1)
        knn_idx = np.argsort(dist)[:params['convex_knn']]
        phi = np.zeros_like(knn_idx)
        if dist[knn_idx[0]] == 0: # check if the min distance is zero
            phi[0] = 1
        else:
            phi = 1 / np.linalg.norm(train_param[knn_idx]-test_param, axis=1)**2
        psi = phi / phi.sum()

        sindy_coeff = np.zeros(test_set_results['sindy_coefficients'][0].shape)
        for i,kidx in enumerate(knn_idx):
            sindy_coeff += psi[i] * test_set_results['sindy_coefficients'][kidx]

    # predict latent-space dynamics given the initial condition of latent-space variables
    z_sim = sindy_simulate(test_set_results['z'][0,idx,:], 
                           test_data['t'].squeeze(), 
                           sindy_coeff, params['poly_order'], 
                           include_sine, include_cosine)
    u_sim = NN(z_sim, test_set_results['decoder_weights'], 
               test_set_results['decoder_biases'], params['activation'])
    return u_sim, idx


def err_map_subset(sess, autoencoder_network, params, test_data=None, err_type=1):
    """
    This function computes errors in random subsets of the parameter space 
    using a speciffied error indicator. The subset size and the threshold for 
    error indicator are adjusted during training.
    inputs:
        test_data: dict, testing data
        err_type: int, types of error indicator. 
                1: max relative error (if test data available)
                2: residual norm (mean)
    outputs:
        err_array: ndarray, errors in the parameter space
        err_max: float, max error
        err_idx: int, index of the max error
        test_data['param'][err_idx]: ndarray, parameters of the case with max error
        sindy_idx: ndarray, indices of local SINDys used for evaluation
    """
    tensorflow_run_tuple = ()
    for key in autoencoder_network.keys():
        tensorflow_run_tuple += (autoencoder_network[key],)
    
    amp = params['test_param'][:,0]
    width = params['test_param'][:,1]
    err_array = np.zeros([amp.size, width.size])
    err_array2 = np.zeros([amp.size, width.size])
    sindy_idx = np.zeros([amp.size, width.size])
    
    # select a random subset for evaluation
    rng = np.random.default_rng()
    a = np.setdiff1d(np.arange(params['num_test']),params['train_idx']) # exclude existing training cases
    rng.shuffle(a)
    subset = a[:params['subsize']]
    
    count = 0
    start_time = time()
    for i,a in enumerate(amp):
        for j,w in enumerate(width):
            if count in subset:
                u_sim, idx = eval_perf(sess, tensorflow_run_tuple, autoencoder_network, params,
                                       test_data['data'][count], test_data['param'][count])
                sindy_idx[i,j] = idx+1
                params['pde']['param'] = [a, w]
                err_array[i,j] = err_indicator(u_sim, params, 
                                               data=test_data['data'][count]['x'], 
                                               err_type=err_type) # residual norm of all time steps
            else:
                sindy_idx[i,j] = -1
                err_array[i,j] = -1
            count += 1
    end_time = time()
    err_max = err_array.max()
    err_idx = np.argmax(err_array)
    print(f"  Time: {end_time-start_time:.2f} s, Case: {test_data['param'][err_idx]}, Tol: {params['tol']:.5f}, Max Error: {err_max:.6f}")
    return err_array, err_max, err_idx, test_data['param'][err_idx], sindy_idx


def update_tol(sess, autoencoder_network, params, training_data, err_type=1):
    """
    This function computes the error indicator and max relative errors of existing 
    training cases in order to update the tolerance threshold of the error indicator.
    inputs:
        training_data: dict, training data
        err_type: int, types of error indicator. 
                1: max relative error (if test data available)
                2: residual norm (mean), 1D Burgers' eqn
                3: residual norm (mean), 2D Burgers' eqn
    outputs:
        tol_new: float, updated tolerance for the error indicator
    """
    tensorflow_run_tuple = ()
    for key in autoencoder_network.keys():
        tensorflow_run_tuple += (autoencoder_network[key],)
    
    num_sindy = params['num_sindy']
    err1 = np.zeros(num_sindy)
    err2 = np.zeros(num_sindy)
    for i in range(num_sindy):
        u_sim,_ = eval_perf(sess, tensorflow_run_tuple, autoencoder_network, params,
                            training_data['data'][i], training_data['param'][i])
        params['pde']['param'] = [training_data['param'][i][0], training_data['param'][i][1]]
        err1[i] = err_indicator(u_sim, params, data=training_data['data'][i]['x'], err_type=err_type) # residual norm
        err2[i] = err_indicator(u_sim, params, data=training_data['data'][i]['x'], err_type=1) # max relative error

    # update tolerance of error indicator
    if params['adaptive'] == 'mean':
        tol_new = (err1 / err2).mean() * params['tol2']
    elif params['adaptive'] == 'last':
        tol_new = (err1[-1] / err2[-1]).mean() * params['tol2']
    else:
        x = err2.reshape(-1,1)
        y = err1.reshape(-1,1)
        reg = LinearRegression().fit(x, y)
        if params['adaptive'] == 'reg_mean':
            tol_new = max(0, reg.coef_[0][0] * params['tol2'] + reg.intercept_[0])
            print(reg.coef_[0][0], reg.intercept_[0])
            
        elif params['adaptive'] == 'reg_max':
            y_diff = y - reg.predict(x)
            tol_new = max(0, reg.coef_[0][0] * params['tol2'] + reg.intercept_[0] + y_diff.max())
            
        elif params['adaptive'] == 'reg_min':
            y_diff = y - reg.predict(x)
            tol_new = max(0, reg.coef_[0][0] * params['tol2'] + reg.intercept_[0] + y_diff.min())

    return tol_new, err2.max()



def print_progress(sess, i, loss, losses, train_dict):
    """
    Print loss function values to keep track of the training progress

    inputs:
        sess: the tensorflow session
        i: the training iteration
        loss: tensorflow object representing the total loss function used in training
        losses: tuple of the individual losses that make up the total loss
        train_dict: feed dictionary of training data

    outputs:
        Tuple of losses calculated on the training set
    """
    training_loss_vals = sess.run((loss,) + tuple(losses.values()), feed_dict=train_dict)
    print("Epoch %d" % i)
    print("  train loss: {0:.4e}, decoder: {1:.4e}, sindy-x: {2:.4e}, sindy-z: {3:.4e}, sindy-reg: {4:.4f}".format(training_loss_vals[0][0], training_loss_vals[1][0], training_loss_vals[2][0], training_loss_vals[3][0], training_loss_vals[4][0]))
    return training_loss_vals



def create_feed_dictionary(data, params, idxs=None):
    """
    Create the feed dictionary for passing into tensorflow

    inputs:
        data: Dictionary object containing the data to be passed in. Must contain input data x,
        along the first (and possibly second) order time derivatives dx (ddx).
        params: Dictionary object containing model and training parameters. The relevant
        parameters are sequential_thresholding (which indicates whether or not
        coefficient thresholding is performed), coefficient_mask (optional if sequential
        thresholding is performed; 0/1 mask that selects the relevant coefficients in the SINDy
        model), and learning rate (float that determines the learning rate).
        idxs: Optional array of indices that selects which examples from the dataset are passed
        into tensorflow. If None, all examples are used.

    outputs:
        feed_dict: Dictionary object containing the relevant data to pass to tensorflow
    """
    num_sindy = params['num_sindy']
    if idxs is None:
        idxs = np.arange(data['data'][0]['x'].shape[0])
    feed_dict = {}
    
    data_x = []
    data_dx = []
    for i in range(num_sindy):
        data_x.append(data['data'][i]['x'][idxs])
        data_dx.append(data['data'][i]['dx'][idxs])
        
    feed_dict['x:0'] = np.stack(data_x, axis=1)   # [batch,num_sindy,input_dim]
    feed_dict['dx:0'] = np.stack(data_dx, axis=1) # [batch,num_sindy,input_dim]
    
    if params['sequential_thresholding']:
        feed_dict['coefficient_mask:0'] = params['coefficient_mask']
        
    feed_dict['learning_rate:0'] = params['learning_rate']
    return feed_dict


def create_feed_dictionary2(data, params, idxs=None):
    """
    Create the feed dictionary for passing into tensorflow
    For testing set that has only 1 case!

    inputs:
        data: Dictionary object containing the data to be passed in. Must contain input data x,
        along the first (and possibly second) order time derivatives dx (ddx).
        params: Dictionary object containing model and training parameters. The relevant
        parameters are sequential_thresholding (which indicates whether or not
        coefficient thresholding is performed), coefficient_mask (optional if sequential
        thresholding is performed; 0/1 mask that selects the relevant coefficients in the SINDy
        model), and learning rate (float that determines the learning rate).
        idxs: Optional array of indices that selects which examples from the dataset are passed
        into tensorflow. If None, all examples are used.

    outputs:
        feed_dict: Dictionary object containing the relevant data to pass to tensorflow
    """
    num_sindy = params['num_sindy']
    feed_dict = {}

    data_x = []    
    data_dx = []
    if 'data' in data.keys():
        if idxs == None:
            idxs = data['data'][0]['x'].shape[0]
        for i in range(num_sindy):
            data_x.append(data['data'][0]['x'][:idxs])   # initial condition
            data_dx.append(data['data'][0]['dx'][:idxs]) # initial condition

    else:
        if idxs == None:
            idxs = data['x'].shape[0]
        for i in range(num_sindy):
            data_x.append(data['x'][:idxs])   # initial condition
            data_dx.append(data['dx'][:idxs]) # initial condition
            
    feed_dict['x:0'] = np.stack(data_x, axis=1)   # [batch,num_sindy,input_dim]
    feed_dict['dx:0'] = np.stack(data_dx, axis=1) # [batch,num_sindy,input_dim]
    
    if params['sequential_thresholding']:
        feed_dict['coefficient_mask:0'] = params['coefficient_mask']
        
    feed_dict['learning_rate:0'] = params['learning_rate']
    return feed_dict


def eval_model(test_data, params, test_param, idx=None, knn=4, calc_dz=False, calc_du=False):
    """
    This function evaluates the trained gLaSDI model.
    """
    timer = []
    
    # Step 1: set up tf graph and load parameters, can be optimized, excluded from ROM computational time
    timer.append(time()) 
    include_sine = False
    include_cosine = False
    if 'include_sine' in params.keys():
        include_sine = params['include_sine']
    if 'include_cosine' in params.keys():
        include_cosine = params['include_cosine']
        
    z_encoder = NN(test_data['x'], params['model_params'][1], params['model_params'][2], params['activation']) # encoder
    u_decoder = NN(z_encoder, params['model_params'][3], params['model_params'][4], params['activation']) # decoder
    
    
    # Step 2: find the nearest neighbor (optional)
    timer.append(time()) 
    if idx == None:
        train_param = np.stack(params['param'])
        idx = np.argmin(np.linalg.norm(train_param-test_param, axis=1))
    
    
    # Step 3: calculate SINDy coefficients
    timer.append(time())        
    if knn == 1:
        print(f"Index of the nearest local SINDy: {idx+1}")
        sindy_coeff = params['model_params'][0][idx]
        
    else: # KNN convex interpolation of coefficients
        dist = np.linalg.norm(train_param-test_param, axis=1)
        knn_idx = np.argsort(dist)[:knn]
        phi = np.zeros_like(knn_idx)
        if dist[knn_idx[0]] == 0: # check if the min distance is zero
            phi[0] = 1
        else:
            phi = 1 / np.linalg.norm(train_param[knn_idx]-test_param, axis=1)**2
        psi = phi / phi.sum()

        sindy_coeff = np.zeros(params['model_params'][0][0].shape)
        for i,kidx in enumerate(knn_idx):
            sindy_coeff += psi[i] * params['model_params'][0][kidx]
            

    # Step 4: lastent-space dynamics prediction and obtain physical dynamics
    timer.append(time())
    z_sim = sindy_simulate(z_encoder[0,:], test_data['t'].squeeze(), 
                           sindy_coeff, params['poly_order'], 
                           include_sine,include_cosine)
    u_sim = NN(z_sim, params['model_params'][3], params['model_params'][4], params['activation'])

    timer.append(time())
    
    if calc_dz:
        dz_encoder = derivative(z_encoder,params['pde']['tstop'])
        dz_sim = derivative(z_sim,params['pde']['tstop'])
    else:
        dz_encoder = 0
        dz_sim = 0
        
    if calc_du:
        du_decoder = derivative(u_decoder,params['pde']['tstop'])
        du_sim = derivative(u_sim,params['pde']['tstop'])
    else:
        du_decoder = 0
        du_sim = 0
    
    timer1 = np.array(timer)
    timer2 = timer1[1:]
    timer_rom = timer2 - timer1[:-1]
    return u_decoder, du_decoder, u_sim, du_sim, z_encoder, dz_encoder, z_sim, dz_sim, idx, timer_rom


# heat map of max relative errors
def max_err_heatmap(max_err, sindy_idx, params, p1_test, p2_test, data_path, idx_list=[], idx_param=[], 
                    xlabel='param1', ylabel='param2', label='Max. Relative Error (%)', dtype='int', scale=1):
    sns.set(font_scale=1.3)
    if dtype == 'int':
        max_err = max_err.astype(int)
        fmt1 = 'd'
    else:
        fmt1 = '.1f'
    rect = []
    for i in range(len(idx_param)):
        print(f"idx: {idx_param[i][0]}, param: {idx_param[i][1]}")
        idd = idx_param[i][0]
        rect.append(patches.Rectangle((idx_list[idd,1], idx_list[idd,0]), 1, 1, linewidth=2, edgecolor='k', facecolor='none'))
    rect2 = deepcopy(rect)
    
    if max_err.size < 100:
        fig = plt.figure(figsize=(5,5))
    else:
        fig = plt.figure(figsize=(9,9))

    fontsize = 14
    if max_err.max() >= 10:
        fontsize = 12
        max_err = max_err.astype(int)
        fmt1 = 'd'
    ax = fig.add_subplot(111)
    cbar_ax = fig.add_axes([0.99, 0.19, 0.02, 0.7])
    if label == 'Residual Norm':
        vmax = params['tol']*scale
    else:
        vmax = max_err.max()*scale
    sns.heatmap(max_err*scale, ax=ax, square=True, 
                xticklabels=p2_test, yticklabels=p1_test, 
                annot=True, annot_kws={'size':fontsize}, fmt=fmt1, 
                cbar_ax=cbar_ax, cbar=True, cmap='vlag', robust=True, vmin=0, vmax=vmax)
        
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
    ax.set_xlabel(xlabel, fontsize=24)
    ax.set_ylabel(ylabel, fontsize=24)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30) 
    
    plt.tight_layout()
    if label == 'Residual Norm':
        plt.savefig(data_path + f'heatmap_resNorm.png', bbox_inches='tight')
    else:
        plt.savefig(data_path + f'heatmap_maxRelErr_glasdi.png', bbox_inches='tight')
    plt.show()