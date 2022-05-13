import numpy as np
import tensorflow.compat.v1 as tf
tf.compat.v1.disable_v2_behavior()
import pickle
from autoencoder_LaSDI import AE_network, DI_network, AE_loss, DI_loss
from sindy_utils import sindy_simulate
import pickle
from time import time
import seaborn as sns
import random
from error_utils import *
import copy


def train_autoencoder(training_data, params):
    tf.reset_default_graph() # clear graph
    tf.set_random_seed(params['seed'])
    np.random.seed(params['seed'])
    random.seed(params['seed'])
    num_batch = params['epoch_size_AE']//params['batch_size']+1
    
    # Create models
    autoencoder = AE_network(params)
    loss = AE_loss(autoencoder)
    learning_rate = tf.placeholder(tf.float32, name='learning_rate')
    train_op = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(loss)
    saver = tf.train.Saver(var_list=tf.get_collection(tf.GraphKeys.GLOBAL_VARIABLES)) 
    
    # Train models
    params['train_losses_AE'] = []
    with tf.Session(config = params['config']) as sess:
        sess.run(tf.global_variables_initializer())

        for i in range(params['max_epochs_AE']):
            for j in range(num_batch):
                if j == num_batch-1:
                    batch_idxs = np.arange(j*params['batch_size'], params['epoch_size_AE'])
                else:
                    batch_idxs = np.arange(j*params['batch_size'], (j+1)*params['batch_size'])
                train_dict = create_feed_dictionary(training_data, params, idxs=batch_idxs)
                sess.run(train_op, feed_dict=train_dict)

            train_dict = create_feed_dictionary(training_data, params)
            train_loss = sess.run(loss, feed_dict=train_dict)
            params['train_losses_AE'].append(train_loss)
            
            if i % params['print_frequency'] == 0:
                print(f'Epoch {i}, loss: {train_loss:.4e}')
                
        tf_run_tuple = ()
        for key in autoencoder.keys():
            tf_run_tuple += (autoencoder[key],)      
        train_dict = create_feed_dictionary(training_data, params)
        tf_results = sess.run(tf_run_tuple, feed_dict=train_dict)
        AE_results = {}
        for i,key in enumerate(autoencoder.keys()):
            AE_results[key] = tf_results[i]
            
        params['AE_params'] = []
        params['AE_params'].append(AE_results['encoder_weights'])
        params['AE_params'].append(AE_results['encoder_biases'])
        params['AE_params'].append(AE_results['decoder_weights'])
        params['AE_params'].append(AE_results['decoder_biases'])
        pickle.dump(params, open(params['fig_path'] + params['save_name'] + '_params.pkl', 'wb'))
        saver.save(sess, params['fig_path'] + params['save_name'])
        
    return params, AE_results['x_decode'], AE_results['z']


def train_DI(training_data, params):
    tf.reset_default_graph() # clear graph
    tf.set_random_seed(params['seed'])
    np.random.seed(params['seed'])
    random.seed(params['seed'])
    num_batch = params['epoch_size_DI']//params['batch_size']+1
    
    # Create models
    DI = DI_network(params)
    loss = DI_loss(DI, params)
    learning_rate = tf.placeholder(tf.float32, name='learning_rate')
    train_op = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(loss)
    saver = tf.train.Saver(var_list=tf.get_collection(tf.GraphKeys.GLOBAL_VARIABLES)) 

    # Train models
    params['train_losses_DI'] = []
    with tf.Session(config = params['config']) as sess:
        sess.run(tf.global_variables_initializer())

        for i in range(params['max_epochs_DI']):
            for j in range(num_batch):
                if j == num_batch-1:
                    batch_idxs = np.arange(j*params['batch_size'], params['epoch_size_DI'])
                else:
                    batch_idxs = np.arange(j*params['batch_size'], (j+1)*params['batch_size'])
                train_dict = create_feed_dictionary2(training_data, params, idxs=batch_idxs)
                sess.run(train_op, feed_dict=train_dict)

            train_dict = create_feed_dictionary2(training_data, params)
            train_loss = sess.run(loss, feed_dict=train_dict)
            params['train_losses_DI'].append(train_loss)
            
            if i % params['print_frequency'] == 0:
                print(f'Epoch {i}, loss: {train_loss:.4e}')

        params['DI_params'] = []
        for i in range(params['num_DI']):
            params['DI_params'].append(sess.run(DI['sindy_coefficients'][i], feed_dict={}))
            
        pickle.dump(params, open(params['fig_path'] + params['save_name'] + '_params.pkl', 'wb'))
        saver.save(sess, params['fig_path'] + params['save_name'])
                    
    return params


def create_feed_dictionary(data, params, idxs=None):
    if idxs is None:
        idxs = np.arange(data.shape[0])
    feed_dict = {}
    
    feed_dict['x:0'] = data[idxs]
    feed_dict['learning_rate:0'] = params['learning_rate']
    return feed_dict


def create_feed_dictionary2(data, params, idxs=None):
    if idxs is None:
        idxs = np.arange(data['data'][0]['z'].shape[0])
    feed_dict = {}
    
    data_z = []
    data_dz = []
    for i in range(len(data['data'])):
        data_z.append(data['data'][i]['z'][idxs])
        data_dz.append(data['data'][i]['dz'][idxs])

    feed_dict['z:0'] = np.stack(data_z, axis=1)   # [batch,num_DI,input_dim]
    feed_dict['dz:0'] = np.stack(data_dz, axis=1) # [batch,num_DI,input_dim]
        
    if params['sequential_thresholding']:
        feed_dict['coefficient_mask:0'] = params['coefficient_mask']
        
    feed_dict['learning_rate:0'] = params['learning_rate']
    return feed_dict
