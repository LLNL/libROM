import tensorflow.compat.v1 as tf
tf.compat.v1.disable_v2_behavior()
from tensorflow.python.ops.parallel_for.gradients import jacobian, batch_jacobian
from time import time
import pickle
import numpy as np


def AE_network(params):
    input_dim = params['input_dim']
    latent_dim = params['latent_dim']
    activation = params['activation']
    model_params = []
    
    x = tf.placeholder(tf.float32, shape=[None, input_dim], name='x')
    z, x_decode, encoder_weights, encoder_biases, decoder_weights, decoder_biases = nonlinear_autoencoder(x, input_dim, latent_dim,
                                                                                                          params['widths'],
                                                                                                          model_params,
                                                                                                          activation=activation)
    network = {}
    network['x'] = x
    network['x_decode'] = x_decode
    network['z'] = z
    network['encoder_weights'] = encoder_weights
    network['encoder_biases'] = encoder_biases
    network['decoder_weights'] = decoder_weights
    network['decoder_biases'] = decoder_biases
    return network
    
    
def DI_network(params):
    latent_dim = params['latent_dim']
    poly_order = params['poly_order']
    num_DI = params['num_DI'] # number of local DIs
    model_params = []
        
    if 'include_sine' in params.keys():
        include_sine = params['include_sine']
    else:
        include_sine = False
        
    if 'include_cosine' in params.keys():
        include_cosine = params['include_cosine']
    else:
        include_cosine = False
        
    library_dim = params['library_dim']
    network = {}

    z = tf.placeholder(tf.float32, shape=[None, num_DI, latent_dim], name='z')
    dz = tf.placeholder(tf.float32, shape=[None, num_DI, latent_dim], name='dz')
    
    # Dynamics identification models
    Theta = [] # matrix of basis functions
    sindy_coefficients = []
    sindy_predict = [] # [batch,num_DI,latent_dim]

    if params['sequential_thresholding']: # 1 mask for all local DIs
        coefficient_mask = tf.placeholder(tf.float32, shape=[library_dim,latent_dim], name='coefficient_mask')
        network['coefficient_mask'] = coefficient_mask
        
    for i in range(num_DI):
        Theta.append(sindy_library_tf(z[:,i,:], latent_dim, poly_order, include_sine, include_cosine))

        if params['coefficient_initialization'] == 'xavier':
            sindy_coefficients.append(tf.get_variable(f'sindy_coefficients{i}', shape=[library_dim,latent_dim], 
                                                      initializer=tf.truncated_normal_initializer()))
        elif params['coefficient_initialization'] == 'specified':
            sindy_coefficients.append(tf.get_variable(f'sindy_coefficients{i}', 
                                                      initializer=params['init_coefficients']))
        elif params['coefficient_initialization'] == 'constant':
            sindy_coefficients.append(tf.get_variable(f'sindy_coefficients{i}', shape=[library_dim,latent_dim], 
                                                      initializer=tf.constant_initializer(1.0)))
        elif params['coefficient_initialization'] == 'normal':
            sindy_coefficients.append(tf.get_variable(f'sindy_coefficients{i}', shape=[library_dim,latent_dim],
                                                      initializer=tf.initializers.random_normal()))
        
        if params['sequential_thresholding']:
            sindy_predict.append(tf.matmul(Theta[i], coefficient_mask * sindy_coefficients[i]))
        else:
            sindy_predict.append(tf.matmul(Theta[i], sindy_coefficients[i])) # [batch,input_dim]
    sindy_predict = tf.stack(sindy_predict, axis=1) # [batch,num_DI,latent_dim]
    
    network['z'] = z   # [batch,num_DI,latent_dim]
    network['dz'] = dz # [batch,num_DI,latent_dim]
    network['dz_predict'] = sindy_predict # [batch,num_DI,latent_dim]
    network['Theta'] = Theta # list
    network['sindy_coefficients'] = sindy_coefficients # list
    
    return network


def AE_loss(network):
    x = network['x']              
    x_decode = network['x_decode']
    return tf.reduce_mean((x - x_decode)**2)


def DI_loss(network, params):
    dz = network['dz']                 # [batch,num_DI,latent_dim]
    dz_predict = network['dz_predict'] # [batch,num_DI,latent_dim]
    loss = 0
    for i in range(params['num_DI']):
        loss += tf.reduce_mean((dz[:,i,:] - dz_predict[:,i,:])**2)
    return loss


def nonlinear_autoencoder(x, input_dim, latent_dim, widths, model_params, activation='elu'):
    """
    Construct a nonlinear autoencoder.

    Arguments:

    Returns:
        z -
        x_decode -
        encoder_weights - List of tensorflow arrays containing the encoder weights
        encoder_biases - List of tensorflow arrays containing the encoder biases
        decoder_weights - List of tensorflow arrays containing the decoder weights
        decoder_biases - List of tensorflow arrays containing the decoder biases
    """
    if activation == 'relu':
        activation_function = tf.nn.relu
    elif activation == 'elu':
        activation_function = tf.nn.elu
    elif activation == 'sigmoid':
        activation_function = tf.sigmoid
    elif activation == 'tanh':
        activation_function = tf.tanh
    elif activation == 'silu':
        activation_function = silu
    else:
        raise ValueError('invalid activation function')
    z,encoder_weights,encoder_biases = build_network_layers(x, input_dim, latent_dim, widths, 
                                                            activation_function, 'encoder', model_params)
    
    x_decode,decoder_weights,decoder_biases = build_network_layers(z, latent_dim, input_dim, widths[::-1], 
                                                                   activation_function, 'decoder', model_params)

    return z, x_decode, encoder_weights, encoder_biases, decoder_weights, decoder_biases


def build_network_layers(input, input_dim, output_dim, widths, activation, name, model_params):
    """
    Construct one portion of the network (either encoder or decoder).

    Arguments:
        input - 2D tensorflow array, input to the network (shape is [?,input_dim])
        input_dim - Integer, number of state variables in the input to the first layer
        output_dim - Integer, number of state variables to output from the final layer
        widths - List of integers representing how many units are in each network layer
        activation - Tensorflow function to be used as the activation function at each layer
        name - String, prefix to be used in naming the tensorflow variables

    Returns:
        input - Tensorflow array, output of the network layers (shape is [?,output_dim])
        weights - List of tensorflow arrays containing the network weights
        biases - List of tensorflow arrays containing the network biases
    """
    weights = []
    biases = []
    last_width = input_dim # output width of last (previous) layer
    
    if len(model_params) > 0: # build layers with existing coefficients
        # build hidden layers
        for i,n_units in enumerate(widths):
            if name == 'encoder':
                W = tf.get_variable(name+'_W'+str(i), initializer=model_params[1][i])
                b = tf.get_variable(name+'_b'+str(i), initializer=model_params[2][i])
            elif name == 'decoder':
                W = tf.get_variable(name+'_W'+str(i), initializer=model_params[3][i])
                b = tf.get_variable(name+'_b'+str(i), initializer=model_params[4][i])
            last_width = n_units
            weights.append(W)
            biases.append(b)

        # build last layer
        if name == 'encoder':
            W = tf.get_variable(name+'_W'+str(len(widths)), initializer=model_params[1][-1])
            b = tf.get_variable(name+'_b'+str(len(widths)), initializer=model_params[2][-1])
        elif name == 'decoder':
            W = tf.get_variable(name+'_W'+str(len(widths)), initializer=model_params[3][-1])
            b = tf.get_variable(name+'_b'+str(len(widths)), initializer=model_params[4][-1])
        weights.append(W)
        biases.append(b)
    
    else: # build layers with randomly initialized coefficients
        # build hidden layers
        for i,n_units in enumerate(widths):
            W = tf.get_variable(name+'_W'+str(i), shape=[last_width,n_units])
            b = tf.get_variable(name+'_b'+str(i), shape=[n_units], initializer=tf.constant_initializer(0.0))
            last_width = n_units
            weights.append(W)
            biases.append(b)
        
        # build last layer
        W = tf.get_variable(name+'_W'+str(len(widths)), shape=[last_width,output_dim])
        b = tf.get_variable(name+'_b'+str(len(widths)), shape=[output_dim],initializer=tf.constant_initializer(0.0))
        weights.append(W)
        biases.append(b)
    
    # forward pass
    for i in range(len(weights)-1):
        input = tf.matmul(input, weights[i]) + biases[i]
        if activation is not None:
            input = activation(input)      
    input = tf.matmul(input, weights[-1]) + biases[-1] # last layer
    
    return input, weights, biases


def sindy_library_tf(z, latent_dim, poly_order, include_sine=False, include_cosine=False):
    """
    Build the SINDy library.

    Arguments:
        z - 2D tensorflow array of the snapshots on which to build the library. Shape is number of
        time points by the number of state variables.
        latent_dim - Integer, number of state variable in z.
        poly_order - Integer, polynomial order to which to build the library. Max value is 5.
        include_sine - Boolean, whether or not to include sine terms in the library. Default False.

    Returns:
        2D tensorflow array containing the constructed library. Shape is number of time points by
        number of library functions. The number of library functions is determined by the number
        of state variables of the input, the polynomial order, and whether or not sines are included.
    """
    library = [tf.ones(tf.shape(z)[0])]
    
    if poly_order > 0:
        for i in range(latent_dim):
            library.append(z[:,i])
    
    if poly_order > 1:
        for i in range(latent_dim):
            for j in range(i,latent_dim):
                library.append(tf.multiply(z[:,i], z[:,j]))
    
    if poly_order > 2:
        for i in range(latent_dim):
            for j in range(i,latent_dim):
                for k in range(j,latent_dim):
                    library.append(z[:,i]*z[:,j]*z[:,k])

    if poly_order > 3:
        for i in range(latent_dim):
            for j in range(i,latent_dim):
                for k in range(j,latent_dim):
                    for p in range(k,latent_dim):
                        library.append(z[:,i]*z[:,j]*z[:,k]*z[:,p])

    if poly_order > 4:
        for i in range(latent_dim):
            for j in range(i,latent_dim):
                for k in range(j,latent_dim):
                    for p in range(k,latent_dim):
                        for q in range(p,latent_dim):
                            library.append(z[:,i]*z[:,j]*z[:,k]*z[:,p]*z[:,q])

    if include_sine:
        for i in range(latent_dim):
            library.append(tf.sin(z[:,i]))
            
    if include_cosine:
        for i in range(latent_dim):
            library.append(tf.cos(z[:,i]))
            
    return tf.stack(library, axis=1)


def sindy_library_tf_order2(z, dz, latent_dim, poly_order, include_sine=False, include_cosine=False):
    """
    Build the SINDy library for a second order system. This is essentially the same as for a first
    order system, but library terms are also built for the derivatives.
    """
    library = [tf.ones(tf.shape(z)[0])]

    z_combined = tf.concat([z, dz], 1)

    if poly_order > 0:
        for i in range(2*latent_dim):
            library.append(z_combined[:,i])

    if poly_order > 1:
        for i in range(2*latent_dim):
            for j in range(i,2*latent_dim):
                library.append(tf.multiply(z_combined[:,i], z_combined[:,j]))

    if poly_order > 2:
        for i in range(2*latent_dim):
            for j in range(i,2*latent_dim):
                for k in range(j,2*latent_dim):
                    library.append(z_combined[:,i]*z_combined[:,j]*z_combined[:,k])

    if poly_order > 3:
        for i in range(2*latent_dim):
            for j in range(i,2*latent_dim):
                for k in range(j,2*latent_dim):
                    for p in range(k,2*latent_dim):
                        library.append(z_combined[:,i]*z_combined[:,j]*z_combined[:,k]*z_combined[:,p])

    if poly_order > 4:
        for i in range(2*latent_dim):
            for j in range(i,2*latent_dim):
                for k in range(j,2*latent_dim):
                    for p in range(k,2*latent_dim):
                        for q in range(p,2*latent_dim):
                            library.append(z_combined[:,i]*z_combined[:,j]*z_combined[:,k]*z_combined[:,p]*z_combined[:,q])

    if include_sine:
        for i in range(2*latent_dim):
            library.append(tf.sin(z_combined[:,i]))

    if include_cosine:
        for i in range(2*latent_dim):
            library.append(tf.cos(z_combined[:,i]))
            
    return tf.stack(library, axis=1)
