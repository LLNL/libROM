# import tensorflow as tf
import tensorflow.compat.v1 as tf
tf.compat.v1.disable_v2_behavior()
from tensorflow.python.ops.parallel_for.gradients import jacobian, batch_jacobian
from time import time
import pickle
import numpy as np

def full_network(params):
    """
    Define the full network architecture.

    Arguments:
        params - Dictionary object containing the parameters that specify the training.
        See README file for a description of the parameters.

    Returns:
        network - Dictionary containing the tensorflow objects that make up the network.
    """
    input_dim = params['input_dim']
    latent_dim = params['latent_dim']
    activation = params['activation']
    poly_order = params['poly_order']
    num_sindy = params['num_sindy'] # number of local SINDys
    model_params = []
    if params['coeff_exist']:
        model_params = pickle.load(open(params['fig_path'] + params['save_name'] + '_params.pkl', 'rb'))['model_params']
        
    if 'include_sine' in params.keys():
        include_sine = params['include_sine']
    else:
        include_sine = False
        
    if 'include_cosine' in params.keys():
        include_cosine = params['include_cosine']
    else:
        include_cosine = False
        
    library_dim = params['library_dim']
    model_order = params['model_order']

    network = {}

    x = tf.placeholder(tf.float32, shape=[None, num_sindy, input_dim], name='x')
    dx = tf.placeholder(tf.float32, shape=[None, num_sindy, input_dim], name='dx')
    if model_order == 2:
        ddx = tf.placeholder(tf.float32, shape=[None, num_sindy, input_dim], name='ddx')

    # Autoencoder
    if activation == 'linear':
        z, x_decode, encoder_weights, encoder_biases, decoder_weights, decoder_biases = linear_autoencoder(x, input_dim, latent_dim, model_params)
    else:
        z, x_decode, encoder_weights, encoder_biases, decoder_weights, decoder_biases = nonlinear_autoencoder(x, input_dim, latent_dim,
                                                                                                              params['widths'],
                                                                                                              model_params,
                                                                                                              activation=activation)

        
    # compute dz/dt
    if params['diff'] == 'symb': # symbolic differentiation
        if model_order == 1:
            dz = z_derivative(x, dx, encoder_weights, encoder_biases, activation=activation) # [batch,num_sindy,latent_dim]
        else:
            dz,ddz = z_derivative_order2(x, dx, ddx, encoder_weights, encoder_biases, activation=activation)
    elif params['diff'] == 'auto': # automatic differentiation
        dzdx_batch = batch_jacobian(z, x)
        dzdx = []
        for i in range(num_sindy):
            dzdx.append(dzdx_batch[:,i,:,i,:]) # [batch,output_dim,input_dim]
        dzdx = tf.stack(dzdx, axis=1) # [batch,num_sindy,output_dim]
        dz = tf.matmul(dzdx, dx[:,:,:,None])[:,:,:,0] # [batch,num_sindy,latent_dim]
        
        
    # SINDy
    Theta = []
    sindy_coefficients = []
    sindy_predict = [] # [batch,num_sindy,latent_dim]

        
    if params['sequential_thresholding']: # 1 mask for all local SINDys
        coefficient_mask = tf.placeholder(tf.float32, shape=[library_dim,latent_dim], name='coefficient_mask')
        network['coefficient_mask'] = coefficient_mask
            

    for i in range(num_sindy):
        if model_order == 1:
            Theta.append(sindy_library_tf(z[:,i,:], latent_dim, poly_order, include_sine, include_cosine))
        else:
            Theta.append(sindy_library_tf_order2(z, dz, latent_dim, poly_order, include_sine, include_cosine))
        
        if params['coeff_exist']:
            if i < len(model_params[0]):
                sindy_coefficients.append(tf.get_variable(f'sindy_coefficients{i}', initializer=model_params[0][i]))
            else:
                print(f"  Existing SINDys: {len(model_params[0])}, Create new SINDy: {i+1}")
        
                # initialize the new local SINDy with the coefficients of the nearest SINDy
                all_param = np.stack(params['param'])
                idx = np.argmin(np.linalg.norm(all_param[:-1]-all_param[-1], axis=1))

                sindy_coefficients.append(tf.get_variable(f'sindy_coefficients{i}', initializer=model_params[0][idx]))  
        else:
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
                
    sindy_predict = tf.stack(sindy_predict, axis=1) # [batch,num_sindy,latent_dim]
    
    # compute dx/dt
    if params['loss_weight_sindy_x'] > 0:
        if params['diff'] == 'symb': # symbolic differentiation
            dx_decode = z_derivative(z, sindy_predict, decoder_weights, decoder_biases, activation=activation) # [batch,num_sindy,input_dim]
        elif params['diff'] == 'auto': # automatic differentiation
            dx_decode = tf.linalg.matmul(dxdz_decode, sindy_predict[:,:,:,None])[:,:,:,0] # [batch,num_sindy,input_dim]
    else:
        dx_decode = tf.zeros(tf.shape(dx))
        
        
    network['x'] = x   # [batch,num_sindy,input_dim]
    network['dx'] = dx # [batch,num_sindy,input_dim]
    network['z'] = z   # [batch,num_sindy,latent_dim]
    network['dz'] = dz # [batch,num_sindy,latent_dim]
    network['x_decode'] = x_decode   # [batch,num_sindy,input_dim]
    network['dx_decode'] = dx_decode # [batch,num_sindy,input_dim]
    network['encoder_weights'] = encoder_weights
    network['encoder_biases'] = encoder_biases
    network['decoder_weights'] = decoder_weights
    network['decoder_biases'] = decoder_biases
    network['Theta'] = Theta # list
    network['sindy_coefficients'] = sindy_coefficients # list
    
    if model_order == 1:
        network['dz_predict'] = sindy_predict # [batch,num_sindy,latent_dim]
    else:
        network['ddz'] = ddz
        network['ddz_predict'] = sindy_predict
        network['ddx'] = ddx
        network['ddx_decode'] = ddx_decode

    return network

    
def define_loss(network, params):
    """
    Create the loss functions.

    Arguments:
        network - Dictionary object containing the elements of the network architecture.
        This will be the output of the full_network() function.
    """
    x = network['x']                # [batch,num_sindy,input_dim]
    x_decode = network['x_decode']  # [batch,num_sindy,input_dim]
    num_sindy = params['num_sindy'] # number of local SINDys
    sindy_coefficients = network['sindy_coefficients']
    
    if params['model_order'] == 1:
        dz = network['dz']                 # [batch,num_sindy,latent_dim]
        dz_predict = network['dz_predict'] # [batch,num_sindy,latent_dim]
        dx = network['dx']                 # [batch,num_sindy,input_dim]
        dx_decode = network['dx_decode']   # [batch,num_sindy,input_dim]
    else:
        ddz = network['ddz']
        ddz_predict = network['ddz_predict']
        ddx = network['ddx']
        ddx_decode = network['ddx_decode']
        
    losses = {}
    losses['decoder'] = tf.zeros((1))
    losses['sindy_x'] = tf.zeros((1))
    losses['sindy_z'] = tf.zeros((1))
    losses['sindy_regularization'] = tf.zeros((1))
    
    for i in range(num_sindy):
        losses['decoder'] += tf.reduce_mean((x[:,i,:] - x_decode[:,i,:])**2)
        if params['model_order'] == 1:
            losses['sindy_x'] += tf.reduce_mean((dx[:,i,:] - dx_decode[:,i,:])**2)
            losses['sindy_z'] += tf.reduce_mean((dz[:,i,:] - dz_predict[:,i,:])**2)
        else:
            losses['sindy_x'] = tf.reduce_mean((ddx - ddx_decode)**2)
            losses['sindy_z'] = tf.reduce_mean((ddz - ddz_predict)**2)
        losses['sindy_regularization'] += tf.reduce_mean(tf.abs(sindy_coefficients[i]))
        
    loss = params['loss_weight_decoder'] * losses['decoder'] \
           + params['loss_weight_sindy_z'] * losses['sindy_z'] \
           + params['loss_weight_sindy_x'] * losses['sindy_x'] \
           + params['loss_weight_sindy_regularization'] * losses['sindy_regularization']

    loss_refinement = params['loss_weight_decoder'] * losses['decoder'] \
                      + params['loss_weight_sindy_z'] * losses['sindy_z'] \
                      + params['loss_weight_sindy_x'] * losses['sindy_x']
    
    return loss, losses, loss_refinement


def linear_autoencoder(x, input_dim, d, model_params):
    z,encoder_weights,encoder_biases = build_network_layers(x, input_dim, latent_dim, [], None, 'encoder', model_params)
    x_decode,decoder_weights,decoder_biases = build_network_layers(z, latent_dim, input_dim, [], None, 'decoder', model_params)

    return z, x_decode, encoder_weights, encoder_biases,decoder_weights,decoder_biases


def silu(x):
    return x * tf.sigmoid(x)

  
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
    z,encoder_weights,encoder_biases = build_network_layers(x, input_dim, latent_dim, widths, activation_function, 'encoder', model_params)
    x_decode,decoder_weights,decoder_biases = build_network_layers(z, latent_dim, input_dim, widths[::-1], activation_function, 'decoder', model_params)

    return z, x_decode, encoder_weights, encoder_biases, decoder_weights, decoder_biases


def nonlinear_autoencoder2(x, input_dim, latent_dim, widths, model_params, activation='elu'):
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
    z,dzdx,encoder_weights,encoder_biases = build_network_layers(x, input_dim, latent_dim, widths, activation_function, 'encoder', model_params)
    x_decode,dxdz_decode,decoder_weights,decoder_biases = build_network_layers(z, latent_dim, input_dim, widths[::-1], activation_function, 'decoder', model_params)

    return z, dzdx, x_decode, dxdz_decode, encoder_weights, encoder_biases, decoder_weights, decoder_biases


def nonlinear_autoencoder3(x, dx, input_dim, latent_dim, widths, model_params, activation='elu'):
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
    z,dzdx,encoder_weights,encoder_biases = build_network_layers(x, input_dim, latent_dim, widths, activation_function, 'encoder', model_params)
    x_decode,dxdz_decode,decoder_weights,decoder_biases = build_network_layers(z, latent_dim, input_dim, widths[::-1], activation_function, 'decoder', model_params)
    
    dz = tf.matmul(dzdx, dx[:,:,:,None])[:,:,:,0] # [batch,num_sindy,latent_dim]
    return z, dz, x_decode, dxdz_decode, encoder_weights, encoder_biases, decoder_weights, decoder_biases


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
    
    if len(model_params) > 0:
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
    
    else:
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
    output = []
    for j in range(input.shape[1]):
        input_j = input[:,j,:]
        for i in range(len(weights)-1):
            input_j = tf.matmul(input_j, weights[i]) + biases[i]
            if activation is not None:
                input_j = activation(input_j)
        output.append(tf.matmul(input_j, weights[-1]) + biases[-1]) # last layer, [batch,output_dim]
    output = tf.stack(output, axis=1) # [batch,num_sindy,output_dim]
    return output, weights, biases


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


def z_derivative(input, dx, weights, biases, activation='elu'):
    """
    Compute the first order time derivatives by propagating through the network.

    Arguments:
        input - 2D tensorflow array, input to the network. Dimensions are number of time points
        by number of state variables.
        dx - First order time derivatives of the input to the network.
        weights - List of tensorflow arrays containing the network weights
        biases - List of tensorflow arrays containing the network biases
        activation - String specifying which activation function to use. Options are
        'elu' (exponential linear unit), 'relu' (rectified linear unit), 'sigmoid',
        or linear.

    Returns:
        dz - Tensorflow array, first order time derivatives of the network output.
    """
    num_sindy = input.shape[1]
    dz = []
    
    if activation == 'elu':
        for j in range(num_sindy):
            input_j = input[:,j,:]
            dz_j = dx[:,j,:]
            for i in range(len(weights)-1):
                input_j = tf.matmul(input_j, weights[i]) + biases[i]
                input_j = tf.nn.elu(input_j)
                dz_j = tf.multiply(tf.minimum(tf.exp(input_j),1.0),
                                      tf.matmul(dz_j, weights[i]))
            dz.append(tf.matmul(dz_j, weights[-1])) # [batch,output_dim]
        dz = tf.stack(dz, axis=1) # [batch,num_sindy,output_dim]
        
    elif activation == 'relu':
        for j in range(num_sindy):
            input_j = input[:,j,:]
            dz_j = dx[:,j,:]
            for i in range(len(weights)-1):
                input_j = tf.matmul(input_j, weights[i]) + biases[i]
                input_j = tf.nn.relu(input_j)
                dz_j = tf.multiply(tf.to_float(input_j > 0), tf.matmul(dz_j, weights[i]))
            dz.append(tf.matmul(dz_j, weights[-1])) # [batch,output_dim]
        dz = tf.stack(dz, axis=1) # [batch,num_sindy,output_dim]
        
    elif activation == 'sigmoid':
        for j in range(num_sindy):
            input_j = input[:,j,:]
            dz_j = dx[:,j,:]
            for i in range(len(weights)-1):
                input_j = tf.matmul(input_j, weights[i]) + biases[i]
                input_j = tf.sigmoid(input_j)
                dz_j = tf.multiply(tf.multiply(input_j, 1-input_j), tf.matmul(dz_j, weights[i]))
            dz.append(tf.matmul(dz_j, weights[-1])) # [batch,output_dim]
        dz = tf.stack(dz, axis=1) # [batch,num_sindy,output_dim]
            
    else:
        for j in range(num_sindy):
            input_j = input[:,j,:]
            dz_j = dx[:,j,:]
            for i in range(len(weights)-1):
                input_j = tf.matmul(input_j, weights[i]) + biases[i]
            dz.append(tf.matmul(dz_j, weights[-1])) # [batch,output_dim]
        dz = tf.stack(dz, axis=1) # [batch,num_sindy,output_dim]
    return dz


def z_derivative_order2(input, dx, ddx, weights, biases, activation='elu'):
    """
    Compute the first and second order time derivatives by propagating through the network.

    Arguments:
        input - 2D tensorflow array, input to the network. Dimensions are number of time points
        by number of state variables.
        dx - First order time derivatives of the input to the network.
        ddx - Second order time derivatives of the input to the network.
        weights - List of tensorflow arrays containing the network weights
        biases - List of tensorflow arrays containing the network biases
        activation - String specifying which activation function to use. Options are
        'elu' (exponential linear unit), 'relu' (rectified linear unit), 'sigmoid',
        or linear.

    Returns:
        dz - Tensorflow array, first order time derivatives of the network output.
        ddz - Tensorflow array, second order time derivatives of the network output.
    """
    dz = dx
    ddz = ddx
    if activation == 'elu':
        for i in range(len(weights)-1):
            input = tf.matmul(input, weights[i]) + biases[i]
            dz_prev = tf.matmul(dz, weights[i])
            elu_derivative = tf.minimum(tf.exp(input),1.0)
            elu_derivative2 = tf.multiply(tf.exp(input), tf.to_float(input<0))
            dz = tf.multiply(elu_derivative, dz_prev)
            ddz = tf.multiply(elu_derivative2, tf.square(dz_prev)) \
                  + tf.multiply(elu_derivative, tf.matmul(ddz, weights[i]))
            input = tf.nn.elu(input)
        dz = tf.matmul(dz, weights[-1])
        ddz = tf.matmul(ddz, weights[-1])
        
    elif activation == 'relu':
        # NOTE: currently having trouble assessing accuracy of 2nd derivative due to discontinuity
        for i in range(len(weights)-1):
            input = tf.matmul(input, weights[i]) + biases[i]
            relu_derivative = tf.to_float(input>0)
            dz = tf.multiply(relu_derivative, tf.matmul(dz, weights[i]))
            ddz = tf.multiply(relu_derivative, tf.matmul(ddz, weights[i]))
            input = tf.nn.relu(input)
        dz = tf.matmul(dz, weights[-1])
        ddz = tf.matmul(ddz, weights[-1])
        
    elif activation == 'sigmoid':
        for i in range(len(weights)-1):
            input = tf.matmul(input, weights[i]) + biases[i]
            input = tf.sigmoid(input)
            dz_prev = tf.matmul(dz, weights[i])
            sigmoid_derivative = tf.multiply(input, 1-input)
            sigmoid_derivative2 = tf.multiply(sigmoid_derivative, 1 - 2*input)
            dz = tf.multiply(sigmoid_derivative, dz_prev)
            ddz = tf.multiply(sigmoid_derivative2, tf.square(dz_prev)) \
                  + tf.multiply(sigmoid_derivative, tf.matmul(ddz, weights[i]))
        dz = tf.matmul(dz, weights[-1])
        ddz = tf.matmul(ddz, weights[-1])
        
    else:
        for i in range(len(weights)-1):
            dz = tf.matmul(dz, weights[i])
            ddz = tf.matmul(ddz, weights[i])
        dz = tf.matmul(dz, weights[-1])
        ddz = tf.matmul(ddz, weights[-1])
    return dz,ddz
