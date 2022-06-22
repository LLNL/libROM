## gLaSDI: Parametric Physics-informed Greedy Latent Space Dynamics Identification
- In the proposed gLaSDI framework, an autoencoder discovers intrinsic nonlinear latent representations of high-dimensional data, while dynamics identification (DI) models capture local latent-space dynamics.  
- An interactive training algorithm is adopted for the autoencoder and local DI models, which enables identification of simple latent-space dynamics and enhances accuracy and efficiency of data-driven reduced-order modeling. 
- To maximize and accelerate the exploration of the parameter space for the optimal model performance, an adaptive greedy sampling algorithm integrated with a physics-informed residual-based error indicator and random-subset evaluation is introduced to search for the optimal training samples on-the-fly.
- Further, to exploit local latent-space dynamics captured by the local DI models for an improved modeling accuracy with a minimum number of local DI models in the parameter space, an efficient k-nearest neighbor convex interpolation scheme is employed.
- The effectiveness of the proposed framework is demonstrated by modeling various nonlinear dynamical problems, including Burgers equations, nonlinear heat conduction, and radial advection. 
- The proposed adaptive greedy sampling outperforms the conventional predefined uniform sampling in terms of accuracy. Compared with the high-fidelity models, gLaSDI achieves 66 to 4,417x speed-up with 1 to 5% relative errors.


## Required Packages
The following versions of packages have been verified to work. Other versions may also work.
- Python: 3.7.10
- TensorFlow: 2.2.0
- Numpy: 1.17.4
- Scipy: 1.4.1
- Sklearn: 0.23.2
- Pandas: 1.1.3
- Matplotlib: 3.4.2
- Seaborn: 0.11.0
- Pickle: 0.7.5

*For LC users*:
1. Configure the `conda` directory to the workspace:
`conda config --add pkgs_dirs /usr/workspace/{OUN}/cache/pkgs`
2. Create a TensorFlow (GPU-version) virtual environment: `bash setup.sh`
3. Activate the TensorFlow virtual environment: `conda activate ~/.conda/envs/tfvenv`
4. If TensorFlow wasn't installed successfually in Step 2, try: `conda install -y tensorflow-gpu`


## Examples
Three examples are provided, including 
- 1D Burgers Equation 
- 2D Burgers Equation
- Time-Dependent Radial Advection ([MFEM Example 9](https://github.com/mfem/mfem/blob/master/examples/ex9.cpp)).

The Python scripts for data generation, model training and evaluation are provided in `gLaSDI/examples/`. 


## Description of Parameters:
- `seed` - int, seed for random number generators; To ensure reproducibility, use the same seed.
- `config` - configuration of TensorFlow; see `tensorflow.ConfigProto`
- `num_sindy` - int, the number of existing local dynamics identification (DI) models in the parameter space
- `param` - the parameter set of sampled parameter cases
- `train_idx` - list, the indices of sampled parameter cases in the discrete parameter space
- `input_dim` - int, input dimension of the auto-encoder
- `latent_dim` - int, latent-space dimension of the auto-encoder
- `poly_order` - int, from 1-5, maximum polynomial order to which to build the DI library 
- `include_sine` - bool, whether or not to include sine functions in the DI library
- `include_cosine` - bool, whether or not to include cosine functions in the DI library
- `include_constant` - bool, whether or not to include constant in the DI library
- `library_dim` - int, total number of basis functions; this is determined based on the `latent_dim`, `model_order`, `poly_order`, `include_sine`, `include_cosine`, and `include_constant`, and can be calculated using the function `library_size` in `sindy_utils.py`
- `sequential_thresholding` - bool, whether or not to perform sequential thresholding on the DI coefficient matrix
- `coefficient_threshold` - float, minimum magnitude of coefficients to keep in the DI coefficient matrix when performing thresholding
- `threshold_frequency` - int, the number of epochs after which to perform thresholding
- `coefficient_mask` - numpy.array, matrix of ones and zeros that determines which coefficients are still included in the DI model; typically initialized to all ones and will be modified by the sequential thresholding procedure
- `coefficient_initialization` - str, how to initialize the DI coefficient matrix; options are 'constant' (initialize as all 1s), 'xavier' (initialize using the xavier initialization approach), 'specified' (pass in an additional parameter init_coefficients that has the values to use to initialize the DI coefficient matrix)
- `loss_weight_decoder` - float, weighting of the auto-encoder reconstruction in the loss function (should keep this at 1.0 and adjust the other weightings proportionally)
- `loss_weight_sindy_z`- float, weighting of the DI prediction in the latent space in the loss function
- `loss_weight_sindy_x` - float, weighting of the DI prediction passed back to the input space in the loss function
- `loss_weight_sindy_regularization` - float, weighting of the L1 regularization on the DI coefficients in the loss function; default: zero
- `diff` - str, 'symb': symbolic differentiation (only for fully connected Autoencoder), 'auto': automatic differentiation; default: 'symb'
- `activation` - str, activation function to be used in the network; options are 'sigmoid', 'relu', 'linear', or 'elu'
- `widths` - list of integers specifying the number of units for each hidden layer of the encoder; decoder widths will be the reverse order of these widths
- `epoch_size` - int, the number of training samples in an epoch
- `batch_size` - int, the number of samples to use in a batch of training
- `learning rate` - float, initial learning rate passed to the Adam optimizer
- `fig_path` - str, path specifying where to save the resulting models
- `print_progress` - bool, whether or not to print updates during training
- `print_frequency` - int, print progress at intervals of this many epochs
- `max_epochs` - int, the maximum number of training epochs
- `update_epoch` - int, greedy sampling frequency, update training set at intervals of this many epochs
- `tol` - float, initial tolerance of the maximum error indicator in the parameter space; it will be updated during training using the prescribed `adaptive` method
- `tol2` - float, initial tolerance of the maximum relative error in the parameter space
- `sindy_max` - int or `None`, the maximum number of local DIs; if tolerance is used as a termination criterion; set it as `None`
- `convex_knn` - int, the number nearest local DIs used for convex interpolation of coefficient matrices during Greedy sampling
- `test_data` - dict, dataset of the discrete parameter space
- `test_param` - numpy.array, parameters of the discrete parameter space
- `num_test` - int, the number of parameter cases of the discrete parameter space
- `coeff_exist` - bool, whether or not to initialize model coefficients with pescribed values; set as `False`
- `retrain` - bool, whether or not to retrain the model; set as `False` for training a new model
- `err_type` - int, 1: max relative error (if test data is available); 2: residual norm for 1D Burgers eqn; 3: residual norm for 2D Burgers eqn; 4: residual norm for time dependent heat conduction (MFEM example 16); 5: residual norm for radial advection (MFEM example 9)
- `subsize` - int, initial random subset size, the number of randomly selected cases for Greedy sampling
- `subsize_max` - int, maximum random subset size in percentage
- `adapative` - str, the method used to update the error tolerance; 'mean': use mean ratios between error indicator and max relative errors; 'reg_mean': use linear regression line; 'reg_max': use linear regression line shifted by std to upper bound; 'reg_min': use linear regression line shifted by std to lower bound, more conservative; 
- `pde` - dict, stores the parameters related to the PDE


## Questions/Comments
Questions and comments should be directed to xih251@eng.uscd.edu