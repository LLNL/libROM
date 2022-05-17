# gLaSDI: Parametric Physics-informed Greedy Latent Space Dynamics Identification
- In the proposed gLaSDI framework, an autoencoder discovers intrinsic nonlinear latent representations of high-dimensional data, while dynamics identification (DI) models capture local latent-space dynamics.  
- An interactive training algorithm is adopted for the autoencoder and local DI models, which enables identification of simple latent-space dynamics and enhances accuracy and efficiency of data-driven reduced-order modeling. 
- To maximize and accelerate the exploration of the parameter space for the optimal model performance, an adaptive greedy sampling algorithm integrated with a physics-informed residual-based error indicator and random-subset evaluation is introduced to search for the optimal training samples on-the-fly.
- Further, to exploit local latent-space dynamics captured by the local DI models for an improved modeling accuracy with a minimum number of local DI models in the parameter space, an efficient k-nearest neighbor convex interpolation scheme is employed.
- The effectiveness of the proposed framework is demonstrated by modeling various nonlinear dynamical problems, including Burgers equations, nonlinear heat conduction, and radial advection. 
- The proposed adaptive greedy sampling outperforms the conventional predefined uniform sampling in terms of accuracy. Compared with the high-fidelity models, gLaSDI achieves 66 to 4,417x speed-up with 1 to 5% relative errors.

## Install (LC Lassen)
(1) The required pakcages are large, so we need to first configure the `conda` directory so that the packages are downloaded and installed in the workspace (replace `UserName` with your OUN):
`conda config --add pkgs_dirs /usr/workspace/{UserName}/cache/pkgs`

(2) Run `bash setup.sh` to create a TensorFlow (GPU-version) virtual environment.

(3) Activate the virtual environment: `conda activate ~/.conda/envs/tfvenv`.

(4) If Step (2）doesn't install TensorFlow successfually, run `conda install -y tensorflow-gpu`

## Examples
Four examples are provided, including **1D Burgers Equation**, **2D Burgers Equation**, **Time Dependent Heat Conduction (MFEM Example 16)**, and **Radial Advection (MFEM Example 9)**. The notebooks for data generation, model training, model evaluation are provided in the `examples` folder.
