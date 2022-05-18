The gLaSDI framework is applied to construct a non-intrusive reduced-order model of the [**MFEM Example 9 Time-Dependent Radial Advection**](https://github.com/mfem/mfem/blob/master/examples/ex9.cpp).


The data can be generated by `generate_data.ipynb`. Alternatively, some datasets are avaiable in `/usr/workspace/he10/data/MFEMex9/`, which are used in the numerical examples in the paper. The datasets are generated from the parameter space constituted by the parameters of the initial condition, including **w1** and **w2**.


For w1 x w2 = [1.5,1.8] x [2.0,2.3]:
- `local4_tstop3.0b.p`: 2x2 parameter cases on the given parameter space; initial parameter cases located at the corners of the parameter space
- `local121_tstop3.0b.p`: 11x11 parameter cases on the given parameter space; a dataset of discrete parameter space for adaptive greedy sampling
- `local441_tstop3.0b.p`: 21x21 parameter cases on the given parameter space; a dataset of discrete parameter space for evaluation


For w1 x w2 = [1.5,2.0] x [2.0,2.5]:
- `local4_tstop3.0.p`: 2x2 parameter cases on the given parameter space; initial parameter cases located at the corners of the parameter space
- `local121_tstop3.0.p`: 11x11 parameter cases on the given parameter space; a dataset of discrete parameter space for adaptive greedy sampling
- `local441_tstop3.0.p`: 21x21 parameter cases on the given parameter space; a dataset of discrete parameter space for evaluation


Note that to enhance the training efficiency, the dataset of a discrete parameter space (e.g., 11x11 parameter cases) for adaptive greedy sampling is generated before training starts, which means the solution of the optimal parameter case (with the maximum error indicator) is retrieved from the dataset without calling the solver of the full-order model during training. The proposed gLaSDI framework can be adapted to a continuous parameter space and the solution of sampled parameter cases can be generated by running full-order model simulations during training.


The [ex9.cpp](https://github.com/mfem/mfem/blob/master/examples/ex9.cpp) is modified to include the parameterized initial condition and an option to compute the residual given the predicted solutions from gLaSDI. The modified `ex9.cpp` and the corresponding executable file (`ex9`) are located at `gLaSDI/src/`. The executable (`ex9`) will be called to compute residual-based error indicator during training of gLaSDI.


Procedure:
- run `generate_data.ipynb` to generate dataset
- run `train_model.ipynb` to train a gLaSDI model by submitting a job using `batch_job.bsub`
- run `test_model.ipynb` to evalute the trained model


If the training time exceeds the maximum allowable wall time, the training can be continued by updating the following parameters in `train_model.ipynb` and submitting additoinal jobs:
- `params['retrain'] = True`
- `save_name` (file name of the trained model) 