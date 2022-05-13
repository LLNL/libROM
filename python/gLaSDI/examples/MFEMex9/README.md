The [ex9.cpp](https://github.com/mfem/mfem/blob/master/examples/ex9.cpp) is modified to include the parameterized initial condition and an option to compute the residual given the predicted solutions from gLaSDI. The modified `ex9.cpp` and the corresponding executable file (`ex9`) have been uploaded to `gLaSDI/src/`. The executable (`ex9`) will be called to compute residual-based error indicator during training of gLaSDI.

Procedure:
- Create the following folders:
    - "data": to store the training data
    - "fig": to store model parameters and figures
    - "temp": to store temporary files created during training, for evaluation of the residual-based error indicator
- run `generate_data.ipynb` in the `gLaSDI/MFEMex9/generate_data/` folder to generate dataset
- run `train_model_sec4-4_case1.ipynb` to train a gLaSDI model by submitting a job using `batch_job.bsub`
- run `test_model.ipynb` to evalute the trained model


If the training time exceeds the maximum allowable wall time, the training can be continued by updating the following parameters in `train_model_sec4-4_case1.ipynb` and submitting additoinal jobs:
- `params['retrain'] = True` in cell [7]
- `save_name` (file name of the trained model) in cell [8]