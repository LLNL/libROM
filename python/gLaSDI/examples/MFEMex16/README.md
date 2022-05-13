The [ex16.cpp](https://github.com/mfem/mfem/blob/master/examples/ex16.cpp) is modified to include the parameterized initial condition and an option to compute the residual given the predicted solutions from gLaSDI. The modified `ex16.cpp` and the corresponding executable file (`ex16`) have been uploaded to `gLaSDI/src/`. The executable (`ex16`) will be called to compute residual-based error indicator during training of gLaSDI.

Procedure:
- Create the following folders:
    - "data": to store the training data
    - "fig": to store model parameters and figures
    - "temp": to store temporary files created during training, for evaluation of the residual-based error indicator
- run `generate_data.ipynb` in the `gLaSDI/MFEMex16/generate_data/` folder to generate dataset
- run `train_model_sec4-3_case2.ipynb` to train a gLaSDI model by submitting a job using `batch_job.bsub`
- run `test_model.ipynb` to evalute the trained model


If the training time exceeds the maximum allowable wall time, the training can be continued by updating the following parameters in `train_model_sec4-3_case2.ipynb` and submitting additoinal jobs:
- `params['retrain'] = True` in cell [9]
- `save_name` (file name of the trained model) in cell [10]