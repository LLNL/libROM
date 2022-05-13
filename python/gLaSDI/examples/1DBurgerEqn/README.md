Procedure:
- run `generate_data.ipynb` to generate dataset; Note: create a folder (named as "data") to store the data
- run `train_model_sec4-1_case2.ipynb` to train a gLaSDI model by submitting a job using `batch_job.bsub`; Note: create a folder (named as "fig") to store model parameters and figures (next step)
- run `test_model_v2.ipynb` to evalute the trained model
- run `plot_heatmap.ipynb` to plot loss history and the greedy sampling process in ternms of the error indicator


If the training time exceeds the maximum allowable wall time, the training can be continued by updating the following parameters in `train_model_sec4-1_case2.ipynb` and submitting additoinal jobs:
- `params['retrain'] = True` in cell [8]
- `save_name` (file name of the trained model) in cell [9]