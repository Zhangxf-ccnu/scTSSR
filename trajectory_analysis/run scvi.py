# -*- coding: utf-8 -*-

# Introduction to single-cell Variational Inference (scVI)---------------------------
# %matplotlib inline
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scvi.dataset import CortexDataset, RetinaDataset
from scvi.models import VAE
from scvi.inference import UnsupervisedTrainer
import torch


# Loading a local dataset 
from scvi.dataset import CsvDataset
save_path = "scTSSR-data/"
local_csv_dataset = CsvDataset("deng_input_for_scvi.csv", 
                               save_path=save_path,
                               new_n_genes=13514)

# Training
n_epochs = 400 

#if n_epochs_all is None else n_epochs_all
lr = 1e-3
use_batches = False
use_cuda = False

vae = VAE(n_input = local_csv_dataset.nb_genes, n_batch=0)
trainer = UnsupervisedTrainer(
    vae,
    local_csv_dataset,
    train_size=0.75,
    use_cuda=use_cuda,
    frequency=5,
)


trainer.train(n_epochs=n_epochs, lr=lr)
torch.save(trainer.model.state_dict(), '%s/vae.pkl' % save_path)
    
    
full = trainer.create_posterior(trainer.model, local_csv_dataset, indices=np.arange(len(local_csv_dataset)))

imputed_values = full.sequential().imputation()   
normalized_values = full.sequential().get_sample_scale()  

# save imputed result as csv
np.savetxt("deng_scvi.csv", imputed_values.T, delimiter=",")

