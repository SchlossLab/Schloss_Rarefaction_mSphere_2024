#!/usr/bin/env python3

import sys
import numpy as np
from biom import Table
import pandas as pd
from gemelli.preprocessing import rclr_transformation
from gemelli.rpca import rpca

np.seterr(divide = 'ignore') # for rclr

shared_file = sys.argv[1]
dist_file = sys.argv[2]

shared_transposed_df = pd.read_csv(shared_file,
                                   sep = "\t")
#shared_transposed_df.index.name = 'otu' # format for R 

# convert to biom-formatted file
shared_transposed_biom = Table(shared_transposed_df.to_numpy(),
                               shared_transposed_df.index,
                               shared_transposed_df.columns) 

rclr_gemelli = rclr_transformation(shared_transposed_biom) # compute RCLR

# Run RPCA
ordination, distance_matrix = rpca(shared_transposed_biom, n_components = 3) 

distance_matrix.to_data_frame().to_csv(dist_file, sep='\t') # save just the distances since that is what is output from the R version
