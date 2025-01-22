import h5py
import anndata
import numpy as np
import argparse
from scipy.sparse import csc_matrix
import os


parser = argparse.ArgumentParser(description='ConvertAnndata')
parser.add_argument('--dataDir', type=str, nargs=1, default='/home/'+os.environ['USER']+'/DeepSpaCE/data', required=True,
                    help='Data directory (default: '+'/home/'+os.environ['USER']+'/DeepSpaCE/data'+')')

args = parser.parse_args()
dataDir = args.dataDir[0]

input_file = dataDir+"cleaned_feature_bc_matrix.h5ad"
output_file = dataDir+"cleaned_feature_bc_matrix.h5"

adata = anndata.read_h5ad(input_file)
X_transposed = adata.X.transpose()

# CSC format check
X_csc = csc_matrix(X_transposed)

num_features, num_barcodes = X_csc.shape

# Write to HDF5 file using 10x format
with h5py.File(output_file, "w") as f:
    
    matrix_grp = f.create_group("matrix")
    
    barcodes = adata.obs.index.to_numpy().astype('S')  # Assuming obs contains barcodes
    matrix_grp.create_dataset("barcodes", data=barcodes)

    features_grp = matrix_grp.create_group("features")
    features_grp.create_dataset("id", data=adata.var.index.to_numpy().astype('S'))
    features_grp.create_dataset("name", data=adata.var["name"].to_numpy().astype('S') if "name" in adata.var else adata.var.index.astype('S'))
    features_grp.create_dataset("feature_type", data=np.array(["Gene Expression"] * num_features, dtype='S'))
    features_grp.create_dataset("genome", data=np.array(["GRCh38"] * num_features, dtype='S'))  # Placeholder for genome

    matrix_grp.create_dataset("data", data=X_csc.data)
    matrix_grp.create_dataset("indices", data=X_csc.indices)
    matrix_grp.create_dataset("indptr", data=X_csc.indptr)
    matrix_grp.create_dataset("shape", data=np.array([num_features, num_barcodes], dtype=np.int64))  # Feature x barcode

