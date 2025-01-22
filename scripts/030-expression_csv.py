import numpy as np
import pandas as pd
import scanpy as sc
import argparse
import os
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='ConvertAnndata')
parser.add_argument('--spacerangerDir', type=str, nargs=1, default='/home/'+os.environ['USER']+'/visium/results/counts', required=True,
                    help='Space Ranger outs directory (default: '+'/home/'+os.environ['USER']+'/visium/results/counts'+')')
parser.add_argument('--spotcleanDir', type=str, nargs=1, default='/home/'+os.environ['USER']+'/visium/results/cleaned', required=True,
                    help='SpotClean outs directory (default: '+'/home/'+os.environ['USER']+'/visium/results/cleaned'+')')
parser.add_argument('--outDir', type=str, nargs=1, default='/home/'+os.environ['USER']+'/visium/results/plots', required=True,
                    help='Plotting output directory (default: '+'/home/'+os.environ['USER']+'/visium/results/plots'+')')
parser.add_argument('-r','--radii', nargs='+', default='', 
                    help='Decontamination radii for SpotClean')
parser.add_argument('-g','--gene', type=str, nargs=1, default='INS', 
                    help='Gene for expression plots')
parser.add_argument('-k','--kernels', nargs='+', default='gaussian',
                    help='Kernels to model decontamination')

args = parser.parse_args()
spacerangerDir = args.spacerangerDir[0]
spotcleanDir = args.spotcleanDir[0]
outDir = args.outDir[0]
radii = args.radii
gene = args.gene[0]
kernels = args.kernels

spaceranger_adata = sc.read_10x_h5(spacerangerDir+'filtered_feature_bc_matrix.h5')

# Set up coordinates based on barcodes in spatial output of Space Ranger run
spatial = pd.read_csv(spacerangerDir+'spatial/tissue_positions.csv', header=None)

x = []
y = []

spaceranger_adata.obs[0] = spaceranger_adata.obs.index
spaceranger_adata.obs = spaceranger_adata.obs.merge(spatial, how = 'inner', on = 0)

x = list(spaceranger_adata.obs[4])
y = list(spaceranger_adata.obs[5])


# Space Ranger data
spaceranger_adata.layers['counts'] = spaceranger_adata.X
spaceranger_adata.layers['log1p_counts'] =np.log(1+pd.DataFrame.sparse.from_spmatrix(spaceranger_adata.layers['counts']))

ins_exp_spaceranger = spaceranger_adata.layers['log1p_counts'].T[list(spaceranger_adata.var.index).index(gene)]
colors_spaceranger = pd.DataFrame.sparse.from_spmatrix(ins_exp_spaceranger)
df = pd.DataFrame([x,y,list(colors_spaceranger.loc[0])]).T
df.to_csv(outDir+'spot_'+gene+'_expression.csv', index=False, header=['x','y','expr'])

# Space Ranger-normalized data
sc.pp.normalize_total(spaceranger_adata,target_sum=1e4, inplace=True)
spaceranger_adata.layers['cp10k'] = spaceranger_adata.X
spaceranger_adata.layers['log1p_cp10k'] =np.log(1+pd.DataFrame.sparse.from_spmatrix(spaceranger_adata.layers['cp10k']))

ins_exp_spaceranger_norm = spaceranger_adata.layers['log1p_cp10k'].T[list(spaceranger_adata.var.index).index(gene)]
colors_spaceranger_norm = pd.DataFrame.sparse.from_spmatrix(ins_exp_spaceranger_norm)
df = pd.DataFrame([x,y,list(colors_spaceranger_norm.loc[0])]).T
df.to_csv(outDir+'spot_'+gene+'_norm_expression.csv', index=False, header=['x','y','expr'])

#SpotClean-corrected data
for r in radii:
    spotclean_adata = sc.read_h5ad(spotcleanDir+'radius='+str(r)+'/cleaned_feature_bc_matrix.h5ad')
    spotclean_adata.layers['scnormalized'] = spotclean_adata.X
    spotclean_adata.layers['log1p_scnormalized'] =np.log(1+pd.DataFrame.sparse.from_spmatrix(spotclean_adata.layers['scnormalized']))

    ins_exp_spotclean = spotclean_adata.layers['log1p_scnormalized'].T[list(spotclean_adata.var.index).index(gene)]
    colors_spotclean = pd.DataFrame.sparse.from_spmatrix(ins_exp_spotclean)
    df = pd.DataFrame([x,y,list(colors_spotclean.loc[0])]).T
    os.makedirs(outDir+'radius='+str(r), exist_ok=True)
    df.to_csv(outDir+'radius='+str(r)+'/spot_'+gene+'_expression.csv', index=False, header=['x','y','expr'])

