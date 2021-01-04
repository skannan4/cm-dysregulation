#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 15:32:15 2020

@author: skannan4
"""

#Trajectory Reconstruction of PSC-CM vs Endogenous CMs
#Chulan Kwon Laboratory
#Primary author: Suraj Kannan
#January 04, 2021
#Document: Velocity Analysis 1.1

###Imports - note, all of these aren't necessary, I borrowed many off of the scanpy vignettes as just in case imports
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
from gprofiler import gprofiler
import scvelo as scv
import rpy2.rinterface_lib.callbacks
import logging
from rpy2.robjects import pandas2ri
import anndata2ri
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)
# Automatically convert rpy2 outputs to pandas dataframes
pandas2ri.activate()
anndata2ri.activate()
plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 3
#sc.set_figure_params(dpi=200s, dpi_save=300)
sc.logging.print_versions()
scv.set_figure_params('scvelo') 

###Set the working directory
os.chdir("/home/skannan4/Documents/Research/BigRef/FinalWorkingFiles/")

###Generate Main adata object
adata = sc.read("python_spliced.csv").transpose()
adata = adata[:,1:29379]

#Obs and vars
phenotype = pd.read_csv("python_pheno.csv", index_col = 0)
genes = pd.read_csv("python_spliced.csv").iloc[:, [0]]
genes.set_index(genes.columns[0], inplace = True)
adata.obs = phenotype
adata.var = genes
adata.obs['timepoint'] = adata.obs['timepoint'].astype('category')
adata.obs['timepoint'].cat.reorder_categories(['e14', 'e18', 'p0', 'p1', 'p4', 'p8', 'p11', 'p14', 'p15', 'p18', 'p22', 'p28', 'p35', 'p56', 'p84', 'D8', 'D10', 'D12', 'D15', 'D18', 'D25', 'D30', 'D45'], inplace=True)

#Spliced and unspliced matrices
unspliced = sc.read("python_unspliced.csv").transpose()[:,1:29379].X
spliced = sc.read("python_spliced.csv").transpose()[:,1:29379].X
adata.layers["spliced"] = spliced
adata.layers["unspliced"] = unspliced

#Add on mnnCorrected PCA plot and UMAP (this is from the cds_branched object in the R code, which contains both the in vivo and in vitro samples)
adata.obsm["X_umap"] = np.genfromtxt("python_umap.csv", delimiter = ",")
adata.obsm["X_pca"] = np.genfromtxt("python_pca.csv", delimiter = ",")

###Analysis for in vivo data

#Construct the adata_invivo object
adata_invivo = adata[adata.obs["group"] == "in vivo"]
adata_invivo.obsm["X_umap"] = np.genfromtxt("python_umap_invivo.csv", delimiter = ",")
adata_invivo.obsm["X_pca"] = np.genfromtxt("python_pca_invivo.csv", delimiter = ",")
#In the default scvelo vignette, you filter genes based on the highest dispersion genes. In my case, however, logically I thought it made more sense to do this based on the genes that are differentially expressed over pseudotime because a) I have this data already and b) it makes more biological sense. I am specifically interested in dynamics associated with genes changing over the maturation trajectory; you could imagine genes with high dispersion but related to some other aspect of biology (or a batch effect) that confounds our results here. Hence, I subset based on differentially expressed genes.
pseudo_ids_invivo = pd.read_table("pseudo_ids_invivo.txt", sep = ",")
adata_invivo = adata_invivo[:,adata_invivo.var.index.isin(pseudo_ids_invivo["Genes"])]

#RNA Velocity analysis
scv.pl.proportions(adata_invivo, groupby = "timepoint", dpi = 200) #Note that the percentage of unspliced decreases steadily over time, as we would expect; by the way, we save this file as Supp 1d
scv.pp.filter_and_normalize(adata_invivo, min_shared_counts=20) #Note - we don't need to identify a gene list here because we have already chosen only the differentially expressed genes
scv.pp.moments(adata_invivo) #Note - we also don't need to compute a PCA because we have already loaded one in. This approach is seemingly better because the PCA we chose is actually the mnnCorrected "Aligned" PCA - thus the neighbours can more correctly be identified without worrying about batch effects
scv.tl.recover_dynamics(adata_invivo, var_names = "all") #This takes a while, but the purpose of recovering all the dynamics is if we want to analyze individual velocities later. I am aware that the default velocity approach doesn't use all the genes, but that's fine.
scv.tl.velocity(adata_invivo, mode='dynamical')
scv.tl.velocity_graph(adata_invivo)
scv.pl.velocity_embedding_stream(adata_invivo, basis='umap', color = "timepoint")
#TEMPORARY APPROACH: So the velocity curves generated are *exactly* backwards. I need to better understand why this happens, and have reached out to the scvelo authors on this. As a temporary measure, I'm simply reversing the direction of the UMAP plot, and replotting. This is TEMPORARY UNTIL FURTHER NOTICE.
adata_invivo.obsm["velocity_umap"] = adata_invivo.obsm["velocity_umap"] * -1
scv.pl.velocity_embedding_stream(adata_invivo, basis='umap', color = "timepoint", legend_loc = "right margin", xlabel = "Component 1", ylabel = "Component 2", palette = ["#F8766D", "#E58700", "#C99800", "#A3A500", "#6BB100", "#00BA38", "#00BF7D", "#00C0AF", "#00BCD8", "#00B0F6" , "#619CFF", "#B983FF", "#E76BF3", "#FD61D1", "#FF67A4"], frameon = True, figsize = (5,5), dpi = 200, linewidth = 0.5, save = "invivo_velocity.svg") #Fig 1e
scv.pl.velocity_embedding_stream(adata_invivo, basis='umap', color = "timepoint")
scv.tl.velocity_confidence(adata_invivo)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata_invivo, c=keys, cmap='coolwarm', perc=[5, 95])
scv.tl.velocity_pseudotime(adata_invivo, root_key = 288, end_key = 903)
scv.tl.latent_time(adata_invivo)
top_genes_invivo = adata_invivo.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.scatter(adata_invivo, basis=top_genes_invivo[:20], ncols=5, frameon=False, color = "timepoint", palette = ["#F8766D", "#E58700", "#C99800", "#A3A500", "#6BB100", "#00BA38", "#00BF7D", "#00C0AF", "#00BCD8", "#00B0F6" , "#619CFF", "#B983FF", "#E76BF3", "#FD61D1", "#FF67A4"])

#Writing files necessary for analysis in R
#In particular - I want the updated obs since it has the velocity lengths
adata_invivo.obs.to_csv("velocity_results_invivo.csv")

###Analysis for in vitro data

#Construct the adata_invitro object
adata_invitro = adata[adata.obs["group"] == "in vitro"]
adata_invitro.obsm["X_umap"] = np.genfromtxt("python_umap_invitro.csv", delimiter = ",")
adata_invitro.obsm["X_pca"] = np.genfromtxt("python_pca_invitro.csv", delimiter = ",")
#As above, I focused on differentially expressed genes.
pseudo_ids_invitro = pd.read_table("pseudo_ids_invitro.txt", sep = ",")
adata_invitro = adata_invitro[:,adata_invitro.var.index.isin(pseudo_ids_invitro["Genes"])]

#RNA Velocity analysis
scv.pl.proportions(adata_invitro, groupby = "timepoint", dpi = 200)
scv.pp.filter_and_normalize(adata_invitro, min_shared_counts=20) #As above - differentially expressed genes
scv.pp.moments(adata_invitro) #As above - using Aligned PCA
scv.tl.recover_dynamics(adata_invitro, var_names = "all")
scv.tl.velocity(adata_invitro, mode='dynamical')
scv.tl.velocity_graph(adata_invitro)
scv.pl.velocity_embedding_stream(adata_invitro, basis='umap', color = "timepoint")
#Unlike the in vivo case, by default the in vitro arrays appear to be going in the correct direction. I need to understand why this is a bit better. Note - the arrow directions are inevitably going to look really messy here. Part of this is due to the incomplete and weak nature of in vitro differentiation, part of it is due to the "broad" trajectory, which is something of an unfortunate consequence from UMAP. A Monocle 2 trajectory sure would have been nice...
scv.pl.velocity_embedding_stream(adata_invitro, basis='umap', color = "timepoint", legend_loc = "right margin", xlabel = "Component 1", ylabel = "Component 2", palette = ["#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC"], frameon = True, figsize = (5,5), dpi = 200, linewidth = 0.5, save = "invitro_velocity.svg") #Fig 1e
scv.tl.velocity_confidence(adata_invitro)
scv.pl.scatter(adata_invitro, c=keys, cmap='coolwarm', perc=[5, 95])
scv.tl.velocity_pseudotime(adata_invitro)
scv.tl.latent_time(adata_invitro)
top_genes_invitro = adata_invitro.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.scatter(adata_invitro, basis=top_genes_invitro[:20], ncols=5, frameon=False, color = "timepoint", palette = ["#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC"])

#Writing files necessary for analysis in R
#In particular - I want the updated obs since it has the velocity lengths
adata_invitro.obs.to_csv("velocity_results_invitro.csv")
