# Trajectory reconstruction identifies dysregulation of perinatal maturation programs in pluripotent stem cell-derived cardiomyocytes

### Introduction
We performed an scRNA-seq study to reconstruct the trajectories of cardiomyocyte (CM) maturation in both endogenous and pluripotent stem cell (PSC)-derived CMs. The goal was to develop a better understanding of the pathways underlying endogenous CM maturation, and investigate how these pathways are dysregulated during in vitro differentiation. Our results can be found in our [preprint](https://www.biorxiv.org/content/10.1101/2021.01.31.428969v1). Here, we share all of the materials needed to reproduce our analysis.

### Data
All of the relevant files to reproduce the analysis can be downloaded at our [Synapse](https://www.synapse.org/#!Synapse:syn23667436/files/). In particular, the following files may be of relevance:

- `gene_mult.mtx`, `gene_mult.genes.txt`, `gene_mult.barcodes.txt`: These are the direct mapped outputs of kallisto|bustools, for those who want to start from the raw mapped counts and pursue their own pipelines. Please note, however, that because of some accidental errors in the sequencing preparation, some barcodes had unexpected collisions and therefore needed to be removed, which can be done by following our pipeline.

- `phenotype_good.txt`: Metadata for the cells in the study. In particular, this includes metadata for only the good barcodes (e.g. not those that had collisions), so those working from `gene_mult.mtx` should use this to subset good barcodes.

- `data_workspace.RData`: A pre-made workspace that contains a lot of the main functions as well as likely objects of interest for users. For the sake of space, this doesn't contain *all* of the objects of interest. For example, the counts tables for the human PSC-CM and perturbation datasets have been removed. However, those can be found in addition workspaces provided on Synapse, such as `human_psccm_data.RData` and `te_data.Rdata`. We recommend `data_workspace.Rdata` for those who quickly want to look for specific findings in our data without having to rerun the entire pipeline.

The following relevant files can be found on Github:


- `bigref_code.R`: An R file containing the entire analysis workflow for the manuscript. All of the files required to reproduce this workflow should be on the Synapse, and a significant chunk of this workflow has been run for `data_workspace.RData`.

- `bigref_figures.R`: A R file containing code necessary to reproduce all figures in the manuscript. Please note that this doesn't run independently of `bigref_code.R`; it was mostly designed just to isolate the plotting functions for figures into one clear place.

- `bigref.py`: A python file containing code necessary to reproduce the RNA Velocity part of the workflow. All of the outputs of this section are also available in the Synapse folder for direct import into R.

### Dependencies
Most of the libraries used in our codebase can be found from CRAN or Bioconductor. However, we additionally make use of the SingleCellNet package from the Cahan lab. Please see [their github](https://github.com/pcahan1/singleCellNet) for instructions on how to install SingleCellNet. Additionally, we make use of scvelo from the Theis lab. Please see [their github](https://github.com/theislab/scvelo) and [documentation](https://scvelo.readthedocs.io/) for more information on installing and running scvelo.

### How to replicate our workflow
To replicate the workflow from the very top, please follow these steps:

1. Download `bigref_code.R`, `bigref_figures.R`, `bigref.py`, and all of the files in the listed Synapse folder.

2. Modify the appropriate lines of both `bigref_code.R` and `bigref_figures.R`: `setwd("~/Documents/Research/BigRef/FinalWorkingFiles/")` to set the working directory to the same working directory you downloaded the Synapse files into it.

3. Similarly, modify the appropriate line of `bigref.py`: `os.chdir("/home/skannan4/Documents/Research/BigRef/FinalWorkingFiles/")` to set the appropriate working directory.

4. Run code for any workflow section or figure of interest. The figures are broken into self-contained code-blocks; at the end of each code block, all temporary objects created for that figure are removed. Each code block should also have details on the figure itself that may be helpful.

Alternatively, note that you can start by loading `data_workspace.RData` into R or an IDE of your choice. As mentioned above, this *should* contain most of the objects of interest for preliminary re-analysis, and can of course be supplemented simply by executing any desired sections in `bigref_code.R`.

Please feel free to email or raise an issue if any of the code doesn't work as claimed!
