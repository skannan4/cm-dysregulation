#Trajectory Reconstruction of PSC-CM vs Endogenous CMs
#Chulan Kwon Laboratory
#Primary author: Suraj Kannan
#December 29, 2020
#Document: Exploratory Code 4.4

#####Load in all necessary libraries
library(ggplot2)
library(Matrix)
library(stringr)
library(monocle3)
library(pheatmap)
library(singleCellNet)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(clipr)
library(DESeq2)

#####Set the working directory
#Note to users - modify this line and download all of the appropriate files from Synapse
setwd("~/Documents/Research/BigRef/FinalWorkingFiles/")

#####Functions to be used

#Note - earlier versions of this code had several helper functions for plotting/calculating first and second derivatives of gene expression. Since the newer version of this code relies on a different approach for gene dynamics, most of those helper functions are removed.

#Function to quickly make a Monocle object from the human data, containing reduced dimension information as well as entropy placed in the proper (pseudotime) slot. Note that this requires some information from the entropy workspace, loaded in below.
makeMonocle = function(data){
  good_cells = combined_datasets[combined_datasets$data == data & combined_datasets$good_cell == TRUE & combined_datasets$include_dataset == TRUE & !combined_datasets$timepoint %in% c("D0", "D2", "D5"), ]$cellname
  tab = get(data)[, good_cells]
  meta.data = combined_datasets[combined_datasets$data == data & combined_datasets$good_cell == TRUE & combined_datasets$include_dataset == TRUE & !combined_datasets$timepoint %in% c("D0", "D2", "D5"), ]
  colnames(tab) = paste(good_cells, data, sep = "_")
  rownames(meta.data) = paste(good_cells, data, sep = "_")
  cds = new_cell_data_set(tab, cell_metadata = meta.data, gene_metadata = data.frame(row.names = rownames(tab), gene1 = rownames(tab), gene2 = rownames(tab)))
  fData(cds)$num_cells_expressed = rowSums(exprs(cds) >= 1)
  cds@colData$timepoint_int = as.numeric(cds@colData$timepoint)
  gene_fits <- fit_models(cds[fData(cds)$num_cells_expressed > (0.25 * ncol(cds)), ], model_formula_str = "~timepoint_int")
  fit_coefs <- coefficient_table(gene_fits)
  global_terms <- fit_coefs %>% filter(term == "timepoint_int")
  global_terms = global_terms[order(global_terms$q_value), ]
  global_ids = as.character(global_terms[global_terms$q_value < 0.05, ]$gene1)
  rm(gene_fits, fit_coefs, global_terms)
  fData(cds)$use = (fData(cds)$gene1 %in% global_ids)
  cds <- preprocess_cds(cds, num_dim = 3, use_genes = global_ids)
  cds <- align_cds(cds)
  cds <- reduce_dimension(cds)
  cds@principal_graph_aux$UMAP$pseudotime = -1 * cds@colData$entropy
  names(cds@principal_graph_aux$UMAP$pseudotime) = colnames(cds)
  cds@colData$pseudotime = -1 * cds@colData$entropy
  return(cds)
}

#A quick gene converter so that we can go between orthologous mouse and human genes; similarly, needs some information from the entropy workspace, loaded in below.
convertGenes = function(genes, startingSpecies = "human"){
  if(startingSpecies == "mouse"){
    ans = oTab[match(genes, oTab$mouse), ]$human
    ans = ans[!is.na(ans)]
    return(ans)
  }
  else if(startingSpecies == "human"){
    ans = oTab[match(genes, oTab$human), ]$mouse
    ans = ans[!is.na(ans)]
    return(ans)
  }
  else {
    print("Pick either mouse or human.")
  }
}

#####Load in working data
#Functions from the entropy project that would be helpful
load("clean_nodatasets_060720.RData")
#Load in the Kallisto|bustools-mapped data
data = readMM("gene_mult.mtx")
rownames(data) = read.table("gene_mult.barcodes.txt", as.is = TRUE)$V1
colnames(data) = read.table("gene_mult.genes.txt", as.is = TRUE)$V1
data = t(data)
rownames(data) = substr(rownames(data), 1, 18)
ercc_data = data[startsWith(rownames(data), "ERCC"), ]
data = rename_genes(data)
#We found that some genes needed to be removed. The first two are well known as genes with significant homology to rRNA. Gm26822 has been identifnied and removed in other studies as well. The rationale for removing Grip2/Zrsr1 is that they are highly expressed only in PSC-CMs only in the kallisto-mapped data (but not in the zUMIs-mapped data, where they are barely expressed at all). This likely indicates some type of mapping artefact, though the nature of the artefact needs to be further clarified.
data = data[!rownames(data) %in% c("Gm42418", "AY036118", "Gm26822", "Grip2", "Zrsr1"), ] 

#Load in pheno_data. We defined this based on our knowledge of the labelings; however, some barcodes had collisions due to mistakes in the library prep. pheno_data removes these. For those working off of the raw FASTQs, please either use this table or contact us for an explanation of which barcodes to discard.
pheno_data = read.table("phenotype_good.txt", as.is = TRUE)
#Cleanup involved with the data, phenotype table, and G_list
data = data[, rownames(pheno_data)]
ercc_data = ercc_data[, rownames(pheno_data)]
timepoint_levels = c("e14", "e18", "p0", "p1", "p4", "p8", "p11", "p14", "p15", "p18", "p22", "p28", "p35", "p56", "p84", "D8", "D10", "D12", "D15" ,"D18", "D25", "D30", "D45")
pheno_data$timepoint = factor(pheno_data$timepoint, levels = timepoint_levels)
pheno_data$group = factor(pheno_data$group, levels = c("in vivo", "in vitro"))
rownames(G_list) = G_list$symbol
#Run our entropy pipeline on the data, which also computes some useful QC metrics (see the entropy manuscript for more details)
data_temp = data_qc(dataset = "data", study = "Kannan Big Ref", timepoint_list = pheno_data$timepoint, scn_calc = TRUE, species = "mouse", sample_type = pheno_data$group, isolation = "LP-FACS", sequencing = "mcSCRB-seq", mapping = "kallisto|bustools", datatype = "UMIs", doi = "None", other_meta = pheno_data$batch)
pheno_data = cbind(pheno_data, data_temp)
pheno_data$ercc = colSums(ercc_data)
colnames(pheno_data)[12] = "timepoint2" #This is because data_qc also adds in a column called "timepoint"
rm(data_temp, ercc_data)

#For all trajectories, I used the cutoffs of top5_norm < 1.8 & depth_norm > -0.7. These are a bit more lax than what we used in the entropy manuscript (top5_norm < 1.3; depth_norm > -0.5), because we could afford to be a bit more relaxed here. I didn't do an aggressive parameter sweep to check these, but from a quick check the output is robust to selection of top5_norm. With regards to depth_norm, the biological interpretation is mostly unchanged in a range between -0.5 to ~1.2, but past ~-0.75, low quality in vivo and in vitro cells get closer to touching one another, likely for technical reasons. I filtered down to -0.7 because that seemed to be a reasonable range for clearing out these technical artefacts without unduly filtering.

#####In vivo CDS

###Construct the cds
cds_invivo <- new_cell_data_set(data[, pheno_data$group == "in vivo" & pheno_data$top5_norm < 1.8 & pheno_data$depth_norm > -0.7], cell_metadata = pheno_data[pheno_data$group == "in vivo" & pheno_data$top5_norm < 1.8 & pheno_data$depth_norm > -0.7, ], gene_metadata = G_list[rownames(data), ])
fData(cds_invivo)$num_cells_expressed = rowSums(exprs(cds_invivo) >= 1) #I turned on multimappers in Kallisto, but I don't feel too bad about considering as "not expressed" if the final count is under 1
cds_invivo <- preprocess_cds(cds_invivo, num_dim = 5)
cds_invivo <- align_cds(cds_invivo, alignment_group = "library")
cds_invivo <- reduce_dimension(cds_invivo, umap.min_dist = 0.2) #Set umap min_dist here just to get a *slightly* neater trajectory (there was an exceedingly mild discontinuity with the default; I don't think setting it makes a huge difference here other than aesthetics)
plot_cells(cds_invivo, label_groups_by_cluster = FALSE, label_cell_groups = FALSE, color_cells_by = "timepoint", cell_size = 1)
cds_invivo <- cluster_cells(cds_invivo, k = 100) #The k is set to aggressively ensure there is only one trajectory - sometimes, I think Monocle gets unduly confused because it *expects* multiple trajectories
cds_invivo <- learn_graph(cds_invivo, learn_graph_control = list(ncenter = 300)) #ncenter set to 300 so that there's more range captured in the pseudotime scores
cds_invivo <- order_cells(cds_invivo) #I selected maybe the first one or two nodes; obviously the values will change a tad bit based on the exact positions, but not enough to make a huge difference
plot_cells(cds_invivo,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5, cell_size = 1)
cds_invivo@colData$pseudotime = cds_invivo@principal_graph_aux$UMAP$pseudotime
ggplot(as.data.frame(cds_invivo@colData), aes(x = timepoint, y = pseudotime, fill = as.factor(batch))) + geom_boxplot() + geom_jitter() + coord_flip()

###Identify differentially regulated genes and plot
gene_fits <- fit_models(cds_invivo[fData(cds_invivo)$num_cells_expressed > 400, ], model_formula_str = "~pseudotime + batch") #Used batch as the term here because I suspected that the batch differences would likely account for more technical variation; also, to ensure that lowly expressed genes didn't dominate up front, I only computed diff gene expression for genes expressed in at least 25% of cells. I was fine with those genes in dimensionality reduction, but not in diff gene expression. However, I felt it would be better to filter before rather than after diff gene computation since it affects q-value computation.
fit_coefs <- coefficient_table(gene_fits)
pseudotime_terms <- fit_coefs %>% filter(term == "pseudotime")
pseudotime_terms = pseudotime_terms[order(pseudotime_terms$q_value), ]
pseudo_ids_invivo = pseudotime_terms[pseudotime_terms$status != "FAIL" & pseudotime_terms$q_value < 0.05, ]$symbol
#Fitting the gene models
model_tbl = fit_models(cds_invivo, model_formula_str = "~ splines::ns(pseudotime, df=3)")
new_data = data.frame(row.names = seq(0, 60), pseudotime = seq(0, 60))
model_expectation = model_predictions(model_tbl, new_data = new_data)
model_expectation = as.data.frame(model_expectation)
model_expectation[is.na(model_expectation)] = 0
model_expectation[model_expectation == Inf] = 0
model_expectation_invivo = model_expectation
#Plot heatmap
pheatmap(model_expectation_invivo[pseudo_ids_invivo, ], show_rownames = FALSE, show_colnames = TRUE, cluster_rows = TRUE, cluster_cols = FALSE, scale = "row", color =  colorRampPalette(brewer.pal(9, "Reds"))(100))
rm(gene_fits, fit_coefs, pseudotime_terms, model_tbl, new_data, model_expectation)

###Identify timing of gene differential expression
#The idea behind this analysis was to understand the timing of differential gene changes. The initial idea was to understand what percentage of genes that are differentially regulated across the entirety of maturation are differentially regulated across the embryonic/neonatal development (pseudotime [0, 20])? Or perinatal development (pseudotime [0, 40])? Of course, this can be done in a continuous way. What we do here is recompute differential gene testing (subsetting only on the in vivo differential gene list identified above) for subsets of pseudotime, e.g. [0, 1], [0, 2], [0, 3], ... , [0, 60]. This can later be used to find out what percentage of the toatl differentially expressed genes are differentially expressed by some pseudotime cutoff.
subset_diff = list()
for(sub in seq(1, 60)){
  cds = cds_invivo[, cds_invivo@colData$pseudotime <= sub]
  gene_fits <- fit_models(cds[pseudo_ids_invivo, ], model_formula_str = "~pseudotime + batch")
  fit_coefs <- coefficient_table(gene_fits)
  pseudotime_terms <- fit_coefs %>% filter(term == "pseudotime")
  pseudotime_terms = pseudotime_terms[order(pseudotime_terms$q_value), ]
  subset_diff[[as.character(sub)]] = pseudotime_terms[pseudotime_terms$status != "FAIL" & pseudotime_terms$q_value < 0.05, ]$symbol
  rm(gene_fits, fit_coefs, pseudotime_terms, cds)
  gc()
  print(sub)
}

###Identify dynamics of individual genes
#The section above enables us to look at the broad timing of differential gene expression. Another idea is to look at how individual genes undergo their activity. The previous version of this code used a rather complicated approach based on derivatives to identify plateau points, but this was overly complicated. Instead, the current approach is to quantify the (pseudo)time required to achieve some percentage of the fold change of the gene. This approach is easy-to-understand, easy-to-compute, and also flexible in that the desired percentage cutoff can be readily changed to yield different pieces of information. For example, a 10% cutoff (e.g. time to 10% FC) tells us, in a sense, how long it takes for a gene to start its activity, while the 95% cutoff (e.g. time to 95% FC) tells us a final plateau time. To define the FC, we consider the maximum FC undergone by a gene in [0, 42] (thereby focusing on the perinatal period). 42 was selected as the upper cutoff because, like all good millenials, I love me a nice pop culture reference.

pseudo_ids_invivo_perinatal = unique(unlist(subset_diff[1:42])) #Getting all genes that are differentially expressed in some subset through the perinatal period"
median_curves = pheatmap(model_expectation_invivo[pseudo_ids_invivo_perinatal, 0:43], show_rownames = FALSE, show_colnames = TRUE, cluster_rows = TRUE, cluster_cols = FALSE, scale = "row", color =  colorRampPalette(brewer.pal(9, "Reds"))(100), cutree_rows = 2) #just a temporary save - we're using this to quickly get which genes are upregulated and which ones are downregulated; I've subset from [0, 42] here to deal with genes like Ttn which otherwise spuriously get marked as downregulated.
test_cut = cutree(median_curves$tree_row, k = 2)
invivo_genes = data.frame(row.names = pseudo_ids_invivo_perinatal)
invivo_genes$direction = "up"
invivo_genes$filler = ""
invivo_genes[names(test_cut)[test_cut == 2], ]$direction = "down"
rm(median_curves, test_cut)

#Compute the time to 10%, 50%, and 95% fold change
model_expectation_invivo_fc = model_expectation_invivo[pseudo_ids_invivo, 1:43] #Set to pseudotime <= 42
model_expectation_invivo_fc = sweep(model_expectation_invivo_fc, 1, model_expectation_invivo_fc[, 1], "/") #This sets the model changes to being measured by fold change over the starting point
invivo_genes$tt10 = 0
invivo_genes$tt50 = 0
invivo_genes$tt95 = 0
for(gene in rownames(invivo_genes)){
  #Need to use separate logic for up and downregulated genes; but same premise
  if(invivo_genes[gene, ]$direction == "up"){ 
    invivo_genes[gene, ]$tt10 = which(model_expectation_invivo_fc[gene, ] >= (1+((max(model_expectation_invivo_fc[gene, ]) - 1) * 0.10)))[1] - 1
    invivo_genes[gene, ]$tt50 = which(model_expectation_invivo_fc[gene, ] >= (1+((max(model_expectation_invivo_fc[gene, ]) - 1) * 0.50)))[1] - 1
    invivo_genes[gene, ]$tt95 = which(model_expectation_invivo_fc[gene, ] >= (1+((max(model_expectation_invivo_fc[gene, ]) - 1) * 0.95)))[1] - 1
  }
  if(invivo_genes[gene, ]$direction == "down"){
    invivo_genes[gene, ]$tt10 = which(model_expectation_invivo_fc[gene, ] <= (1-((1-min(model_expectation_invivo_fc[gene, ])) * 0.10)))[1] - 1
    invivo_genes[gene, ]$tt50 = which(model_expectation_invivo_fc[gene, ] <= (1-((1-min(model_expectation_invivo_fc[gene, ])) * 0.50)))[1] - 1
    invivo_genes[gene, ]$tt95 = which(model_expectation_invivo_fc[gene, ] <= (1-((1-min(model_expectation_invivo_fc[gene, ])) * 0.95)))[1] - 1
  }
}
invivo_genes$filler = NULL
#Clustering invivo_genes using kmeans on the calculated parameters: I determined the number of clusters in a bit of a manual way. I added clusters so long as a population with new dynamics emerged, and I stopped adding clusters when a new cluster had identical parameters to an already existing cluster. It ended up being that 5 clusters was appropriate.
set.seed(11)
invivo_genes$cluster = kmeans(invivo_genes[, c("tt10", "tt50", "tt95")], 5, nstart = 10)$cluster
remap_clusters = c(4, 3, 2, 5, 1)
names(remap_clusters) = seq(1:5)
invivo_genes$cluster = remap_clusters[invivo_genes$cluster]
cluster_names = c("EA1", "EA2", "EA3", "LA1", "LA2")
invivo_genes$cluster_name = cluster_names[invivo_genes$cluster]
#Plot kinetic heatmaps
m_melt = melt(t(scale(t(model_expectation_invivo[pseudo_ids_invivo_perinatal, 1:43]))))
m_melt$cluster_name = invivo_genes[m_melt$Var1, ]$cluster_name
m_melt$direction = invivo_genes[m_melt$Var1, ]$direction
ggplot(m_melt, aes(x = Var2, y = value)) + geom_density_2d() + facet_wrap(direction~cluster_name, scale = "free_y") + labs(title = "Kinetic maps for each gene cluster") + theme(plot.title = element_text(hjust = 0.5))
#Plot parameters
invivo_genes_melt = melt(as.matrix(invivo_genes[, !colnames(invivo_genes) %in% c("direction", "cluster_name", "cluster")]))
invivo_genes_melt$direction = invivo_genes[invivo_genes_melt$Var1, ]$direction
invivo_genes_melt$cluster_name = invivo_genes[invivo_genes_melt$Var1, ]$cluster_name
supp.labs = c("Time to 10% FC", "Time to 50% FC", "Time to 95% FC")
names(supp.labs) = c("tt10", "tt50", "tt95")
invivo_genes_melt$prettylabel = supp.labs[invivo_genes_melt$Var2] 
ggplot(invivo_genes_melt, aes(x = as.factor(cluster_name), y = value, fill = direction)) + geom_boxplot() + facet_wrap(~prettylabel, scale = 'free_y', labeller = labeller(sample_type = supp.labs))
#Plot smoothed curves
median_curves = with(m_melt, tapply(value, list(Var2, paste(direction, cluster_name)), median))
median_curves = sweep(median_curves, 2, apply(median_curves, 2, min))
median_curves = sweep(median_curves, 2, apply(median_curves, 2, max), "/")
median_curves = melt(median_curves)
rm(supp.labs, model_expectation_invivo_fc, cluster_names, remap_clusters)

###Identify TF regulators of each cluster
#We use overrepresentation analysis through WebGestalt to identify TFs whose target genes are enriched in each cluster. WebGestalt largely draws from the MSigDB from GSEA for TF analysis. Since cluster 5 (the late-activating switch genes) is so small, we combine with cluster 4 (early-activating switch genes) for analysis. Because of the nature of MSigDB's annotations, this section required some tweaking. Firstly, we download the list of all enriched TFs based on FDR < 0.05. We then merge the lists together for the four gene lists. Next, we manually annotate the TFs, combining redundant TFs or large families. Then, we combine the redundant TFs by selecting the max enrichment score for each TF for each list. The results can then be plotted.

clus1 = read.table("TF_new/cluster1/results.txt", as.is = TRUE, header = TRUE, row.names = 1)
clus2 = read.table("TF_new/cluster2/results.txt", as.is = TRUE, header = TRUE, row.names = 1)
clus3 = read.table("TF_new/cluster3/results.txt", as.is = TRUE, header = TRUE, row.names = 1)
clus45 = read.table("TF_new/cluster45/results.txt", as.is = TRUE, header = TRUE, row.names = 1)
#Merge together the results
clus = merge(data.frame(row.names = rownames(clus1), clus1 = clus1$enrichmentRatio), data.frame(row.names = rownames(clus2), clus2 = clus2$enrichmentRatio), by = "row.names", all = TRUE)
rownames(clus) = clus$Row.names
clus$Row.names = NULL
clus = merge(clus, data.frame(row.names = rownames(clus3), clus3 = clus3$enrichmentRatio), by = "row.names", all = TRUE)
rownames(clus) = clus$Row.names
clus$Row.names = NULL
clus = merge(clus, data.frame(row.names = rownames(clus45), clus45 = clus45$enrichmentRatio), by = "row.names", all = TRUE)
rownames(clus) = clus$Row.names
clus$Row.names = NULL
#Yup, this is how I annotated them. Wish there were a better way.
clus$tf = c("UNKNOWN", "NFMUE1", "SRF", "E2F1", "UNKNOWN", "ARNT", "YY1", "UNKNOWN", "CMYB", "UNKNOWN", "UNKNOWN", "SRF", "NGFIC", "UNKNOWN", "GABP", "ELK1", "UNKNOWN", "NRF1", "YY1", "UNKNOWN", "SREBP1", "ERR1", "SRF", "SREBP1", "ZF5", "SRF", "MAZ", "SOX9", "SOX9", "UNKNOWN", "YY1", "UNKNOWN", "EVI1", "UNKNOWN", "SF1", "GABP", "UNKNOWN", "AP2ALPHA", "E2F", "HSF", "NRF2", "DR3", "E2F", "USF", "E2F", "EF21DP1RB", "SPZ1", "RREB1", "XBP1", "ETS", "ATF3", "E2F", "AP1", "NRF1", "E2F1", "UNKNOWN", "TEF1", "CP2", "ETS2", "AP2REP", "AP2GAMMA", "SMAD", "E2F", "UNKNOWN", "EGR1", "E2F", "TEL2", "AP1FJ", "HSF1", "UNKNOWN", "SREBP1", "AP1", "AP2", "E4F1", "YY1", "CREB", "SREBP1", "ATF4", "PAX4", "MAX", "MYB", "NFKB", "TCF11", "USF", "E4F1", "UNKNOWN", "UNKNOWN", "SRF", "UNKNOWN", "MYC", "UNKNOWN", "UNKNOWN", "ETS2", "PITX2", "AP1", "UNKNOWN", "UNKNOWN", "TAXCREB", "UNKNOWN", "RSRFC4", "PAX8", "UNKNOWN", "UNKNOWN", "UNKNOWN", "UNKNOWN", "HTF", "UNKNOWN", "PAX4", "RSRFC4", "PAX3", "SMAD4", "GATA2", "FXR", "AHR", "UNKNOWN", "ATF6", "UNKNOWN", "UNKNOWN", "GRE", "ATF1", "E2F", "SP1", "UNKNOWN", "MEF2", "E2F1", "EGR2", "ELF1", "NFKAPPAB65", "UNKNOWN", "UNKNOWN", "UNKNOWN", "UNKNOWN", "CACBINDINGPROTEIN", "E2F1", "VDR", "SRY", "CDC5", "RSRFC4", "PAX6", "IRF", "SP1", "SP1", "CP2", "NERF", "UNKNOWN", "SP1", "E2F", "CEBPGAMMA", "CREL", "USF", "USF", "TFIII", "UNKNOWN", "SRF", "MTF1", "E2F", "MYC", "MYCMAX", "ELK1", "E2F1DP2", "CEBPDELTA", "E2F1DP1", "E2F1DP2", "E2F4DP2", "E2F4DP1", "AP1", "PUI", "MZF1", "NRF2", "NFKB", "ZIC2", "PEA3", "CEBPB", "ZIC3", "TBP", "TEF1", "AML1", "AML1", "HS", "CETS1P54", "OLF1", "MZF1", "CHOP", "PTF1BETA", "PAX4", "ELK1", "MYCMAX", "LFA1", "MYB", "CREB", "STAT1", "HSF2", "ATF", "USF", "HIF1", "FOXO3", "AP2", "USF2", "TCF11MAFG", "NFAT", "E2F1", "TTF1", "HIF1", "STAT5A", "CREB", "LEF1", "OCT1", "STAT1", "NF1", "FR", "LYF1", "ARNT", "TITF1", "EGR", "TFIIA", "E2F1", "HNF4", "CEBP", "ERR1", "MYCMAX", "AP1", "CACCCBINDINGFACTOR", "MAF", "E4BP4", "UNKNOWN", "GR", "DR1", "FOXO1", "LBP1", "FOXO4", "ETS1", "WHN", "AP4", "UNKNOWN", "SF1", "CIZ", "AP4", "SP1", "CREB", "GATA1", "MYOGENIN", "TATA", "CREB", "PUI", "NMYC", "PITX2", "CEBPB", "UNKNOWN", "UNKNOWN", "UNKNOWN", "MEIS1", "HNF3", "MYOD", "NF1", "MEF2", "NFY", "LEF1", "FREAC2")
#Merging redundant TFs. In some cases, I merged if there was a family of TFs that I didn't mind being treated that way (such as the E2F or CEBP families). In other cases, the changes were specific - e.g. RSRFC4 is actually MEF2A.
clus$tf_family = clus$tf
clus[startsWith(clus$tf, "AP1"), ]$tf_family = "AP1"
clus[startsWith(clus$tf, "AP2"), ]$tf_family = "AP2"
clus[startsWith(clus$tf, "ATF"), ]$tf_family = "ATF"
clus[startsWith(clus$tf, "CEBP"), ]$tf_family = "CEBP"
clus[startsWith(clus$tf, "E2F"), ]$tf_family = "E2F"
clus[startsWith(clus$tf, "ETS"), ]$tf_family = "ETS"
clus[startsWith(clus$tf, "HSF"), ]$tf_family = "HSF"
clus[clus$tf == "NFKAPPAB65", ]$tf_family = "NFKB"
clus[clus$tf == "RSRFC4", ]$tf_family = "MEF2"
clus[startsWith(clus$tf, "TCF11"), ]$tf_family = "TCF11"
clus[startsWith(clus$tf, "USF"), ]$tf_family = "USF"
#Aggregate by TF family, selecting the maximum score
clus_agg = as.data.frame(aggregate(as.matrix(clus[, 1:4]), by=list(clus$tf_family), max))
rownames(clus_agg) = clus_agg$Group.1
clus_agg$Group.1 = NULL
clus_agg = clus_agg[order(rowMeans(clus_agg), decreasing = TRUE), ]
colnames(clus_agg) = c("EA1", "EA2", "EA3", "LA")
#Visualize with heatmap
pheatmap(t(clus_agg[2:26, ]), show_rownames = TRUE, show_colnames = TRUE, cluster_rows = FALSE, cluster_cols = FALSE, color =  colorRampPalette(brewer.pal(9, "Reds"))(100), cellwidth = 10, cellheight = 10)

#####In vitro CDS
#The in vitro data is a bit more complicated. We found at first that the difference across timepoints is small relative to difference due to noise, and thus, Monocle doesn't actually recover the correct timepoint-based trajectory. Instead, we used a semi-supervised approach to get that trajectory. In this section, we first create a cds_invitro, which we will rename cds_invitro_unbiased. We then used the semi-supervised approach to get the cds_invitro object we will use in future sections.

###Construct the initial cds
cds_invitro <- new_cell_data_set(data[ , pheno_data$group == "in vitro" & pheno_data$top5_norm < 1.8 & pheno_data$depth_norm > -0.7], cell_metadata = pheno_data[pheno_data$group == "in vitro" & pheno_data$top5_norm < 1.8 & pheno_data$depth_norm > -0.7, ], gene_metadata = G_list[rownames(data), ])
fData(cds_invitro)$num_cells_expressed = rowSums(exprs(cds_invitro) >= 1)
cds_invitro <- preprocess_cds(cds_invitro, num_dim = 3)
cds_invitro <- align_cds(cds_invitro, alignment_group = "library")
cds_invitro <- reduce_dimension(cds_invitro, umap.min_dist = 0.5) #Again, the setting of the umap min_dist helps with the aesthetics, in particular because there are relatively fewer cells.
plot_cells(cds_invitro, label_groups_by_cluster = FALSE, label_cell_groups = FALSE, color_cells_by = "timepoint", cell_size = 1)
cds_invitro_unbiased = cds_invitro

#Now, it is immediately obvious that this trajectory doesn't correspond to timepoint (which can be further checked by looking at the principal components more carefully). The semi-supervised approach is as follows - we first identify genes that are differentially expressed between D8 and D45 cells. We then re-run the preprocessing step, but limit it to only those differentially expressed genes.
cds_invitro@colData$timepoint_int = as.numeric(cds_invitro@colData$timepoint)-15 #This step isn't necessary, it just helps make things a bit neater
gene_fits <- fit_models(cds_invitro[fData(cds_invitro)$num_cells_expressed > 165, cds_invitro@colData$timepoint %in% c("D8", "D45")], model_formula_str = "~timepoint_int + library") #Only one batch so used library as the confounding term
fit_coefs <- coefficient_table(gene_fits)
timepoint_terms <- fit_coefs %>% filter(term == "timepoint_int")
timepoint_terms = timepoint_terms[order(timepoint_terms$q_value), ]
timepoint_ids_invitro = timepoint_terms[timepoint_terms$status != "FAIL" & timepoint_terms$q_value < 0.05, ]$symbol
rm(gene_fits, fit_coefs, timepoint_terms)

#We can now port the timepoint_ids_invitro list into our dimensionality reduction. There's no need to recreate the whole object here, but it's so small it doesn't hurt. Once we have our dimensionality reduction we can proceed with trajectory reconstruction.
cds_invitro <- new_cell_data_set(data[ , pheno_data$group == "in vitro" & pheno_data$top5_norm < 1.8 & pheno_data$depth_norm > -0.7], cell_metadata = pheno_data[pheno_data$group == "in vitro" & pheno_data$top5_norm < 1.8 & pheno_data$depth_norm > -0.7, ], gene_metadata = G_list[rownames(data), ])
fData(cds_invitro)$num_cells_expressed = rowSums(exprs(cds_invitro) >= 1)
cds_invitro <- preprocess_cds(cds_invitro, num_dim = 3, use_genes = timepoint_ids_invitro) #use_genes to select here
cds_invitro <- align_cds(cds_invitro, alignment_group = "library")
cds_invitro <- reduce_dimension(cds_invitro, umap.min_dist = 1) 
plot_cells(cds_invitro, label_groups_by_cluster = FALSE, label_cell_groups = FALSE, color_cells_by = "timepoint", cell_size = 1)
cds_invitro <- cluster_cells(cds_invitro, k = 100)
cds_invitro <- learn_graph(cds_invitro, learn_graph_control = list(ncenter = 300, minimal_branch_len = 30)) #So - here is one of the challenges of Monocle 3 as opposed to 2, particularly with fewer cells. UMAP does a poor job making a "condensed" trajectory (it looks like a big cloud), which in turn will affect Monocle's ability to identify the graph. There are two trade-offs to balance here - more centers means that the pseudotime is more spread out, capturing a more dynamic range, but also means the chance for more spurrious branches. I found two sets of settings I liked here - ncenter = 100, and ncenter = 300, minimal_branch_length = 30. The latter seems more unnecessarily complicated but I just like it a bit more (though the downstream results in both cases are almost identical). If someone can give me a more compelling set of parameters, I'd use them here.
cds_invitro <- order_cells(cds_invitro)
plot_cells(cds_invitro,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5, cell_size = 1)
cds_invitro@colData$pseudotime = cds_invitro@principal_graph_aux$UMAP$pseudotime
ggplot(as.data.frame(cds_invitro@colData), aes(x = timepoint, y = pseudotime)) + geom_boxplot() + geom_jitter() + coord_flip()

###Identify differentially regulated genes and plot
gene_fits <- fit_models(cds_invitro[fData(cds_invitro)$num_cells_expressed > 165, ], model_formula_str = "~pseudotime + library") #Only one batch so used library as the confounding term
fit_coefs <- coefficient_table(gene_fits)
pseudotime_terms <- fit_coefs %>% filter(term == "pseudotime")
pseudotime_terms = pseudotime_terms[order(pseudotime_terms$q_value), ]
pseudo_ids_invitro = pseudotime_terms[pseudotime_terms$status != "FAIL" & pseudotime_terms$q_value < 0.05, ]$symbol
#Fitting the gene models
model_tbl = fit_models(cds_invitro, model_formula_str = "~ splines::ns(pseudotime, df=3)")
new_data = data.frame(row.names = seq(0, 34), pseudotime = seq(0, 34))
model_expectation = model_predictions(model_tbl, new_data = new_data)
model_expectation = as.data.frame(model_expectation)
model_expectation[is.na(model_expectation)] = 0
model_expectation[model_expectation == Inf] = 0
model_expectation_invitro = model_expectation
#Plot heatmap
pheatmap(model_expectation_invitro[pseudo_ids_invitro, ], show_rownames = FALSE, show_colnames = TRUE, cluster_rows = TRUE, cluster_cols = FALSE, scale = "row", color =  colorRampPalette(brewer.pal(9, "Reds"))(100))
rm(gene_fits, fit_coefs, pseudotime_terms, model_tbl, new_data, model_expectation)

#####Combined CDS
#We take several approaches here to combining the in vitro and in vivo trajectories. The first is the most appropriate, "conventional" method. In this approach, we combine the two sets of cells. By setting the libraries from batch 3 (which contain in vivo and in vitro cells) "first" in line, mnnCorrect treats these as the reference batch, thereby correctly aligning subsequent batches. As can be seen, this results in two separate trajectories for the in vitro and in vivo cells.

###Construct combined CDS
cells = c(colnames(cds_invivo), colnames(cds_invitro))
cds_combined = new_cell_data_set(data[ , cells], cell_metadata = pheno_data[cells, ], gene_metadata = G_list[rownames(data), ])
colData(cds_combined)$library = -1*(colData(cds_combined)$library - 37) #This places the libraries in the correct order, such that the libraries with both in vitro and in vivo cells come first
#Set a couple of factors that could be helpful later - the number of cells expressing each gene and their average level in each group
fData(cds_combined)$num_cells_expressed = rowSums(exprs(cds_combined) >= 1)
fData(cds_combined)$num_cells_expressed_invivo = rowSums(exprs(cds_combined[, cds_combined@colData$group == "in vivo"]) >= 1)
fData(cds_combined)$level_invivo = rowMeans(normalized_counts(cds_combined[, cds_combined@colData$group == "in vivo"], norm_method =  "size_only"))
fData(cds_combined)$num_cells_expressed_invitro = rowSums(exprs(cds_combined[, cds_combined@colData$group == "in vitro"]) >= 1)
fData(cds_combined)$level_invitro = rowMeans(normalized_counts(cds_combined[, cds_combined@colData$group == "in vitro"], norm_method =  "size_only"))
#Generating the trajectory
cds_combined = preprocess_cds(cds_combined, num_dim = 15) #After some trial and error, dim = 15 seemed appropriate here
cds_combined = align_cds(cds_combined, alignment_group = "library")
cds_combined = reduce_dimension(cds_combined)
plot_cells(cds_combined, label_groups_by_cluster = FALSE, label_cell_groups = FALSE, color_cells_by = "timepoint", cell_size = 1)
#There isn't necessarily a point in computing pseudotime here, since there are two separate trajectories that more or less completely correspond to the individual trajectories from the previous sections. Thus, we load in those pseudotimes here, though of course they aren't really needed
#individual_pseudotime = scaled (so that in vitro and in vivo both run from 0 to 1) while full_pseudotime = unscaled
colData(cds_combined)$individual_pseudotime = 0
colData(cds_combined)[colnames(cds_invivo)[colnames(cds_invivo) %in% colnames(cds_combined)], ]$individual_pseudotime = colData(cds_invivo)[colnames(cds_invivo)[colnames(cds_invivo) %in% colnames(cds_combined)], ]$pseudotime/max(colData(cds_invivo)[colnames(cds_invivo)[colnames(cds_invivo) %in% colnames(cds_combined)], ]$pseudotime)
colData(cds_combined)[colnames(cds_invitro)[colnames(cds_invitro) %in% colnames(cds_combined)], ]$individual_pseudotime = colData(cds_invitro)[colnames(cds_invitro)[colnames(cds_invitro) %in% colnames(cds_combined)], ]$pseudotime/max(colData(cds_invitro)[colnames(cds_invitro)[colnames(cds_invitro) %in% colnames(cds_combined)], ]$pseudotime)
colData(cds_combined)$full_pseudotime = 0
colData(cds_combined)[colnames(cds_invivo)[colnames(cds_invivo) %in% colnames(cds_combined)], ]$full_pseudotime = colData(cds_invivo)[colnames(cds_invivo)[colnames(cds_invivo) %in% colnames(cds_combined)], ]$pseudotime
colData(cds_combined)[colnames(cds_invitro)[colnames(cds_invitro) %in% colnames(cds_combined)], ]$full_pseudotime = colData(cds_invitro)[colnames(cds_invitro)[colnames(cds_invitro) %in% colnames(cds_combined)], ]$pseudotime
plot_cells(cds_combined, label_groups_by_cluster = FALSE, label_cell_groups = FALSE, color_cells_by = "individual_pseudotime", cell_size = 1)

#The next approach to combining the trajectories makes use of the way mnnCorrect performs batch correction. If we have two celltypes A and B, it makes sense to set the reference batch to one containing both A and B. However, if we set as reference a batch containing *only* A, mnnCorrect will effectively try to project celltype B onto A, which is what we are trying to do here. We thus leave the batches containing only in vivo cells first, and proceed (using otherwise the same settings as above).

###Construct single CDS
cds_single = cds_combined #quick way to not have to rerun all that code
colData(cds_single)$library = -1*(colData(cds_single)$library - 37) #This once again reverses the libraries
cds_single = preprocess_cds(cds_single, num_dim = 15)
cds_single = align_cds(cds_single, alignment_group = "library")
cds_single = reduce_dimension(cds_single)
plot_cells(cds_single, label_groups_by_cluster = FALSE, label_cell_groups = FALSE, color_cells_by = "timepoint", cell_size = 1)
#This time, we compute pseudotimes
cds_single = cluster_cells(cds_single, k = 100)
cds_single = learn_graph(cds_single, learn_graph_control = list(ncenter = 200)) #As with cds_combined, this produces some spurious branches but improves the dynamic range
cds_single = order_cells(cds_single) #The tricky part with these weird trajectories is deciding the root - you can bias your results pretty significantly based on that. With this version, I selected two roots - where the e14 cells are and that lower left tip of cells (mostly in vitro).
plot_cells(cds_single,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5, cell_size = 1)
cds_single@colData$pseudotime = cds_single@principal_graph_aux$UMAP$pseudotime
ggplot(as.data.frame(cds_single@colData), aes(x = timepoint, y = pseudotime, fill = as.factor(batch))) + geom_boxplot() + geom_jitter() + coord_flip() 

#We try another approach as a supplement to further test the forced alignment of the in vivo and in vitro trajectories. In this approach, we treat in vivo vs in vitro group as a batch effect to be corrected by mnnCorrect. Since we'd like both batch and library to be corrected using nonlinear methods (rather than passing one to reduced_model_formula_str), we do a bit of tweaking. First, batch correction is run with library as before. Then, the new aligned coordinates are passed in the PCA slot, and a second round of batch correction is run with group. Since the mnnCorrect alignment draws from the PCA coordinates, it effectively does the second alignment on top of the first one. The end result is functionally similar to cds_single.

###Construct single CDS 2
cds_single2 = cds_combined
cds_single2 = preprocess_cds(cds_single2, num_dim = 15)
cds_single2 = align_cds(cds_single2, alignment_group = "library")
cds_single2@reducedDims$PCA = cds_single2@reducedDims$Aligned #Move the aligned cds into the PCA slot
cds_single2@colData$group_int = 0 #This is a quick way to put the in vivo samples in front
cds_single2@colData[cds_single2@colData$group == "in vitro", ]$group_int = 1
cds_single2 = align_cds(cds_single2, alignment_group = "group_int")
cds_single2 = reduce_dimension(cds_single2, umap.min_dist = 0.4) #Here, without the umap min_dist parameter, the in vitro/neonatal cells separate out from the later stage in vivo cells; this parameter just helps connect the trajectory
plot_cells(cds_single2, label_groups_by_cluster = FALSE, label_cell_groups = FALSE, color_cells_by = "timepoint", cell_size = 1)
#This time, we compute pseudotimes
cds_single2 = cluster_cells(cds_single2, k = 100)
cds_single2 = learn_graph(cds_single2, learn_graph_control = list(ncenter = 200))
cds_single2 = order_cells(cds_single2) #As above, root selection is a bit tricky - I selected the midpoint node where the e14 cells are
plot_cells(cds_single2,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5, cell_size = 1)
cds_single2@colData$pseudotime = cds_single2@principal_graph_aux$UMAP$pseudotime
ggplot(as.data.frame(cds_single2@colData), aes(x = timepoint, y = pseudotime, fill = as.factor(batch))) + geom_boxplot() + geom_jitter() + coord_flip() 

#Now, as a first step of comparison, we look at the fold changes between in vitro and in vivo. There are two particular pseudotimes of interest to us - in vivo pseudotime = 20 (where the in vitro cells appear to be arrested) and in vivo pseudotime = 40 (where the end of main maturation-related trends seems to occur). Hence, we compute the in vivo fold changes at both points to compare to in vitro.
compared_fc = data.frame(row.names = pseudo_ids_invivo_perinatal, invivo_fc_40 = log2((model_expectation_invivo[pseudo_ids_invivo_perinatal, 41] + 0.00001)/(model_expectation_invivo[pseudo_ids_invivo_perinatal, 1] + 0.00001)), invivo_fc_20 = log2((model_expectation_invivo[pseudo_ids_invivo_perinatal, 21] + 0.00001)/(model_expectation_invivo[pseudo_ids_invivo_perinatal, 1] + 0.00001)), invitro_fc = log2((model_expectation_invitro[pseudo_ids_invivo_perinatal, 35] + 0.00001)/(model_expectation_invitro[pseudo_ids_invivo_perinatal, 1] + 0.00001)), diff_invitro = pseudo_ids_invivo_perinatal %in% pseudo_ids_invitro)
#Since we have fold changes, we can now compute genes that are differentially expressed in the correct direction
pseudo_ids_invitro_good = c(rownames(compared_fc)[compared_fc$diff_invitro == TRUE & compared_fc$invitro_fc > 0 & rownames(compared_fc) %in% rownames(invivo_genes)[invivo_genes$direction == "up"]], rownames(compared_fc)[compared_fc$diff_invitro == TRUE & compared_fc$invitro_fc < 0 & rownames(compared_fc) %in% rownames(invivo_genes)[invivo_genes$direction == "down"]])
#There are actually some contexts later in this analysis (e.g. the Venn Diagram plots and IPA analysis) where we actually need the fold changes of all of the differentially expressed genes in vitro, even those that don't overlap with the differentially expressed genes in vivo. We compute these fold changes below.
pseudo_ids_invitro_fc = log2((model_expectation_invitro[pseudo_ids_invitro, 35] + 0.00001)/(model_expectation_invitro[pseudo_ids_invitro, 1] + 0.00001))
names(pseudo_ids_invitro_fc) = pseudo_ids_invitro

#####Mouse Global Gene Dysregulation Analysis
#One observation is that, when using standard batch correction methods, in vivo and in vitro cells separate out completely. This suggests that there are global gene expression differences between these groups of cells. In this section, we identify those genes, looking both at comparison between all in vivo CMs vs in vitro CMs, and just early stage (e.g. in vivo pseudotime < 15) in vivo CMs vs all in vitro CMs. This latter comparison allows us to compare in vitro CMs to their earliest in vivo neighbours.

###All in vivo CM comparison
gene_fits <- fit_models(cds_combined[fData(cds_combined)$num_cells_expressed_invivo > 400 | fData(cds_combined)$num_cells_expressed_invitro > 165, ], model_formula_str = "~group") #Set to 25% of cells expressing either in vitro or in vivo
fit_coefs <- coefficient_table(gene_fits)
global_terms <- fit_coefs %>% filter(term == "groupin vitro")
global_terms = global_terms[order(global_terms$q_value), ]
global_ids = global_terms[global_terms$q_value < 0.05, ]$symbol
global_ids_clean = global_ids[!startsWith(global_ids, "Rps") & !startsWith(global_ids, "Rpl")] #The ribosomal protein coding genes are fine, I just wanted to remove them so I could look at everything else
rm(gene_fits, fit_coefs)

#For the sake of visualization, I wanted to make heatmaps that were binned by pseudotime. In this case, I did so by binning each group by their respective pseudotimes into bins of 5. A couple of observations here. Firstly, in retrospect, I could have done this by using the modeled trend lines, as I did when computing fold change, though the difference is probably negligible. Secondly, for binning, I used the pseudotimes from each individual trajectory, rather than from any combined methods. The reason for this is a) I trust those pseudotimes more and b) they are a bit more spread out, allowing better temporal granularity. But this also means that bin 5 in vivo doesn't necessarily correspond to bin 5 in vitro. I think this is fine - this is for visualizing trends.
normalized_data = normalized_counts(cds_combined) #Get the normalized counts
cds_combined@colData$pseudotime_bin = paste(cds_combined@colData$group, 5*(trunc(cds_combined@colData$full_pseudotime/5)+1))
cds_combined@colData$pseudotime_bin =  factor(cds_combined@colData$pseudotime_bin, levels = c("in vivo 5", "in vivo 10", "in vivo 15", "in vivo 20", "in vivo 25", "in vivo 30", "in vivo 35", "in vivo 40", "in vivo 45", "in vivo 50", "in vivo 55", "in vivo 60", "in vitro 5", "in vitro 10", "in vitro 15", "in vitro 20", "in vitro 25", "in vitro 30", "in vitro 35"))
global_data = t(aggregate(t(as.matrix(normalized_data)), by=list(cds_combined@colData$pseudotime_bin), mean)) #Bins the normalized counts
colnames(global_data) = global_data[1, ]
global_data = global_data[2:nrow(global_data), ]
global_data_numeric = apply(global_data, 2, as.numeric) #The data comes out as character, so this is to convert
rownames(global_data_numeric) = rownames(global_data)
#Plot the globally dysregulated genes
pheatmap(global_data_numeric[global_ids, ], show_rownames = FALSE, show_colnames = TRUE, cluster_rows = TRUE, cluster_cols= FALSE, scale = "row")

###Compared to only in vivo CMs < pseudotime 15
gene_fits <- fit_models(cds_combined[fData(cds_combined)$num_cells_expressed_invivo > 400 | fData(cds_combined)$num_cells_expressed_invitro > 165, (cds_combined@colData$group == "in vivo" & cds_combined@colData$full_pseudotime <= 15) | (cds_combined@colData$group == "in vitro")], model_formula_str = "~group") #Set to 25% of cells expressing either in vitro or in vivo
fit_coefs <- coefficient_table(gene_fits)
global_terms_early <- fit_coefs %>% filter(term == "groupin vitro")
global_terms_early = global_terms_early[order(global_terms_early$q_value), ]
global_ids_early = global_terms_early[global_terms_early$q_value < 0.05, ]$symbol
global_ids_early_clean = global_ids_early[!startsWith(global_ids_early, "Rps") & !startsWith(global_ids_early, "Rpl")]
rm(gene_fits, fit_coefs)

#Looking at the GO terms of globally upregulated and globally downregulated genes suggests, at a glance, that the PSC-CMs are generally immature (cell cycle terms are upregulated, sarcomeric/contraction/metabolic genes are downregulated). But a strange observation is that in many cases, the gene levels don't just fall to the "immature" side of the spectrum - they fall entirely outside the normal spectrum of in vivo development. To visually make this point, we scaled every in vivo differentially regulated gene from 0 to 1, where 0 is the lowest level over maturation and 1 is the highest level. (Note that in some cases, 0 will be the "immature" level and in others the "mature" level, depending on whether the gene is upregulated or downregulated). We apply this same scaling to the mean in vitro gene level as well - and find that many fall outside of the [0, 1] window.
global_vals = data.frame(row.names = pseudo_ids_invivo_perinatal, min_level = rowMins(global_data_numeric[pseudo_ids_invivo_perinatal, 1:9]), max_level = rowMaxs(global_data_numeric[pseudo_ids_invivo_perinatal, 1:9]), invitro_level = as.numeric(rowMeans(normalized_data[pseudo_ids_invivo_perinatal, colnames(cds_invitro)])))
global_vals = sweep(global_vals, 1, global_vals$min_level) #Set the lowest level to 0
global_vals = sweep(global_vals, 1, global_vals$max_level, "/") #Set the highest level to 1
ggplot(global_vals, aes(x = invitro_level)) + geom_point(aes(y = 0), size = 0.5) + geom_line(data = data.frame(x = c(0, 0), y = c(-0.1, 0.1)), aes(x = x, y = y)) + geom_line(data = data.frame(x = c(1, 1), y = c(-0.1, 0.1)), aes(x = x, y = y)) + theme_classic() + geom_density() + xlim(-5, 5)  +  coord_cartesian(clip = 'off') + theme(axis.title = element_text(size = 0))

#####Human dataset analysis
#In this section, we generate results for several human PSC-CM datasets. In particular, we analyzed four datasets - Churko et al., Friedman, Nguyen, and Lukowski et al., Gerbin, Grancharova, Donovan-Maiye, and Hendershott et al., and Ruan and Liao et al. (see the alldata object for more metadata on each of these studies). We were surprised to note that we couldn't really generate trajectories for all of the datasets in Monocle 3 (or even Monocle 2) owing to what I am assuming were batch effects in the data. In specific, unlike our data where a neat trajectory is formed in vivo and in vitro, with these datasets the timepoints continuously separate out into clusters (regardless of whether gene filtering is done, though I'm not showing that here). I assume this happens because each timepoint was prepared separately (e.g. on a separate 10x lane or iCell8 chip). To handle this issue, we instead computed gene expression changes over *entropy* (based on our previous manuscript - see Kannan 2020). Entropy is relatively resistant to batch effects and also captures similar gene expressions to pseudotime, albeit while losing interesting neighbourhood information that may come from dimensionality reduction. At any rate, this works as an approximation to allow us to study the human in vitro datasets as well.

#Load in the datasets. I have prepared a workspace with just the data - it comes from our earlier entropy work.
load("human_psccm_data.RData")

#Generate Monocle objects for each of the human datasets
cds_churko = makeMonocle("churko_data")
cds_friedman = makeMonocle("friedman_data")
cds_gerbin = makeMonocle("gerbin_data")
cds_ruan = makeMonocle("ruan_data")

#For each dataset, identify the differentially expressed genes. As with our analysis, a gene is considered expressed if it is expressed in >25% of cells
human_genes = list()
human_genes_up = list()
human_genes_down = list()
human_expressed = list()
for(cds_name in c("cds_churko", "cds_friedman", "cds_gerbin", "cds_ruan")){
  cds = get(cds_name)
  fData(cds)$num_cells_expressed = rowSums(exprs(cds) >= 1)
  gene_fits <- fit_models(cds[fData(cds)$num_cells_expressed > (0.25 * ncol(cds)), ], model_formula_str = "~entropy")
  fit_coefs <- coefficient_table(gene_fits)
  entropy_terms <- fit_coefs %>% filter(term == "entropy")
  diff = entropy_terms[entropy_terms$q_value < 0.05, ]$gene1
  diff_up = entropy_terms[entropy_terms$q_value < 0.05 & entropy_terms$normalized_effect < 0, ]$gene1 #because entropy proceeds int the opposite direction to pseudotime, the normalized effects will also be reversed for identifying up and downregulated genes
  diff_down = entropy_terms[entropy_terms$q_value < 0.05 & entropy_terms$normalized_effect > 0, ]$gene1 #as above
  human_genes[[cds_name]] = diff
  human_genes_up[[cds_name]] = diff_up
  human_genes_down[[cds_name]] = diff_down
  human_expressed[[cds_name]] = entropy_terms$gene1
  print(cds_name)
  rm(cds, gene_fits, fit_coefs, entropy_terms, diff, diff_up, diff_down)
  gc()
}

#We need to work with orthologous genes. It will be easier to do our comparisons using one species. It could be either mouse or human, of course, but for sake of consistency we use mouse.
human_expressed_mouse = lapply(human_expressed, convertGenes)
human_genes_mouse = lapply(human_genes, convertGenes)
human_genes_up_mouse = lapply(human_genes_up, convertGenes)
human_genes_down_mouse = lapply(human_genes_down, convertGenes)
#Now note - not all of our in vivo genes have orthologs! So we'll subset the pseudo_ids_invivo to only those with orthologs
pseudo_ids_invivo_ortho = pseudo_ids_invivo_perinatal[pseudo_ids_invivo_perinatal %in% oTab$mouse]
#We will need these lists for up and downregulated genes for the dysregulation analysis

#For the sake of plotting and exploring, I also generate a list here of genes that are differentially expressed in vitro and in vivo in the correct direction. A consensus list is made by taking genes that are correctly differentially expressed in at least two studies.
human_correct_mouse = list(c(human_genes_up_mouse[["cds_churko"]][human_genes_up_mouse[["cds_churko"]] %in% rownames(invivo_genes)[invivo_genes$direction == "up"]], human_genes_down_mouse[["cds_churko"]][human_genes_down_mouse[["cds_churko"]] %in% rownames(invivo_genes)[invivo_genes$direction == "down"]]), c(human_genes_up_mouse[["cds_friedman"]][human_genes_up_mouse[["cds_friedman"]] %in% rownames(invivo_genes)[invivo_genes$direction == "up"]], human_genes_down_mouse[["cds_friedman"]][human_genes_down_mouse[["cds_friedman"]] %in% rownames(invivo_genes)[invivo_genes$direction == "down"]]), c(human_genes_up_mouse[["cds_gerbin"]][human_genes_up_mouse[["cds_gerbin"]] %in% rownames(invivo_genes)[invivo_genes$direction == "up"]], human_genes_down_mouse[["cds_gerbin"]][human_genes_down_mouse[["cds_gerbin"]] %in% rownames(invivo_genes)[invivo_genes$direction == "down"]]), c(human_genes_up_mouse[["cds_ruan"]][human_genes_up_mouse[["cds_ruan"]] %in% rownames(invivo_genes)[invivo_genes$direction == "up"]], human_genes_down_mouse[["cds_ruan"]][human_genes_down_mouse[["cds_ruan"]] %in% rownames(invivo_genes)[invivo_genes$direction == "down"]]))
names(human_correct_mouse) = names(human_genes)
human_good = names(table(unlist(human_correct_mouse)))[table(unlist(human_correct_mouse)) >= 2]
human_good_new = names(table(unlist(human_correct_mouse)))[table(unlist(human_correct_mouse)) >= 1]
human_good_new = human_good_new[!human_good_new %in% dysreg]

#####Dysregulated Gene Analysis
#In this section, we identify dysregulated genes in the mouse and human datasets. The logic is as follows. We focused on genes that are differentialy regulated in vivo but not in vitro, as these are evidence of failure to complete the in vivo maturation program. Note - there are some genes that are differentially regulated in vitro but not in vivo, but we didn't pay much attention to them (from my initial glance, a large number are cardiac differentiation genes that were probably already differentially regulated in vivo prior to the start of our analysis - since the in vitro analysis might start a bit before the in vivo trajectories). Thus, we break the in vivo genes into up and downregulated, and check whether these genes are in the corresponding in vitro lists (this also inherently captures genes that are differentially regulated in different directions across both). One major caveat however - we make the assumption that if a gene is downregulated in vivo but not expressed (e.g. expressed < 25% of cells) in vitro, then we can ignore that gene. This is perhaps a debatable criterion - for example, temporal expression followed by downregulation might be required. However, to hone in on the most important genes, we ignore that. (Also note that another major assumption here is that human and mouse should have the same directionality of changes - there are some major caveats to that as well, but based on the available data, it is a good approximation).

#P.S. I'm very aware of how bad these variable names suck. Sorry.

###Mouse dyregulation
pseudo_ids_invitro_down = rownames(compared_fc[rownames(compared_fc) %in% pseudo_ids_invitro & compared_fc$invitro_fc < 0, ])
pseudo_ids_invitro_up = rownames(compared_fc[rownames(compared_fc) %in% pseudo_ids_invitro & compared_fc$invitro_fc > 0, ])
dysreg_up = rownames(invivo_genes)[invivo_genes$direction == "up" & !rownames(invivo_genes) %in% pseudo_ids_invitro_up]
dysreg_down = rownames(invivo_genes)[fData(cds_invitro)[rownames(invivo_genes), ]$num_cells_expressed > 165 & !rownames(invivo_genes) %in% pseudo_ids_invitro_down & invivo_genes$direction == "down"]
dysreg_mouse = c(dysreg_up, dysreg_down)
#We can now categorize genes as being correctly differentially regulated or dysregulated in mouse - we use these categories for some of our plotting.
invivo_genes_melt$mouse = "n/a"
invivo_genes_melt[invivo_genes_melt$Var1 %in% pseudo_ids_invitro_good, ]$mouse = "Correct"
invivo_genes_melt[invivo_genes_melt$Var1 %in% dysreg_mouse, ]$mouse = "Dysregulated"
invivo_genes_melt$mouse = factor(invivo_genes_melt$mouse, levels = c("n/a", "Correct", "Dysregulated"))

###Human dysregulation
#As above, only take ones with orthologous genes
pseudo_ids_invivo_up_ortho = rownames(invivo_genes)[rownames(invivo_genes) %in% pseudo_ids_invivo_ortho & invivo_genes$direction == "up"]
pseudo_ids_invivo_down_ortho = rownames(invivo_genes)[rownames(invivo_genes_down) %in% pseudo_ids_invivo_ortho & invivo_genes$direction == "down"]
human_dysreg_up = lapply(names(human_genes_up), function(x) {pseudo_ids_invivo_up_ortho[!pseudo_ids_invivo_up_ortho %in% human_genes_up_mouse[[x]]]})
names(human_dysreg_up) = names(human_genes_up)
human_dysreg_down = lapply(names(human_genes_down), function(x) {pseudo_ids_invivo_down_ortho[!pseudo_ids_invivo_down_ortho %in% human_genes_down_mouse[[x]] & pseudo_ids_invivo_down_ortho %in% human_expressed_mouse[[x]]]})
names(human_dysreg_down) = names(human_genes_down)
#Our selection rule was that to be selected as a consensus human dysregulated gene, it needed to appear in the dysregulated list for at least 3/4 studies
human_consensus_up = names(table(unlist(human_dysreg_up)))[table(unlist(human_dysreg_up)) >= 3]
human_consensus_down = names(table(unlist(human_dysreg_down)))[table(unlist(human_dysreg_down)) >= 3]
dysreg_human = c(human_consensus_up, human_consensus_down)
invivo_genes_melt$human = "n/a"
invivo_genes_melt[invivo_genes_melt$Var1 %in% human_good, ]$human = "Correct"
invivo_genes_melt[invivo_genes_melt$Var1 %in% dysreg_human, ]$human = "Dysregulated"
invivo_genes_melt$human = factor(invivo_genes_melt$human, levels = c("n/a", "Correct", "Dysregulated"))

#Lastly, the final dysregulated list is the intersection of the human and mouse lists
dysreg = intersect(dysreg_mouse, dysreg_human)

###TF Enrichment Analysis of Dysregulated Genes
#As with the in vivo data, we use overrepresentation analysis through WebGestalt to identify TFs whose targets are enriched in our dysregulated gene list. We handle redundancy a bit differently here - since we are looking for best candidates, we use WebGestalt's affinity propagation (effectively a form of clustering) to select best TFs. My understanding of AP is that it uses a clustering algorithm to group together similar lists. This works nicely, but in some cases the TF that AP selects as the representative of the group isn't the best - either because it is lowly expressed/activated in our tissues, or is of unknown status. In these cases, I manually selected another TF from the group to represent that cluster. This is obviously a tad bit biased, so I try to explain my logic below - but I think my choices are justifiable.

#Load in the results
dysreg_tf_all = read.table("TF_new/dysreg/results.txt", as.is = TRUE, row.names = 1, header = TRUE)
dysreg_tf_ap = read.csv("TF_new/dysreg/ap_results.txt", as.is = TRUE, fill = TRUE, header = FALSE, sep = "\t")
dysreg_tf_ap[6,1] = dysreg_tf_ap[6,5] #NFE2 isn't expressed in vivo or in vitro. However, NRF2/NFE2L2, in the same group, is not only expressed highly in vivo (higher than several of the others in the category), but also has a high IPA activation score. Therefore, it seemed reasonable.
dysreg_tf_ap[7,1] = dysreg_tf_ap[7,4] #The selected TF is "unknown." Of the remaining TFs in the cluster, Mef2a is probably the most well-known and has previous evidence suggesting impact on CM maturation.
dysreg_tf_ap[8, 1] = dysreg_tf_ap[8, 2] #ER isn't expressed in vivo, but Err1 is, and is also supported by a reasonable amount of literature.
dysreg_tf_good = dysreg_tf_all[dysreg_tf_ap$V1, ]
dysreg_tf_good$clean_name = c("NRF1", "SRF", "YY1", "SOX9", "AP1FJ", "NRF2", "MEF2", "ERR1", "PPAR")
dysreg_tf_good = dysreg_tf_good[order(dysreg_tf_good$enrichmentRatio, decreasing = TRUE), ]

###IPA Analysis
#Because of our concerns of global gene expression differences, we also looked at IPA predicted activity for these TFs. To get the IPA results, we need to write fold changes to text files that can be uploaded into IPA. There are several comparisons of interest, written below.
write.table(data.frame(genes = rownames(compared_fc), fc = compared_fc$invivo_fc_20), "IPA_invivo_0_20.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(data.frame(genes = rownames(compared_fc), fc = compared_fc$invivo_fc_40), "IPA_invivo_0_40.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(data.frame(genes = rownames(compared_fc), fc = log2((model_expectation_invivo[pseudo_ids_invivo_perinatal, 41] + 0.00001)/(model_expectation_invivo[pseudo_ids_invivo_perinatal, 21] + 0.00001))), "IPA_invivo_20_40.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(data.frame(genes = rownames(compared_fc), fc = compared_fc$invitro_fc), "IPA_invitro_0_40.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(data.frame(genes = pseudo_ids_invitro, fc = pseudo_ids_invitro_fc), "IPA_invitro_diff.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

#Load in IPA results
dysreg_ipa = read.table("TF_new/dysreg/dysreg_ipa_results.txt", as.is = TRUE, row.names = 1, header = TRUE, sep = "\t")
dysreg_ipa = as.data.frame(data.matrix(dysreg_ipa))

###ATAC-seq Analysis
#The purpose of this section is to assess the chromatin accessibility of dysregulated vs correctly regulated genes in PSC-CMs. We drew from three datasets of PSC-CMs in the literature, spanning from D15-D30.

#For the purpose of this analysis, we want to remove mitochondrially encoded genes, as I'm pretty sure ATAC-seq (which uses nuclei) wouldn't be able to detect these guys at all.
human_good_clean = human_good[!startsWith(human_good), "mt-"]

#Load in the HOMER annotated peaks. Obviously, these file locations will need to be modified to the user's setup.
#Study: Liu et al. (DOI: 10.1161/CIRCRESAHA.116.310456)
#Timepoint: D30
#Data: GSE85330
c15_1 = read.csv("psc-atac/c15_1.txt", sep = "\t")
c15_2 = read.csv("psc-atac/c15_2.txt", sep = "\t")
c20_1 = read.csv("psc-atac/c20_1.txt", sep = "\t")
c20_2 = read.csv("psc-atac/c20_2.txt", sep = "\t")
h1_1 = read.csv("~psc-atac/h1_1.txt", sep = "\t")
h1_2 = read.csv("psc-atac/h1_2.txt", sep = "\t")
h9_1 = read.csv("psc-atac/h9_1.txt", sep = "\t")
h9_2 = read.csv("psc-atac/h9_2.txt", sep = "\t")
#Study: Bertero and Fields et al. (DOI: 10.1038/s41467-019-09483-5)
#Timepoint: D14
#Data: GSE106689
rues2_1 = read.csv("psc-atac/rues2_d14_1.txt", sep = "\t")
rues2_2 = read.csv("psc-atac/rues2_d14_2.txt", sep = "\t")
#Study: Greenwald and Li et al. (DOI: 10.1038/s41467-019-08940-5)
#Timepoint: D25
#Data: GSE125540
pscscore_1 = read.csv("psc-atac/pscscore_1.txt", sep = "\t")
pscscore_2 = read.csv("psc-atac/pscscore_2.txt", sep = "\t")
pscscore_3 = read.csv("psc-atac/pscscore_3.txt", sep = "\t")
pscscore_4 = read.csv("psc-atac/pscscore_4.txt", sep = "\t")
pscscore_5 = read.csv("psc-atac/pscscore_5.txt", sep = "\t")
pscscore_6 = read.csv("psc-atac/pscscore_6.txt", sep = "\t")
pscscore_7 = read.csv("psc-atac/pscscore_7.txt", sep = "\t")
pscscore_8 = read.csv("psc-atac/pscscore_8.txt", sep = "\t")
pscscore_9 = read.csv("psc-atac/pscscore_9.txt", sep = "\t")
pscscore_10 = read.csv("psc-atac/pscscore_10.txt", sep = "\t")

#Extract only the promoter-TSS peaks
c15_1_tss = c15_1[startsWith(as.character(c15_1$Annotation), "promoter-TSS"), ]
c15_2_tss = c15_2[startsWith(as.character(c15_2$Annotation), "promoter-TSS"), ]
c20_1_tss = c20_1[startsWith(as.character(c20_1$Annotation), "promoter-TSS"), ]
c20_2_tss = c20_2[startsWith(as.character(c20_2$Annotation), "promoter-TSS"), ]
h1_1_tss = h1_1[startsWith(as.character(h1_1$Annotation), "promoter-TSS"), ]
h1_2_tss = h1_2[startsWith(as.character(h1_2$Annotation), "promoter-TSS"), ]
h9_1_tss = h9_1[startsWith(as.character(h9_1$Annotation), "promoter-TSS"), ]
h9_2_tss = h9_2[startsWith(as.character(h9_2$Annotation), "promoter-TSS"), ]
rues2_1_tss = rues2_1[startsWith(as.character(rues2_1$Annotation), "promoter-TSS"), ]
rues2_2_tss = rues2_2[startsWith(as.character(rues2_2$Annotation), "promoter-TSS"), ]
pscscore_1_tss = pscscore_1[startsWith(as.character(pscscore_1$Annotation), "promoter-TSS"), ]
pscscore_2_tss = pscscore_2[startsWith(as.character(pscscore_2$Annotation), "promoter-TSS"), ]
pscscore_3_tss = pscscore_3[startsWith(as.character(pscscore_3$Annotation), "promoter-TSS"), ]
pscscore_4_tss = pscscore_4[startsWith(as.character(pscscore_4$Annotation), "promoter-TSS"), ]
pscscore_5_tss = pscscore_5[startsWith(as.character(pscscore_5$Annotation), "promoter-TSS"), ]
pscscore_6_tss = pscscore_6[startsWith(as.character(pscscore_6$Annotation), "promoter-TSS"), ]
pscscore_7_tss = pscscore_7[startsWith(as.character(pscscore_7$Annotation), "promoter-TSS"), ]
pscscore_8_tss = pscscore_8[startsWith(as.character(pscscore_8$Annotation), "promoter-TSS"), ]
pscscore_9_tss = pscscore_9[startsWith(as.character(pscscore_9$Annotation), "promoter-TSS"), ]
pscscore_10_tss = pscscore_10[startsWith(as.character(pscscore_10$Annotation), "promoter-TSS"), ]

#For each sample, compute the percentage of genes in the correctly regulated list and dysregulated list that have promoter-TSS peaks
liu_atac = as.data.frame(matrix(nrow = 8, ncol = 2))
rownames(liu_atac) = c("c15_1", "c15_2", "c20_1", "c20_2", "h1_1", "h1_2", "h9_1", "h9_2")
colnames(liu_atac) = c("Correct", "Dysregulated")
liu_atac$Correct = c(table(toupper(human_good_clean) %in% c15_1_tss$Gene.Name)[2]/length(human_good_clean), table(toupper(human_good_clean) %in% c15_2_tss$Gene.Name)[2]/length(human_good_clean), table(toupper(human_good_clean) %in% c20_1_tss$Gene.Name)[2]/length(human_good_clean), table(toupper(human_good_clean) %in% c20_2_tss$Gene.Name)[2]/length(human_good_clean), table(toupper(human_good_clean) %in% h1_1_tss$Gene.Name)[2]/length(human_good_clean), table(toupper(human_good_clean) %in% h1_2_tss$Gene.Name)[2]/length(human_good_clean), table(toupper(human_good_clean) %in% h9_1_tss$Gene.Name)[2]/length(human_good_clean), table(toupper(human_good_clean) %in% h9_2_tss$Gene.Name)[2]/length(human_good_clean))
liu_atac$Dysregulated= c(table(toupper(dysreg) %in% c15_1_tss$Gene.Name)[2]/length(dysreg), table(toupper(dysreg) %in% c15_2_tss$Gene.Name)[2]/length(dysreg), table(toupper(dysreg) %in% c20_1_tss$Gene.Name)[2]/length(dysreg), table(toupper(dysreg) %in% c20_2_tss$Gene.Name)[2]/length(dysreg), table(toupper(dysreg) %in% h1_1_tss$Gene.Name)[2]/length(dysreg), table(toupper(dysreg) %in% h1_2_tss$Gene.Name)[2]/length(dysreg), table(toupper(dysreg) %in% h9_1_tss$Gene.Name)[2]/length(dysreg), table(toupper(dysreg) %in% h9_2_tss$Gene.Name)[2]/length(dysreg))

bertero_atac = as.data.frame(matrix(nrow = 2, ncol = 2))
rownames(bertero_atac) = c("rues2_1", "rues2_2")
colnames(bertero_atac) = c("Correct", "Dysregulated")
bertero_atac$Correct = c(table(toupper(human_good_clean) %in% rues2_1_tss$Gene.Name)[2]/length(human_good_clean), table(toupper(human_good_clean) %in% rues2_2_tss$Gene.Name)[2]/length(human_good_clean))
bertero_atac$Dysregulated = c(table(toupper(dysreg) %in% rues2_1_tss$Gene.Name)[2]/length(dysreg), table(toupper(dysreg) %in% rues2_2_tss$Gene.Name)[2]/length(dysreg))

greenwald_atac = as.data.frame(matrix(nrow = 10, ncol = 2))
rownames(greenwald_atac) = c("pscscore_1", "pscscore_2", "pscscore_3", "pscscore_4", "pscscore_5", "pscscore_6", "pscscore_7", "pscscore_8", "pscscore_9", "pscscore_10")
colnames(greenwald_atac) = c("Correct", "Dysregulated")
greenwald_atac$Correct = c(table(toupper(human_good_clean) %in% pscscore_1_tss$Gene.Name)[2]/length(human_good_clean), table(toupper(human_good_clean) %in% pscscore_2_tss$Gene.Name)[2]/length(human_good_clean), table(toupper(human_good_clean) %in% pscscore_3_tss$Gene.Name)[2]/length(human_good_clean), table(toupper(human_good_clean) %in% pscscore_4_tss$Gene.Name)[2]/length(human_good_clean), table(toupper(human_good_clean) %in% pscscore_5_tss$Gene.Name)[2]/length(human_good_clean), table(toupper(human_good_clean) %in% pscscore_6_tss$Gene.Name)[2]/length(human_good_clean), table(toupper(human_good_clean) %in% pscscore_7_tss$Gene.Name)[2]/length(human_good_clean), table(toupper(human_good_clean) %in% pscscore_8_tss$Gene.Name)[2]/length(human_good_clean), table(toupper(human_good_clean) %in% pscscore_9_tss$Gene.Name)[2]/length(human_good_clean), table(toupper(human_good_clean) %in% pscscore_10_tss$Gene.Name)[2]/length(human_good_clean))
greenwald_atac$Dysregulated = c(table(toupper(dysreg) %in% pscscore_1_tss$Gene.Name)[2]/length(dysreg), table(toupper(dysreg) %in% pscscore_2_tss$Gene.Name)[2]/length(dysreg), table(toupper(dysreg) %in% pscscore_3_tss$Gene.Name)[2]/length(dysreg), table(toupper(dysreg) %in% pscscore_4_tss$Gene.Name)[2]/length(dysreg), table(toupper(dysreg) %in% pscscore_5_tss$Gene.Name)[2]/length(dysreg), table(toupper(dysreg) %in% pscscore_6_tss$Gene.Name)[2]/length(dysreg), table(toupper(dysreg) %in% pscscore_7_tss$Gene.Name)[2]/length(dysreg), table(toupper(dysreg) %in% pscscore_8_tss$Gene.Name)[2]/length(dysreg), table(toupper(dysreg) %in% pscscore_9_tss$Gene.Name)[2]/length(dysreg), table(toupper(dysreg) %in% pscscore_10_tss$Gene.Name)[2]/length(dysreg))

#Clean up and combine the studies
liu_atac$study = "Liu (D30)"
liu_atac$sample = rownames(liu_atac)
bertero_atac$study = "Bertero and Fields (D14)"
bertero_atac$sample = rownames(bertero_atac)
greenwald_atac$study = "Greenwald and Li (D25)"
greenwald_atac$sample = rownames(greenwald_atac)
combined_atac = melt(rbind(liu_atac, bertero_atac, greenwald_atac))

#####Cell and Tissue Engineering Dataset Analysis
#Our next goal was to look at how effectively cellular and tissue engineering methods could improve PSC-CM maturation in a biomimetic manner. We identified datasets that compared perturbed PSC-CMs to controls, with the assumption that any gene changes should occur in the same direction as in vivo CM maturation.

###LOAD IN DATASETS HERE
#Please note that this section is to show users how I processed input data, in case the users want to repeat the steps from the raw data sources. However, I'm not reproviding the files to do this (it would be a bit redundant). However, I do provide a workspace with all of the datasets that can be loaded in below.

temp = list.files(path = "~/Downloads/GSE116574_RAW/", pattern="*.gz", full.names = TRUE)
myfiles = lapply(temp, read.delim)
branco_data = data.frame(row.names = myfiles[[1]][, 1], agg_1 = myfiles[[1]][, 2], agg_2 = myfiles[[2]][, 2], agg_3 = myfiles[[3]][, 2], ctrl_1 = myfiles[[4]][, 2], ctrl_2 = myfiles[[5]][, 2], ctrl_3 = myfiles[[6]][, 2])
#Zhao et al.
zhao_data = read.table("~/Downloads/GSE114976_count_data.csv.gz", as.is = TRUE, sep = ",", row.names = 1, header = TRUE)
zhao_data = zhao_data[!duplicated(zhao_data$SYMBOL), ]
zhao_data = zhao_data[!is.na(zhao_data$SYMBOL), ]
rownames(zhao_data) = zhao_data$SYMBOL
zhao_data = zhao_data[, 4:10]
zhao_data = round(zhao_data)
#Kuppusamy et al.
kuppusamy_data = read.table("~/Downloads/let7.tsv.gz", as.is= TRUE, sep = "\t", header = TRUE, row.names = 1)
kuppusamy_data = kuppusamy_data[, 1:11]
kuppusamy_pheno = data.frame(row.names = colnames(kuppusamy_data), condition = paste(sapply(strsplit(colnames(kuppusamy_data), "[.]"), "[[", 2), "_", sapply(strsplit(colnames(kuppusamy_data), "[.]"), "[[", 3), sep = ""))
#Feyen et al.
feyen_data = read.table("~/Downloads/GSE151279_counts.csv.gz", as.is = TRUE, sep = ",", row.names = 1, header = TRUE)
#Lam et al. (please beware as this one may be big)
lam_ctrl1 = readMM("~/Downloads/GSE157157_RAW/GSM4756813_Control1hCM_matrix.mtx.gz")
colnames(lam_ctrl1) = read.table("~/Downloads/GSE157157_RAW/GSM4756813_Control1hCM_barcodes.tsv.gz", as.is = TRUE)$V1
rownames(lam_ctrl1) = read.table("~/Downloads/GSE157157_RAW/GSM4756813_Control1hCM_features.tsv.gz", as.is = TRUE)$V1
lam_hcas1 = readMM("~/Downloads/GSE157157_RAW/GSM4756814_Control1hCAS_matrix.mtx.gz")
colnames(lam_hcas1) = read.table("~/Downloads/GSE157157_RAW/GSM4756814_Control1hCAS_barcodes.tsv.gz", as.is = TRUE)$V1
rownames(lam_hcas1) = read.table("~/Downloads/GSE157157_RAW/GSM4756814_Control1hCAS_features.tsv.gz", as.is = TRUE)$V1
lam_hcts1 = readMM("~/Downloads/GSE157157_RAW/GSM4756815_Control1hCTS_matrix.mtx.gz")
colnames(lam_hcts1) = read.table("~/Downloads/GSE157157_RAW/GSM4756815_Control1hCTS_barcodes.tsv.gz", as.is = TRUE)$V1
rownames(lam_hcts1) = read.table("~/Downloads/GSE157157_RAW/GSM4756815_Control1hCTS_features.tsv.gz", as.is = TRUE)$V1
lam_ctrl2 = readMM("~/Downloads/GSE157157_RAW/GSM4756816_Control2hCM_matrix.mtx.gz")
colnames(lam_ctrl2) = read.table("~/Downloads/GSE157157_RAW/GSM4756816_Control2hCM_barcodes.tsv.gz", as.is = TRUE)$V1
rownames(lam_ctrl2) = read.table("~/Downloads/GSE157157_RAW/GSM4756816_Control2hCM_features.tsv.gz", as.is = TRUE)$V1
lam_hcas2 = readMM("~/Downloads/GSE157157_RAW/GSM4756817_Control2hCAS_matrix.mtx.gz")
colnames(lam_hcas2) = read.table("~/Downloads/GSE157157_RAW/GSM4756817_Control2hCAS_barcodes.tsv.gz", as.is = TRUE)$V1
rownames(lam_hcas2) = read.table("~/Downloads/GSE157157_RAW/GSM4756817_Control2hCAS_features.tsv.gz", as.is = TRUE)$V1
lam_hcts2 = readMM("~/Downloads/GSE157157_RAW/GSM4756818_Control2hCTS_matrix.mtx.gz")
colnames(lam_hcts2) = read.table("~/Downloads/GSE157157_RAW/GSM4756818_Control2hCTS_barcodes.tsv.gz", as.is = TRUE)$V1
rownames(lam_hcts2) = read.table("~/Downloads/GSE157157_RAW/GSM4756818_Control2hCTS_features.tsv.gz", as.is = TRUE)$V1
lam_names = Reduce(intersect, list(rownames(lam_ctrl1), rownames(lam_ctrl2), rownames(lam_hcas1), rownames(lam_hcas2), rownames(lam_hcts1), rownames(lam_hcts2)))
colnames(lam_ctrl1) = paste(colnames(lam_ctrl1), "ctrl.1", sep = "_")
colnames(lam_ctrl2) = paste(colnames(lam_ctrl2), "ctrl.2", sep = "_")
colnames(lam_hcas1) = paste(colnames(lam_hcas1), "hcas.1", sep = "_")
colnames(lam_hcas2) = paste(colnames(lam_hcas2), "hcas.2", sep = "_")
colnames(lam_hcts1) = paste(colnames(lam_hcts1), "hcts.1", sep = "_")
colnames(lam_hcts2) = paste(colnames(lam_hcts2), "hcts.2", sep = "_")
lam_data = cbind(lam_ctrl1[lam_names, ], lam_ctrl2[lam_names, ], lam_hcas1[lam_names, ], lam_hcas2[lam_names, ], lam_hcts1[lam_names, ], lam_hcts2[lam_names, ])
lam_pheno = data.frame(row.names = colnames(lam_data), timepoint = "D30", full_condition = sapply(strsplit(colnames(lam_data), "_"), "[[", 2))
lam_pheno$condition = sapply(strsplit(as.character(lam_pheno$full_condition), "[.]"), "[[", 1)
lam_pheno$sample = sapply(strsplit(as.character(lam_pheno$full_condition), "[.]"), "[[", 2)
lam_data = rename_genes(lam_data, species = "human")
lam_temp = data_qc("lam_data", "Lam", timepoint_list = lam_pheno$timepoint, sample_type = "directed differentiation", species = "human", isolation = "10x Chromium V2", sequencing = "10x Chromium", mapping = "CellRanger", datatype = "UMIs", doi = "doi:10.1161/JAHA.120.016528", other_meta = lam_pheno$condition)

#Giacomelli et al.
#For this dataset, I drew the data from our previous database (see the entropy paper).

#Or, just load in the data from this workspace.
load(te_data.RData)

#For all the bulk datasets, for consistency, I set a threshold of q-val < 0.05 and abs(log2FC) > 0.5. For the single cell datasets I did not use a fold change threshold.
threshold = 0.5

###3d aggregate Culture
#Study: Branco et al. (DOI:10.1038/s41598-019-45047-9)
#Design: Microwell 3D aggregation vs. monolayer culture, compared at D20
#Data: GSE116574

branco_dds = DESeqDataSetFromMatrix(countData = branco_data, colData = data.frame(row.names = colnames(branco_data), condition = sapply(strsplit(colnames(branco_data), "_"), "[[", 1)), design = ~condition)
branco_dds$condition = factor(branco_dds$condition, levels = c("ctrl", "agg"))
branco_dds = DESeq(branco_dds)
branco_res = results(branco_dds)
branco_res = branco_res[order(branco_res$padj), ]
branco_res = branco_res[!is.na(branco_res$padj), ]

branco_down = rownames(branco_res)[branco_res$padj < 0.05 & branco_res$log2FoldChange <= (-1 * threshold)]
branco_down_mouse = convertGenes(branco_down)
branco_up = rownames(branco_res)[branco_res$padj < 0.05 & branco_res$log2FoldChange > (threshold)]
branco_up_mouse = convertGenes(branco_up)

table(c(branco_down_mouse, branco_up_mouse) %in% pseudo_ids_invivo_perinatal)
#455
table(branco_down_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "down"])
#293
table(branco_up_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "up"])
#16
branco_good = c(branco_down_mouse[branco_down_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "down"]], branco_up_mouse[branco_up_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "up"]])
#So 309/455 genes are expressed in correct direction
table(branco_good %in% dysreg) #32
table(branco_good %in% c(pseudo_ids_invitro_good, human_good)) #151
table(branco_good %in% c(human_good_new)) #191
plot(unlist(lapply(subset_diff, function(x) {sum(branco_good %in% x)/length(x)})))
table(invivo_genes[branco_good, ]$cluster)/table(invivo_genes$cluster)

###Biowire (suspended microtissues with electrical stimulation)
#Study: Zhao et al. (DOI:10.1016/j.cell.2018.11.042)
#Design: Stimulated ventricular CMs on biowire vs unstimulated control; not exactly sure when compared, I think D20
#Data: GSE114976

zhao_dds = DESeqDataSetFromMatrix(countData = zhao_data, colData = data.frame(row.names = colnames(zhao_data), group = c(rep("Stim", 3), rep("Ctrl", 4))), design = ~group)
zhao_dds = DESeq(zhao_dds)
zhao_res = results(zhao_dds)
zhao_res = zhao_res[order(zhao_res$padj), ]
zhao_res = zhao_res[!is.na(zhao_res$padj), ]

zhao_down = rownames(zhao_res)[zhao_res$padj < 0.05 & zhao_res$log2FoldChange <= (-1 * threshold)]
zhao_down_mouse = convertGenes(zhao_down)
zhao_up = rownames(zhao_res)[zhao_res$padj < 0.05 & zhao_res$log2FoldChange > (threshold)]
zhao_up_mouse = convertGenes(zhao_up)

table(c(zhao_down_mouse, zhao_up_mouse) %in% pseudo_ids_invivo_perinatal)
#762
table(zhao_down_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "down"])
#281
table(zhao_up_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "up"])
#175
zhao_good = c(zhao_down_mouse[zhao_down_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "down"]], zhao_up_mouse[zhao_up_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "up"]])
#So 456/762 genes are expressed in correct direction
table(zhao_good %in% dysreg) #107
table(zhao_good %in% c(pseudo_ids_invitro_good, human_good)) #176
table(zhao_good %in% c(human_good_new)) #1915
plot(unlist(lapply(subset_diff, function(x) {sum(zhao_good %in% x)/length(x)})))
table(invivo_genes[zhao_good, ]$cluster)/table(invivo_genes$cluster)

###Let7 Overexpression
#Study: Kuppusamy et al. (DOI:10.1073/pnas.1424042112)
#Design: Let7 OE PSC-CMs vs Ctrl CMs with EV, compared at D30
#Data: GSE62913

kupp_let_dds = DESeqDataSetFromMatrix(countData = kuppusamy_data[, 7:11], colData = kuppusamy_pheno[7:11, , drop = FALSE], design = ~condition)
kupp_let_dds = DESeq(kupp_let_dds)
kupp_let_res = results(kupp_let_dds)
kupp_let_res = kupp_let_res[order(kupp_let_res$padj), ]
kupp_let_res = kupp_let_res[!is.na(kupp_let_res$padj), ]

kupp_let_down = rownames(kupp_let_res)[kupp_let_res$padj < 0.05 & kupp_let_res$log2FoldChange <= (-1*threshold)]
kupp_let_down_mouse = convertGenes(kupp_let_down)
kupp_let_up = rownames(kupp_let_res)[kupp_let_res$padj < 0.05 & kupp_let_res$log2FoldChange >= threshold]
kupp_let_up_mouse = convertGenes(kupp_let_up)

table(c(kupp_let_down_mouse, kupp_let_up_mouse) %in% pseudo_ids_invivo_perinatal)
#554
table(kupp_let_down_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "down"])
#236
table(kupp_let_up_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "up"])
#147
kupp_let_good = c(kupp_let_down_mouse[kupp_let_down_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "down"]], kupp_let_up_mouse[kupp_let_up_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "up"]])
#So 383/554 genes are expressed in correct direction
table(kupp_let_good %in% dysreg) #60
table(kupp_let_good %in% c(pseudo_ids_invitro_good, human_good)) #209
table(kupp_let_good %in% c(human_good_new)) #227
plot(unlist(lapply(subset_diff, function(x) {sum(kupp_let_good %in% x)/length(x)})))
table(invivo_genes[kupp_let_good, ]$cluster)/table(invivo_genes$cluster)

###Long-term Culture
#Study: Kuppusamy et al. (DOI:10.1073/pnas.1424042112)
#Design: 1 year vs D20 PSC-CMs
#Data: GSE62913

kupp_age_dds = DESeqDataSetFromMatrix(countData = kuppusamy_data[, 1:6], colData = kuppusamy_pheno[1:6, , drop = FALSE], design = ~condition)
kupp_age_dds$condition = factor(kupp_age_dds$condition, levels = c("h7_day20", "h7_1yr"))
kupp_age_dds = DESeq(kupp_age_dds)
kupp_age_res = results(kupp_age_dds)
kupp_age_res = kupp_age_res[order(kupp_age_res$padj), ]
kupp_age_res = kupp_age_res[!is.na(kupp_age_res$padj), ]

#The default test has d20 on top
kupp_age_down = rownames(kupp_age_res)[kupp_age_res$padj < 0.05 & kupp_age_res$log2FoldChange <= (-1*threshold)]
kupp_age_down_mouse = convertGenes(kupp_age_down)
kupp_age_up = rownames(kupp_age_res)[kupp_age_res$padj < 0.05 & kupp_age_res$log2FoldChange >= threshold]
kupp_age_up_mouse = convertGenes(kupp_age_up)

table(c(kupp_age_down_mouse, kupp_age_up_mouse) %in% pseudo_ids_invivo_perinatal)
#754
table(kupp_age_down_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "down"])
#243
table(kupp_age_up_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "up"])
#102
kupp_age_good = c(kupp_age_down_mouse[kupp_age_down_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "down"]], kupp_age_up_mouse[kupp_age_up_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "up"]])
#So 345/754 genes are expressed in correct direction
table(kupp_age_good %in% dysreg) #70
table(kupp_age_good %in% c(pseudo_ids_invitro_good, human_good)) #140
table(kupp_age_good %in% c(human_good_new)) #168
plot(unlist(lapply(subset_diff, function(x) {sum(kupp_age_good %in% x)/length(x)})))
table(invivo_genes[kupp_age_good, ]$cluster)/table(invivo_genes$cluster)

###Metabolic Media
#Study: Feyen, McKeithan, and Bruynell et al. (DOI:10.1016/j.celrep.2020.107925)
#Design: Metabolic Media vs Control (RPMI/B27 Media) at D40
#Data: GSE151279

feyen_dds = DESeqDataSetFromMatrix(countData = feyen_data[, 4:10], colData = data.frame(row.names = colnames(feyen_data)[4:10], condition = sapply(strsplit(colnames(feyen_data)[4:10], "_"), "[[", 1)), design = ~condition)
feyen_dds$condition = factor(feyen_dds$condition, levels = c("W3RC", "W3MC"))
feyen_dds = DESeq(feyen_dds)
feyen_res = results(feyen_dds)
feyen_res = feyen_res[order(feyen_res$padj), ]
feyen_res = feyen_res[!is.na(feyen_res$padj), ]

feyen_down = rownames(feyen_res)[feyen_res$padj < 0.05 & feyen_res$log2FoldChange <= (-1 * threshold)]
feyen_down_mouse = convertGenes(feyen_down)
feyen_up = rownames(feyen_res)[feyen_res$padj < 0.05 & feyen_res$log2FoldChange > (threshold)]
feyen_up_mouse = convertGenes(feyen_up)

table(c(feyen_down_mouse, feyen_up_mouse) %in% pseudo_ids_invivo_perinatal)
#264
table(feyen_down_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "down"])
#76
table(feyen_up_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "up"])
#127
feyen_good = c(feyen_down_mouse[feyen_down_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "down"]], feyen_up_mouse[feyen_up_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "up"]])
#So 203/264 genes are expressed in correct direction
table(feyen_good %in% dysreg) #53
table(feyen_good %in% c(pseudo_ids_invitro_good, human_good)) #96
table(feyen_good %in% c(human_good_new)) #100
plot(unlist(lapply(subset_diff, function(x) {sum(feyen_good %in% x)/length(x)})))
table(invivo_genes[feyen_good, ]$cluster)/table(invivo_genes$cluster)

###Engineered Microtissue
#Study: Lam et al. (DOI:10.1161/JAHA.120.016528)
#Design: PSC-CMs in ECTs vs PSC-CMs in monolayer, compared at D30
#Data: GSE157157

gene_fits <- fit_models(cds_lam[fData(cds_lam)$num_cells_expressed > 400, cds_lam@colData$other_meta %in% c("ctrl", "hcts")], model_formula_str = "~other_meta")
fit_coefs <- coefficient_table(gene_fits)
hcts_terms <- fit_coefs %>% filter(term == "other_metahcts")
hcts_terms = hcts_terms[order(hcts_terms$q_value), ]
hcts_ids= hcts_terms[hcts_terms$status != "FAIL" & hcts_terms$q_value < 0.05, ]$symbol
hcts_ids_mouse = convertGenes(hcts_ids)

hcts_ids_up= hcts_terms[hcts_terms$status != "FAIL" & hcts_terms$q_value < 0.05 & hcts_terms$normalized_effect >= 0, ]$symbol
hcts_ids_up_mouse = convertGenes(hcts_ids_up)
hcts_ids_down = hcts_terms[hcts_terms$status != "FAIL" & hcts_terms$q_value < 0.05 & hcts_terms$normalized_effect <= 0, ]$symbol
hcts_ids_down_mouse = convertGenes(hcts_ids_down)

table(c(hcts_ids_mouse) %in% pseudo_ids_invivo_perinatal)
#1012
table(hcts_ids_down_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "down"])
#536
table(hcts_ids_up_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "up"])
#142
hcts_good = c(hcts_ids_down_mouse[hcts_ids_down_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "down"]], hcts_ids_up_mouse[hcts_ids_up_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "up"]])
#So 678/1012 genes are expressed in correct direction
table(hcts_good %in% dysreg) #147
table(hcts_good %in% c(pseudo_ids_invitro_good, human_good)) #385
table(hcts_good %in% c(human_good_new)) #448
plot(unlist(lapply(subset_diff, function(x) {sum(hcts_good %in% x)/length(x)})))
table(invivo_genes[hcts_good, ]$cluster)/table(invivo_genes$cluster)

###Anistropic Sheets
#Study: Lam et al. (DOI:10.1161/JAHA.120.016528)
#Design: PSC-CMs in patterned (anisotropic) sheets vs PSC-CMs in monolayer, compared at D30
#Data: GSE157157

cds_lam <- new_cell_data_set(lam_data[ , lam_temp$good_cell == TRUE], cell_metadata = lam_temp[lam_temp$good_cell == TRUE, ], gene_metadata = G_list_human[rownames(lam_data), ])
fData(cds_lam)$num_cells_expressed = rowSums(exprs(cds_lam) >= 1)
gene_fits <- fit_models(cds_lam[fData(cds_lam)$num_cells_expressed > 400, cds_lam@colData$other_meta %in% c("ctrl", "hcas")], model_formula_str = "~other_meta")
fit_coefs <- coefficient_table(gene_fits)
hcas_terms <- fit_coefs %>% filter(term == "other_metahcas")
hcas_terms = hcas_terms[order(hcas_terms$q_value), ]
hcas_ids= hcas_terms[hcas_terms$status != "FAIL" & hcas_terms$q_value < 0.05, ]$symbol
hcas_ids_mouse = convertGenes(hcas_ids)

hcas_ids_up= hcas_terms[hcas_terms$status != "FAIL" & hcas_terms$q_value < 0.05 & hcas_terms$normalized_effect >= 0, ]$symbol
hcas_ids_up_mouse = convertGenes(hcas_ids_up)
hcas_ids_down = hcas_terms[hcas_terms$status != "FAIL" & hcas_terms$q_value < 0.05 & hcas_terms$normalized_effect <= 0, ]$symbol
hcas_ids_down_mouse = convertGenes(hcas_ids_down)

table(c(hcas_ids_mouse) %in% pseudo_ids_invivo_perinatal)
#990
table(hcas_ids_down_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "down"])
#507
table(hcas_ids_up_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "up"])
#83
hcas_good = c(hcas_ids_down_mouse[hcas_ids_down_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "down"]], hcas_ids_up_mouse[hcas_ids_up_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "up"]])
#So 590/990 genes are expressed in correct direction
table(hcas_good %in% dysreg) #86
table(hcas_good %in% c(pseudo_ids_invitro_good, human_good)) #345
table(hcas_good %in% c(human_good_new)) #425
plot(unlist(lapply(subset_diff, function(x) {sum(hcas_good %in% x)/length(x)})))
table(invivo_genes[hcas_good, ]$cluster)/table(invivo_genes$cluster)

###Co-culture
#Study: Giacomelli, Meraviglia, and Campostini et al. (DOI:10.1016/j.stem.2020.05.004)
#Design: PSC-CMs co-cultued with PSC-FBs and PSC-ECs (possibly with some type of 3d architecture?) vs PSM-CMs in monolayer, compared at D21
#Data: GSE147694

giacomelli_pheno = combined_datasets[combined_datasets$data == "giacomelli_cm_data", ]
rownames(giacomelli_pheno) = giacomelli_pheno$cellname
cds_giacomelli <- new_cell_data_set(giacomelli_cm_data[ , giacomelli_pheno$good_cell == TRUE], cell_metadata = giacomelli_pheno[giacomelli_pheno$good_cell == TRUE, ], gene_metadata = G_list_human[rownames(giacomelli_cm_data), ])
fData(cds_giacomelli)$num_cells_expressed = rowSums(exprs(cds_giacomelli) >= 1)
gene_fits <- fit_models(cds_giacomelli[fData(cds_giacomelli)$num_cells_expressed > 400, ], model_formula_str = "~other_meta")
fit_coefs <- coefficient_table(gene_fits)
terms <- fit_coefs %>% filter(term == "other_metaCMEC")
terms = terms[order(terms$q_value), ]
giacomelli_ids= terms[terms$status != "FAIL" & terms$q_value < 0.05, ]$symbol
giacomelli_ids_mouse = convertGenes(giacomelli_ids)

giacomelli_ids_up= terms[terms$status != "FAIL" & terms$q_value < 0.05 & terms$normalized_effect >= 0, ]$symbol
giacomelli_ids_up_mouse = convertGenes(giacomelli_ids_up)
giacomelli_ids_down = terms[terms$status != "FAIL" & terms$q_value < 0.05 & terms$normalized_effect <= 0, ]$symbol
giacomelli_ids_down_mouse = convertGenes(giacomelli_ids_down)

table(c(giacomelli_ids_mouse) %in% pseudo_ids_invivo_perinatal)
#701
table(giacomelli_ids_down_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "down"])
#262
table(giacomelli_ids_up_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "up"])
#100
giacomelli_good = c(giacomelli_ids_down_mouse[giacomelli_ids_down_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "down"]], giacomelli_ids_up_mouse[giacomelli_ids_up_mouse %in% rownames(invivo_genes)[invivo_genes$direction == "up"]])
#So 362/701 genes are expressed in correct direction
table(giacomelli_good %in% dysreg) #99
table(giacomelli_good %in% c(pseudo_ids_invitro_good, human_good)) #187
table(giacomelli_good %in% c(human_good_new)) #199
plot(unlist(lapply(subset_diff, function(x) {sum(giacomelli_good %in% x)/length(x)})))
table(invivo_genes[giacomelli_good, ]$cluster)/table(invivo_genes$cluster)

###TF Analysis
#For TF Analysis, I used WebGestalt for enrichment of downstream targets, as in previous sections (using the same settings). Because our primary interest is in the dysregulated genes, I focused only on the enrichment of those TFs, using the max enrichment ratio for each TF. How did I go about getting those max enrichments? Manually, because to be honest, at this exact moment, I'm too lazy to come up with something smarter.

te_tf = matrix(nrow = 9, ncol = 8)
rownames(te_tf) = c("Srf", "Yy1", "Jun", "Nfe2l2", "Sox9", "Ppara", "Esrra", "Nrf1", "Mef2a")
colnames(te_tf) = c("3D Agg.", "ES", "Let7 OE", "LTC", "MM", "EMT", "AIS", "CC")
te_agg = read.table("TF_TE_new/3dagg/results.txt", as.is = TRUE, row.names = 1, header = TRUE)
te_agg = te_agg[order(te_agg$enrichmentRatio, decreasing = TRUE), ]
te_es = read.table("TF_TE_new/es/results.txt", as.is = TRUE, row.names = 1, header = TRUE)
te_es = te_es[order(te_es$enrichmentRatio, decreasing = TRUE), ]
te_let7oe = read.table("TF_TE_new/let7oe/results.txt", as.is = TRUE, row.names = 1, header = TRUE)
te_let7oe = te_let7oe[order(te_let7oe$enrichmentRatio, decreasing = TRUE), ]
te_ltc = read.table("TF_TE_new/ltc/results.txt", as.is = TRUE, row.names = 1, header = TRUE)
te_ltc = te_ltc[order(te_ltc$enrichmentRatio, decreasing = TRUE), ]
te_mm = read.table("TF_TE_new/mm/results.txt", as.is = TRUE, row.names = 1, header = TRUE)
te_mm = te_mm[order(te_mm$enrichmentRatio, decreasing = TRUE), ]
te_emt = read.table("TF_TE_new/emt/results.txt", as.is = TRUE, row.names = 1, header = TRUE)
te_emt = te_emt[order(te_emt$enrichmentRatio, decreasing = TRUE), ]
te_ais = read.table("TF_TE_new/ais/results.txt", as.is = TRUE, row.names = 1, header = TRUE)
te_ais = te_ais[order(te_ais$enrichmentRatio, decreasing = TRUE), ]
te_cc = read.table("TF_TE_new/cc/results.txt", as.is = TRUE, row.names = 1, header = TRUE)
te_cc = te_cc[order(te_cc$enrichmentRatio, decreasing = TRUE), ]

#Here is a quick function to jut get all of the results for each TF quickly.
dysreg_ident = function(enrichment_results){
  srf = min(which(grepl("SRF", rownames(enrichment_results)) & !grepl("RSRF", rownames(enrichment_results))))
  yy1 = min(which(grepl("YY1", rownames(enrichment_results))))
  jun = min(which(grepl("AP1FJ", rownames(enrichment_results))))
  nrf2 = min(which(grepl("NRF2", rownames(enrichment_results))))
  sox9 = min(which(grepl("SOX9", rownames(enrichment_results))))
  ppara = min(which(grepl("PPAR", rownames(enrichment_results))))
  esrra = min(which(grepl("ERR1", rownames(enrichment_results))))
  nrf1 = min(which(grepl("NRF1", rownames(enrichment_results))))
  mef2a = min(which(grepl("MEF2", rownames(enrichment_results)) | grepl("RSRF", rownames(enrichment_results))))
  return(c(enrichment_results[srf, ]$enrichmentRatio, enrichment_results[yy1, ]$enrichmentRatio, enrichment_results[jun, ]$enrichmentRatio, enrichment_results[nrf2, ]$enrichmentRatio, enrichment_results[sox9, ]$enrichmentRatio, enrichment_results[ppara, ]$enrichmentRatio, enrichment_results[esrra, ]$enrichmentRatio, enrichment_results[nrf1, ]$enrichmentRatio, enrichment_results[mef2a, ]$enrichmentRatio))
}

te_tf[, 1] = dysreg_ident(te_agg)
te_tf[, 2] = dysreg_ident(te_es)
te_tf[, 3] = dysreg_ident(te_let7oe)
te_tf[, 4] = dysreg_ident(te_ltc)
te_tf[, 5] = dysreg_ident(te_mm)
te_tf[, 6] = dysreg_ident(te_emt)
te_tf[, 7] = dysreg_ident(te_ais)
te_tf[, 8] = dysreg_ident(te_cc)
te_tf[is.na(te_tf)] = 0

#####RNA Velocity Analysis
#To do the RNA Velocity analysis, we use the scvelo package (Bergen 2019). This package is absolutely great, but also only in Python, so we need to export a couple of files from here to Python. Please check the separate Python code for how the velocity analysis was done. Once done, we can re-import some of the data here for plotting (mostly for consistency - though as of now, some plots can only be made in Python).  
#Our approach to RNA Velocity analysis uses the UMAP trajectories generated above, as well as the differential gene expression lists identified above (rather than using scvelo to identify "high dispersion genes"). I think this approach makes more intuitive sense since we are interested in the dynamics of genes that are differentially regulated during this period.

###Generating spliced/unspliced counts tables for export
spliced = readMM("spliced.mtx")
unspliced = readMM("unspliced.mtx")
rownames(spliced) = read.table("spliced.barcodes.txt")$V1
rownames(unspliced) = read.table("unspliced.barcodes.txt")$V1
colnames(spliced) = read.table("spliced.genes.txt")$V1
colnames(unspliced) = read.table("unspliced.genes.txt")$V1
ttg = read.table("mouse_velocity_t2g.txt", as.is = TRUE)
gene_names = ttg[match(colnames(spliced), ttg$V2), ]$V3
colnames(spliced) = gene_names
colnames(unspliced) = gene_names
spliced = spliced[, !duplicated(colnames(spliced))]
unspliced = unspliced[, !duplicated(colnames(unspliced))]
barcodes = colnames(cds_combined)
spliced = spliced[barcodes, ]
unspliced = unspliced[barcodes, ]
spliced = t(mito_correct(t(spliced)))
unspliced = t(mito_correct(t(unspliced)))
spliced = spliced[, colnames(spliced) %in% G_list[G_list$gene_biotype %in% c("lincRNA", "antisense", "protein_coding"), ]$symbol]
unspliced = unspliced[, colnames(unspliced) %in% G_list[G_list$gene_biotype %in% c("lincRNA", "antisense", "protein_coding"), ]$symbol]
spliced = spliced[, !colnames(spliced) %in% c("Gm42418", "AY036118", "Gm26822", "Grip2", "Zrsr1")]
unspliced = unspliced[, !colnames(unspliced) %in% c("Gm42418", "AY036118", "Gm26822", "Grip2", "Zrsr1")]
write.table(as.matrix(t(spliced)), "python_spliced.csv", row.names = TRUE, col.names = FALSE, sep = ",")
write.table(as.matrix(t(unspliced)), "python_unspliced.csv", row.names = TRUE, col.names = FALSE, sep = ",")
rm(spliced, unspliced, gene_names, barcodes)

###Exporting files needed for the velocity analysis
#Differentially expressed genes
write(pseudo_ids_invivo, "pseudo_ids_invivo.txt", sep = ",")
write(pseudo_ids_invitro, "pseudo_ids_invitro.txt", sep = ",")

#mnnCorrected Aligned PCA plots
write.table(cds_combined@reducedDims$Aligned, "python_pca.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(cds_invivo@reducedDims$Aligned, "python_pca_invivo.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(cds_invitro@reducedDims$Aligned, "python_pca_invitro.csv", sep = ",", row.names = FALSE, col.names = FALSE)

#UMAP trajectory Plots
write.table(cds_combined@reducedDims$UMAP, "python_umap.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(cds_invivo@reducedDims$UMAP, "python_umap_invivo.csv", sep = ",", row.names = FALSE, col.names = FALSE)
write.table(cds_invitro@reducedDims$UMAP, "python_umap_invitro.csv", sep = ",", row.names = FALSE, col.names = FALSE)

#Phenotype table
write.csv(as.data.frame(cds_combined@colData), "~python_pheno.csv")

###Importing files needed for the velocity analysis in R
invivo_velocity = read.csv("velocity_results_invivo.csv", as.is = TRUE, row.names = 1)
invivo_velocity$timepoint = factor(invivo_velocity$timepoint, levels = timepoint_levels)
invitro_velocity = read.csv("velocity_results_invitro.csv", as.is = TRUE, row.names = 1)
invitro_velocity$timepoint = factor(invitro_velocity$timepoint, levels = timepoint_levels)

KwonLab378
tf_table = read.table("~/Documents/Research/BigRef/Mus_musculus_TF.txt", as.is = TRUE, header = TRUE, sep = "\t")
