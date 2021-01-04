#Trajectory Reconstruction of PSC-CM vs Endogenous CMs
#Chulan Kwon Laboratory
#Primary author: Suraj Kannan
#December 29, 2020
#Document: Figures 3.3

#Load in all necessary libraries
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
library(scales)
library(gridExtra)
library(VennDiagram)
library(lubridate)


#Please note - this file is meant to accompany the corresponding code file for this project; it assumes that all of the code there has been run before any of this code is run. This code is simply intended to produce manuscript-ready figures.

#So initially, the figures here were made with my computer screen in mind. As it turns out, with svg formatting, it is better to start small and then go bigger, rather than the other way around. Thus, the size of everything needed to be changed. The fig_factor below allows for this - if set to 1, the figures will be appropriate for a computer screen/ppt. If set to 2.5, they will be appropriate for the manuscript.
fig_factor = 1

#####Fig 01/Supp 01/Supp 02
#The following section outlines all of the panels associated with Figure 1, which details the in vivo trajectory reconstruction and the quantification of in vivo gene dynamics. Also includes Supplementary Figure 1 (additional details about in vivo trajectory) and Supplementary Figure 2 (GO terms for the various gene clusters).

###Fig 01a
#This is a workflow image and thus there is nothing to plot.

###Fig 01b
#Reconstructed trajectory of the in vivo CMs
plot_cells(cds_invivo, label_groups_by_cluster = FALSE, label_cell_groups = FALSE, color_cells_by = "timepoint", cell_size = 1/fig_factor, show_trajectory_graph = FALSE) + xlab("Component 1") + ylab("Component 2") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1/fig_factor, 'lines')) + guides(color = guide_legend(title = "Timepoint", override.aes = list(size = 4/fig_factor)))

###Fig 01c
#Boxplot of pseudotime scores assigned by Monocle to in vivo trajectory by timepoint
ggplot(as.data.frame(cds_invivo@colData), aes(x = timepoint, y = pseudotime, fill = timepoint)) + geom_boxplot(lwd = 1/fig_factor, outlier.size = 0) + geom_jitter(size = 1.5/fig_factor, alpha = 0.3) + coord_flip() + xlab("Timepoint") + ylab("Pseudotime") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.position = "none")

###Fig 01d
#Lineplot showing what percent of in vivo differentially expressed genes become differentially expressed by each pseudotime upper bound. Lets say diff_genes() outputs the number of genes differentially expressed in some pseudotime interval. Then here, we are plotting diff_genes([0,n])/diff_genes([0,60]), where n is the pseudotime upper bound plotted on the x-axis.
lineplot_temp = data.frame(pseudotime = seq(1, 60), kannan = unlist(lapply(subset_diff, function(x){length(x)/length(pseudo_ids_invivo)})))
ggplot(lineplot_temp, aes(x = as.numeric(pseudotime), y = kannan * 100)) + geom_point(size = 2/fig_factor) + geom_line(lwd = 1/fig_factor)+ xlab("Pseudotime Upper Bound n") + ylab("% Total Diff. Exp Genes") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + ylim(0, 100)
rm(lineplot_temp)

###Fig 01e
#RNA Velocity plot of differentially expressed genes across the in vivo trajectory. Note - this had to be made in Python, so please see the corresponding velo_code.py.

###Fig 01f
#Velocity lengths over pseudotime for the in vivo trajectory
ggplot(invivo_velocity, aes(x = individual_pseudotime * max(cds_invivo@colData$pseudotime), y = velocity_length)) + geom_point(aes(fill = timepoint), colour="black", pch=21, size = 4/fig_factor) + xlab("Pseudotime") + ylab("Velocity Length") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1/fig_factor, 'lines')) + guides(fill = guide_legend(title = "Timepoint", override.aes = list(size = 4/fig_factor)))

###Fig 01g
#Workflow figure for computing time to x% FC for an individual gene trajectory. This figure is mostly for example - I selected Cox5a as it appeared to be visually appealing for demonstrating the procedure.
norm_data = as.data.frame(as.matrix(t(normalized_counts(cds_invivo, norm_method = "size_only"))))
norm_data$pseudotime = cds_invivo@colData$pseudotime
norm_data$timepoint = cds_invivo@colData$timepoint
ggplot(norm_data, aes(x = pseudotime, y = Cox5a)) + geom_point(aes(color = timepoint), alpha = 0.4) + geom_line(data = data.frame(pseudotime = seq(0, 60), Cox5a = as.numeric(model_expectation_invivo["Cox5a", ])), size = 1.5) + ylim(0, 60) + theme_linedraw() + theme(legend.position = "none", axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), axis.text.y = element_text(size = 0)) + xlab("Pseudotime") + ylab("Gene Expression")
rm(norm_data)

###Fig 01h
#Smooth curves for the identified downregulated and upregulated clusters
#Downregulated
ggplot(median_curves[startsWith(as.character(median_curves$Var2), "down"), ], aes(x = Var1, y = value, color = sapply(strsplit(as.character(Var2), " "), "[[", 2))) + geom_smooth(size = 1.5/fig_factor, linetype = "longdash") + xlab("Pseudotime") + ylab("Scaled Expression") + theme_classic() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), axis.text.y = element_text(size = 0), axis.ticks.y = element_blank(), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), panel.background = element_rect(color = "black", size = 1/fig_factor), legend.justification=c(1,1), legend.position=c(0.98,0.98), legend.key.size = unit(1/fig_factor, 'lines')) + labs(color = "Gene Cluster")
#Upregulated
ggplot(median_curves[startsWith(as.character(median_curves$Var2), "up"), ], aes(x = Var1, y = value, color = sapply(strsplit(as.character(Var2), " "), "[[", 2))) + geom_smooth(size = 1.5/fig_factor, linetype = "longdash") + xlab("Pseudotime") + ylab("Scaled Expression") + theme_classic() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), axis.text.y = element_text(size = 0), axis.ticks.y = element_blank(), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), panel.background = element_rect(color = "black", size = 1/fig_factor), legend.justification=c(1,1), legend.position=c(0.28,0.98), legend.key.size = unit(1/fig_factor, 'lines')) + labs(color = "Gene Cluster")

###Fig 01i
#Parameters used to cluster the genes - the time to 10%, 50%, and 95% - plotted per cluster
ggplot(invivo_genes_melt, aes(x = as.factor(cluster_name), y = value, fill = cluster_name)) + geom_boxplot() + geom_jitter(alpha = 0.03) + facet_wrap(~prettylabel, scale = 'free_y', labeller = labeller(sample_type = supp.labs)) + xlab("Cluster") + ylab("Pseudotime") + theme_bw() + theme(axis.title = element_text(size = 18/fig_factor), axis.text = element_text(size = 10/fig_factor), legend.position = "none", strip.text = element_text(size = 16/fig_factor))

###Fig 01j
#Chart showing number of up and downregulated genes in each cluster. The chart was made manually, but the values can be found by running this:
table(invivo_genes$cluster_name, invivo_genes$direction)

###Fig01k
#Heatmap of top 25 TFs enriched in the various clusters in vivo.
pheatmap(t(clus_agg[2:26, ]), show_rownames = TRUE, show_colnames = TRUE, cluster_rows = FALSE, cluster_cols = FALSE, color =  rev(colorRampPalette(brewer.pal(9, "RdBu"))(100)), cellwidth = 10, cellheight = 10)

###Supp 01a
#Heatmap of the SingleCellNet scores for the in vivo and in vitro cells. Unfortunately, because of the way my code written, it is difficult to independently call this figure. It is the output of the data_qc function and can be recovered by running that code.

###Supp 01b
#Gene ontology for up and downregulated genes over in vivo pseudotime. The following parameters were used to visualize: Genes were input to Gene Ontology website with default settings, and results downloaded. Top 150 terms by fold enrichment were loaded into Revigo with corresponding q-value. Condensed list was set to small.

###Supp 01c
#Bar plot of the proportions of spliced to unspliced counts for the in vivo cells per timepoint. Note - for convenience, this was made in Python, so please see the corresponding velo_code.py.

###Supp 02a-d - Kinetic heatmap and GO/Pathway terms for each gene cluster. This code is modular so just replace the first lines with the appropriate cluster and the cluster/annotation in the kinetic heatmap
#Settings from WebGestalt - BP terms, size limit 850, top 25 terms, weighted set cover expecting 10

go_terms_down = read.csv("GO_new/clus45_down.txt", as.is = TRUE, header=  TRUE, sep = "\t")
go_terms_down_sub = read.csv("GO_new/clus45_down_sub.txt", as.is = TRUE, header=  TRUE, sep = "\t")
go_terms_down = go_terms_down[go_terms_down$geneSet %in% go_terms_down_sub[, 1], ]
go_terms_down$full = str_wrap(go_terms_down$description, 30)
go_terms_up = read.csv("GO_new/clus45_up.txt", as.is = TRUE, header=  TRUE, sep = "\t")
go_terms_up_sub = read.csv("GO_new/clus45_up_sub.txt", as.is = TRUE, header=  TRUE, sep = "\t")
go_terms_up = go_terms_up[go_terms_up$geneSet %in% go_terms_up_sub[, 1], ]
go_terms_up$full = str_wrap(go_terms_up$description, 30)
g1 = ggplot(m_melt[m_melt$cluster_name == "EA3" & m_melt$direction == "down", ], aes(x = Var2, y = value)) + geom_density_2d(lwd = 0.5/fig_factor, color = "#00BF7D") + theme_linedraw() + xlab("Pseudotime") + ylab("Scaled Expression") + theme(axis.title = element_text(size = 18/fig_factor), axis.text = element_text(size = 14/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + annotate("text", x=Inf, y = Inf, label = "812 genes", vjust=2, hjust=1.1, size = 7/fig_factor) #hjust 1.1 for downregulated
g2 = ggplot(m_melt[m_melt$cluster_name == "EA3" & m_melt$direction == "up", ], aes(x = Var2, y = value)) + geom_density_2d(lwd = 0.5/fig_factor, color = "#00BF7D") + theme_linedraw() + xlab("Pseudotime") + ylab("Scaled Expression") + theme(axis.title = element_text(size = 18/fig_factor), axis.text = element_text(size = 14/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + annotate("text", x=Inf, y = Inf, label = "124 genes", vjust=2, hjust=2.05, size = 7/fig_factor) #hjust 2.05-2.35 for upregulated
g3 = ggplot(go_terms_down, aes(x = reorder(description, enrichmentRatio), y = enrichmentRatio)) + geom_bar(stat = "identity", alpha = 0.1, fill = "blue") + theme_classic() +  theme(axis.title.y = element_blank(), axis.title = element_text(size = 18/fig_factor), axis.text = element_text(size = 14/fig_factor), axis.text.y  = element_blank(), panel.background = element_rect(color = "black", size = 1/fig_factor)) + scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) + geom_text(aes(y = 0.1 ,label = full), hjust = "left", size = 4.75/fig_factor, lineheight = 2/fig_factor) + ylab("Enrichment Ratio") + coord_flip()
g4 = ggplot(go_terms_up, aes(x = reorder(description, enrichmentRatio), y = enrichmentRatio)) + geom_bar(stat = "identity", alpha = 0.1, fill = "red") + theme_classic() +  theme(axis.title.y = element_blank(), axis.title = element_text(size = 18/fig_factor), axis.text = element_text(size = 14/fig_factor), axis.text.y  = element_blank(), panel.background = element_rect(color = "black", size = 1/fig_factor)) + scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) + geom_text(aes(y = 0.1 ,label = full), hjust = "left", size = 4.75/fig_factor, lineheight = 2/fig_factor) + ylab("Enrichment Ratio") + coord_flip()
grid.arrange(g2, g1, g4, g3, layout_matrix = rbind(c(1, 2), c(1, 2), c(3, 4), c(3,4), c(3,4)))#Use for figures

#Note that for clusters 4+5, since we combined them together for the sake of GO term analysis, the kinetic maps have to be made slightly differently. The code for g1 + g2 can be found here.
g1 = ggplot() + geom_density_2d(data = m_melt[m_melt$cluster_name == "LA1" & m_melt$direction == "down", ], aes(x = Var2, y = value), lwd = 0.5/fig_factor, color = "#00B0F6") + geom_density_2d(data = m_melt[m_melt$cluster_name == "LA2" & m_melt$direction == "down", ], aes(x = Var2, y = value), lwd = 0.5/fig_factor, color = "#E76BF3") + theme_linedraw() + xlab("Pseudotime") + ylab("Scaled Expression") + theme(axis.title = element_text(size = 18/fig_factor), axis.text = element_text(size = 14/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + annotate("text", x=Inf, y = Inf, label = "321 genes\n77 genes", vjust=1.2, hjust=1.1, size = 7/fig_factor) 
g2 = ggplot() + geom_density_2d(data = m_melt[m_melt$cluster_name == "LA1" & m_melt$direction == "up", ], aes(x = Var2, y = value), lwd = 0.5/fig_factor, color = "#00B0F6") + geom_density_2d(data = m_melt[m_melt$cluster_name == "LA2" & m_melt$direction == "up", ], aes(x = Var2, y = value), lwd = 0.5/fig_factor, color = "#E76BF3") + theme_linedraw() + xlab("Pseudotime") + ylab("Scaled Expression") + theme(axis.title = element_text(size = 18/fig_factor), axis.text = element_text(size = 14/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + annotate("text", x=Inf, y = Inf, label = "170 genes\n14 genes", vjust=1.2, hjust=2.05, size = 7/fig_factor) #hjust 2.05-2.35 for upregulated

rm(g1, g2, g3, g4, go_terms_down, go_terms_down_sub, go_terms_up, go_terms_up_sub)

#####Fig 02/Supp 03
#The following section outlines all of the panels associated with Figure 2, which details the in vitro trajectory reconstruction and alignments with the invivo datasets. This also includes Supplementary Figure 3, which provides additional details about the in vitro trajectory.

###Fig 02a
#This is a workflow image and thus there is nothing to plot.

###Fig 02b
#This is an immunofluorescence image and thus there is nothing to plot.

###Fig 02c
#Reconstruction of the in vitro trajectory using supervised feature selection
plot_cells(cds_invitro, label_groups_by_cluster = FALSE, label_cell_groups = FALSE, color_cells_by = "timepoint", cell_size = 2/fig_factor, show_trajectory_graph = FALSE) + xlab("Component 1") + ylab("Component 2") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1/fig_factor, 'lines')) + guides(color = guide_legend(title = "Timepoint", override.aes = list(size = 4/fig_factor)))

###Fig 02d
#Boxplot of pseudotime scores assigned by Monocle to in vitro trajectory by timepoint
ggplot(as.data.frame(cds_invitro@colData), aes(x = timepoint, y = pseudotime, fill = timepoint)) + geom_boxplot(lwd = 1/fig_factor, outlier.size = 0) + geom_jitter(size = 1.5/fig_factor, alpha = 0.3) + coord_flip() + xlab("Timepoint") + ylab("Pseudotime") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.position = "none")

###Fig 02e 
#Gene ontology for differentially regulated genes in vitro over pseudotime. The following parameters were used to visualize: Genes were input to Gene Ontology website with default settings, and results downloaded. Top 150 terms by fold enrichment were loaded into Revigo with corresponding q-value. Condensed list was set to small.

###Fig 02f
#Combined trajectory (using appropriate batch correction) by timepoint
plot_cells(cds_combined, label_groups_by_cluster = FALSE, label_cell_groups = FALSE, color_cells_by = "timepoint", cell_size = 1.5/fig_factor, show_trajectory_graph = FALSE) + xlab("Component 1") + ylab("Component 2") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1/fig_factor, 'lines')) + guides(color = guide_legend(title = "Timepoint", override.aes = list(size = 4), nrow = 15)) + scale_color_manual(values = c(colorRampPalette(brewer.pal(9, "YlOrRd"))(15), brewer.pal(9, "Blues")))

###Fig 02g
#Combined trajectory, treating group as a batch effect, by timepoint
plot_cells(cds_single2, label_groups_by_cluster = FALSE, label_cell_groups = FALSE, color_cells_by = "timepoint", cell_size = 1.5/fig_factor, show_trajectory_graph = FALSE) + xlab("Component 1") + ylab("Component 2") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1/fig_factor, 'lines')) + guides(color = guide_legend(title = "Timepoint", override.aes = list(size = 4), nrow = 15)) + scale_color_manual(values = c(colorRampPalette(brewer.pal(9, "YlOrRd"))(15), brewer.pal(9, "Blues")))

###Fig 02h
#Trajectory boxplot by timepoint for combined trajectory treating group as a batch effect
ggplot(as.data.frame(cds_single2@colData), aes(x = timepoint, y = pseudotime, fill = timepoint)) + geom_boxplot(lwd = 1/fig_factor, outlier.size = 0) + geom_jitter(size = 1.5/fig_factor, alpha = 0.3) + coord_flip() + xlab("Timepoint") + ylab("Pseudotime") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.position = "none") + scale_fill_manual(values = c(colorRampPalette(brewer.pal(9, "YlOrRd"))(15), brewer.pal(9, "Blues")))

###Fig 02i
#Combined trajectory, treating in vivo only batch as reference batch, by timepoint
plot_cells(cds_single, label_groups_by_cluster = FALSE, label_cell_groups = FALSE, color_cells_by = "timepoint", cell_size = 1.5/fig_factor, show_trajectory_graph = FALSE) + xlab("Component 1") + ylab("Component 2") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1/fig_factor, 'lines')) + guides(color = guide_legend(title = "Timepoint", override.aes = list(size = 4), nrow = 15)) + scale_color_manual(values = c(colorRampPalette(brewer.pal(9, "YlOrRd"))(15), brewer.pal(9, "Blues")))

###Fig 02j
#Trajectory boxplot by timepoint for combined trajectory treating in vivo only batch as reference batch
ggplot(as.data.frame(cds_single@colData), aes(x = timepoint, y = pseudotime, fill = timepoint)) + geom_boxplot(lwd = 1/fig_factor, outlier.size = 0) + geom_jitter(size = 1.5/fig_factor, alpha = 0.3) + coord_flip() + xlab("Timepoint") + ylab("Pseudotime") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.position = "none") + scale_fill_manual(values = c(colorRampPalette(brewer.pal(9, "YlOrRd"))(15), brewer.pal(9, "Blues")))

###Fig 02k
#Entropy scores across in vivo and in vitro cells
ggplot(as.data.frame(cds_combined@colData), aes(x = timepoint, y = entropy, fill = timepoint)) + geom_boxplot(lwd = 1/fig_factor, outlier.size = 0) + geom_jitter(size = 1/fig_factor, alpha = 0.3) + xlab("Timepoint") + ylab("Entropy") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), axis.text.x = element_text(angle = 45), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.position = "none") + scale_fill_manual(values = c(colorRampPalette(brewer.pal(9, "YlOrRd"))(15), brewer.pal(9, "Blues")))

###Supp 03a
#This is an flow cytometry image and thus there is nothing to plot.

###Supp 03b
#Boxplot of pseudotime scores assigned by Monocle to in vitro trajectory by timepoint, but further labeled by line
ggplot(as.data.frame(cds_invitro@colData), aes(x = timepoint, y = pseudotime, fill = line)) + geom_boxplot(lwd = 1/fig_factor, outlier.size = 0) + geom_jitter(size = 1.5/fig_factor) + coord_flip() + xlab("Timepoint") + ylab("Pseudotime") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + guides(fill = guide_legend(title = "Line", override.aes = list(size = 1/fig_factor)))

###Supp 03c
#Boxplot of pseudotime scores assigned by Monocle to in vitro trajectory by timepoint, but futher labeled by size
ggplot(as.data.frame(cds_invitro@colData[cds_invitro@colData$timepoint %in% c("D25", "D30", "D45"), ]), aes(x = timepoint, y = pseudotime, fill = size)) + geom_boxplot(lwd = 1/fig_factor, outlier.size = 0) + geom_jitter(size = 1.5/fig_factor) + coord_flip() + xlab("Timepoint") + ylab("Pseudotime") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + guides(fill = guide_legend(title = "Size", override.aes = list(size = 1/fig_factor)))

###Supp 03d
#RNA Velocity plot of differentially expressed genes across the in vitro trajectory. Note - this had to be made in Python, so please see the corresponding velo_code.py.

###Supp 03e
#Velocity lengths over pseudotime for the in vitro trajectory
ggplot(invitro_velocity, aes(x = individual_pseudotime * max(cds_invitro@colData$pseudotime), y = velocity_length)) + geom_point(aes(fill = timepoint), colour="black", pch=21, size = 4/fig_factor) + xlab("Pseudotime") + ylab("Velocity Length") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1/fig_factor, 'lines')) + guides(fill = guide_legend(title = "Timepoint", override.aes = list(size = 4/fig_factor)))

#####Fig 03/Supp 04
#The following section outlines all of the panels associated with Figure 3, which details the dysregulation of perinatal maturation programs in PSC-CMs. It also includes Supplementary Figure 4, which contains additional details about the comparison of in vivo and in vitro datasets

###Fig 03a
#Identification of global dysregulation of perinatal maturation genes. This is largely a schematic, but the values presented can be found by:
table(pseudo_ids_invivo_perinatal %in% global_ids)
table(pseudo_ids_invivo_perinatal %in% global_ids_early)

###Fig 03b
#Heatmap for various categories of globally dysregulated genes (sarcomeric, mitochondrial, ribosomal)
temp_cols = c("red3", "blue")
names(temp_cols) = c("In Vivo", "In Vitro")
pheatmap(global_data_numeric[c("Myh6", "Myl2", "Myl3", "Myl4", "Myl7", "Ttn", "Myom1", "Myom2", "Tnnt2", "Tnni1", "Tnni3", "Actn2", "mt-Co1", "mt-Co2", "mt-Nd1", "mt-Nd4", "Ndufa8", "Ndufa10", "Cox8a", "Cox17", "Atp5g2", "Atp5a1", "Atp5b", "Atp2a2", global_ids[!global_ids %in% global_ids_clean][1:12]), ], show_rownames = TRUE, show_colnames = TRUE, cluster_rows = FALSE, cluster_cols= FALSE, scale = "row", annotation = data.frame(row.names = colnames(global_data_numeric), group = c(rep("In Vivo", 12), rep("In Vitro", 7))), annotation_colors = list(group = temp_cols), annotation_row = data.frame(row.names = c("Myh6", "Myl2", "Myl3", "Myl4", "Myl7", "Ttn", "Myom1", "Myom2", "Tnnt2", "Tnni1", "Tnni3", "Actn2", "mt-Co1", "mt-Co2", "mt-Nd1", "mt-Nd4", "Ndufa8", "Ndufa10", "Cox8a", "Cox17", "Atp5g2", "Atp5a1", "Atp5b", "Atp2a2", global_ids[!global_ids %in% global_ids_clean][1:12]), category = c(rep("Sarcomeric", 12), rep("Mitochondrial", 12), rep("Ribosomal", 12))), color =  rev(colorRampPalette(brewer.pal(9, "RdBu"))(100)), fontsize = 10/fig_factor)

###Fig 03c
#Histogram of mean in vitro expression levels of perinatal maturation genes after scaling by min and max levels in vivo.
ggplot(global_vals, aes(x = invitro_level)) + geom_point(aes(y = 0), size = 1.5/fig_factor) + geom_line(data = data.frame(x = c(0, 0), y = c(0, 0.1)), aes(x = x, y = y)) + geom_line(data = data.frame(x = c(1, 1), y = c(0, 0.1)), aes(x = x, y = y)) + theme_classic() + geom_density(lwd = 2/fig_factor)  +  coord_cartesian(clip = 'off') + theme(axis.title = element_text(size = 0), axis.text = element_text(size = 18/fig_factor)) + scale_y_continuous(limits = c(0,1.1), expand = c(0, 0)) + scale_x_continuous(breaks=seq(-3,4,1), limits = c(-3, 4))
#The value presented in the chart can be found by:
table(global_vals$invitro_level < 0 | global_vals$invitro_level > 1)

###Fig 03d
#Compared fold changes between in vitro and in vivo cells for two times in the in vivo trajectory - pseudotime 20 and 40, corresponding to critical periods in the progression of maturation.

break_at_2 = function(limits) {
  seq(ceiling(limits[1]), floor(limits[2]), 2)
}
#In vivo pseudotime = 20
ggplot(compared_fc[subset_diff[[20]], ], aes(x = invivo_fc_20, y = invitro_fc, color = diff_invitro)) + geom_hline(yintercept = 0, alpha = 0.5) + geom_vline(xintercept = 0, alpha = 0.5) +  geom_point(data = compared_fc[compared_fc$diff_invitro == FALSE, ], alpha = 0.3) + geom_point(data = compared_fc[compared_fc$diff_invitro == TRUE, ], alpha = 0.4) + coord_fixed() + scale_x_continuous(breaks = break_at_2) + scale_y_continuous(breaks = break_at_2) + xlab(expression(log[2](FC) ~ "In Vivo")) + ylab(expression(log[2](FC) ~ "In Vitro")) + theme_linedraw() + scale_color_manual(name = "Diff. Expressed \nIn Vitro", values = c("black", "red")) + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor))
#In vivo pseudotime = 40
ggplot(compared_fc[subset_diff[[40]], ], aes(x = invivo_fc_40, y = invitro_fc, color = diff_invitro)) + geom_hline(yintercept = 0, alpha = 0.5) + geom_vline(xintercept = 0, alpha = 0.5) +  geom_point(data = compared_fc[compared_fc$diff_invitro == FALSE, ], alpha = 0.3) + geom_point(data = compared_fc[compared_fc$diff_invitro == TRUE, ], alpha = 0.4) + coord_fixed() + scale_x_continuous(breaks = break_at_2) + scale_y_continuous(breaks = break_at_2) + xlab(expression(log[2](FC) ~ "In Vivo")) + ylab(expression(log[2](FC) ~ "In Vitro")) + theme_linedraw() + scale_color_manual(name = "Diff. Expressed \nIn Vitro", values = c("black", "red")) + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor))

###Fig 03e
#Venn diagram of differentially expressed genes between in vivo and in vitro cells
#Upregulated genes
v1 =venn.diagram(x = list(rownames(invivo_genes)[invivo_genes$direction == "up"], names(pseudo_ids_invitro_fc)[pseudo_ids_invitro_fc > 0]), category.names = c("In Vivo", "In Vitro"), filename = NULL, output = TRUE, col = "black",
                 lty = "dotted",
                 lwd = 4/fig_factor,
                 fill = c("red", "cornflowerblue"),
                 alpha = 0.50,
                 label.col = c("darkred", "darkorchid4", "darkblue"),
                 cex = 1.5/fig_factor,
                 fontfamily = "sans",
                 fontface = "bold",
                 cat.col = c("darkred", "darkblue"),
                 cat.cex = 1.5/fig_factor,
                 cat.fontfamily = "sans", margin = 0.05, cat.dist = 0.05, cat.pos = c(-30, 30))
grid.newpage()
grid.draw(v1)
rm(v1)
#Downregulated genes
v1 =venn.diagram(x = list(rownames(invivo_genes)[invivo_genes$direction == "down"], names(pseudo_ids_invitro_fc)[pseudo_ids_invitro_fc < 0]), category.names = c("In Vivo", "In Vitro"), filename = NULL, output = TRUE, col = "black",
                 lty = "dotted",
                 lwd = 4/fig_factor,
                 fill = c("red", "cornflowerblue"),
                 alpha = 0.50,
                 label.col = c("darkred", "darkorchid4", "darkblue"),
                 cex = 1.5/fig_factor,
                 fontfamily = "sans",
                 fontface = "bold",
                 cat.col = c("darkred", "darkblue"),
                 cat.cex = 1.5/fig_factor,
                 cat.fontfamily = "sans", margin = 0.05, cat.dist = 0.05, cat.pos = c(-30, 30))
grid.newpage()
grid.draw(v1)
rm(v1)

###Fig 03f
#Lineplot showing the % of in vivo genes that are correctly differentially regulated in vitro at each pseudotime bin (for mouse data). Note that in a sense, we are showing *cumulative* percentages.
lineplot_temp = data.frame(pseudotime = seq(1, 42), kannan = unlist(lapply(subset_diff[1:42], function(x){sum(pseudo_ids_invitro_good %in% x)/length(x)})))
ggplot(lineplot_temp, aes(x = as.numeric(pseudotime), y = kannan * 100)) + geom_point(size = 2/fig_factor) + geom_line(lwd = 1/fig_factor)+ xlab("In Vivo Pseudotime Upper Bound n") + ylab("% Genes Recapitulated in PSC-CMs") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + ylim(0, 60)
rm(lineplot_temp)

###Fig 03g
#Barplot showing % of in vivo genes in each cluster that are correctly differentially regulated in vitro (for mouse data).
bar_temp = as.data.frame(table(invivo_genes[pseudo_ids_invitro_good, ]$cluster_name)/table(invivo_genes$cluster_name))
ggplot(bar_temp, aes(x = Var1, y = Freq * 100, fill = Var1)) + geom_bar(stat = "identity", alpha = 0.7) + geom_text(aes(label = paste(round(Freq * 100, 2), "%", sep ="")), size = 8/fig_factor) + xlab("In Vivo Gene Cluster") + ylab("% Genes Recapitulated in PSC-CMs") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.position = "none")
rm(bar_temp)

###Fig 03h
#Lineplot showing the % of in vivo genes that are correctly differentially regulated in vitro at each pseudotime bin (for human data). As above, we are showing *cumulative* percentages.
lineplot_temp = data.frame(pseudotime = seq(1, 42), churko = unlist(lapply(subset_diff[1:42], function(x){sum(human_correct_mouse[["cds_churko"]] %in% x)/length(x[x %in% pseudo_ids_invivo_ortho])})), friedman = unlist(lapply(subset_diff[1:42], function(x){sum(human_correct_mouse[["cds_friedman"]] %in% x)/length(x[x %in% pseudo_ids_invivo_ortho])})), gerbin = unlist(lapply(subset_diff[1:42], function(x){sum(human_correct_mouse[["cds_gerbin"]] %in% x)/length(x[x %in% pseudo_ids_invivo_ortho])})), ruan = unlist(lapply(subset_diff[1:42], function(x){sum(human_correct_mouse[["cds_ruan"]] %in% x)/length(x[x %in% pseudo_ids_invivo_ortho])})))
lineplot_temp = melt(lineplot_temp, id.vars = "pseudotime")
ggplot(lineplot_temp, aes(x = as.numeric(pseudotime), y = value*100, color = variable)) + geom_point(size = 2/fig_factor) + geom_line(lwd = 1/fig_factor) + xlab("In Vivo Pseudotime Upper Bound n") + ylab("% Genes Recapitulated in PSC-CMs") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.position = "none") + scale_color_discrete(name = "Study", labels = c("C", "FNL", "GGDH", "RL")) + ylim(0, 82)
rm(lineplot_temp)

###Fig 03i
#Barplot showing % of in vivo genes in each cluster that are correctly differentially regulated in vitro (for human data).
bar_temp = rbind(as.data.frame(table(invivo_genes[human_correct_mouse[["cds_churko"]], ]$cluster_name)/table(invivo_genes$cluster_name)), as.data.frame(table(invivo_genes[human_correct_mouse[["cds_friedman"]], ]$cluster_name)/table(invivo_genes$cluster_name)), as.data.frame(table(invivo_genes[human_correct_mouse[["cds_gerbin"]], ]$cluster_name)/table(invivo_genes$cluster_name)), as.data.frame(table(invivo_genes[human_correct_mouse[["cds_ruan"]], ]$cluster_name)/table(invivo_genes$cluster_name)))
bar_temp$study = c(rep("C", 5), rep("FNL", 5), rep("GGDH", 5), rep("RL", 5))
ggplot(bar_temp, aes(x = study, y = Freq * 100, fill = Var1)) + geom_bar(stat = "identity", alpha = 0.7, position = "dodge") + xlab("Study") + ylab("% Genes Recapitulated in PSC-CMs") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor))
rm(bar_temp)

###Supp 04a 
#Gene ontology for globally differentially upregulated and downregulated genes between in vivo and in vitro CMs. The following parameters were used to visualize: Genes were input to Gene Ontology website with default settings, and results downloaded. Top 150 terms by fold enrichment were loaded into Revigo with corresponding q-value. Condensed list was set to small.

###Supp 04b
#UMAP reduced dimensionality for each of the human in vitro PSC-CM datasets. We constructed these in a semi-supervised approach (though similar results are seen with the unbiased approach) - we found that rather than forming clear trajectories, the timepoints typically completely separated out.
#Churko
plot_cells(cds_churko, label_groups_by_cluster = FALSE, label_cell_groups = FALSE, color_cells_by = "timepoint", cell_size = 1/fig_factor, show_trajectory_graph = FALSE) + xlab("Component 1") + ylab("Component 2") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1/fig_factor, 'lines')) + guides(color = guide_legend(title = "Timepoint", override.aes = list(size = 4/fig_factor)))
#Friedman, Nguyen, and Lukowski
plot_cells(cds_friedman, label_groups_by_cluster = FALSE, label_cell_groups = FALSE, color_cells_by = "timepoint", cell_size = 1/fig_factor, show_trajectory_graph = FALSE) + xlab("Component 1") + ylab("Component 2") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1/fig_factor, 'lines')) + guides(color = guide_legend(title = "Timepoint", override.aes = list(size = 4/fig_factor)))
#Gerbin, Grancharova, Donovan-Maiye, and Hendershott
plot_cells(cds_gerbin, label_groups_by_cluster = FALSE, label_cell_groups = FALSE, color_cells_by = "timepoint", cell_size = 1/fig_factor, show_trajectory_graph = FALSE) + xlab("Component 1") + ylab("Component 2") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1/fig_factor, 'lines')) + guides(color = guide_legend(title = "Timepoint", override.aes = list(size = 4/fig_factor)))
#Ruan and Liao
plot_cells(cds_ruan, label_groups_by_cluster = FALSE, label_cell_groups = FALSE, color_cells_by = "timepoint", cell_size = 1/fig_factor, show_trajectory_graph = FALSE) + xlab("Component 1") + ylab("Component 2") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1/fig_factor, 'lines')) + guides(color = guide_legend(title = "Timepoint", override.aes = list(size = 4/fig_factor)))

###Supp 04c
#Time to 10%, 50%, and 95% plots for the correctly differentially expressed and dysregulated genes in mouse and human datasets. For simplicity, with the human datasets we used the consensus differentially expressed and dysregulated gene lists. To be correctly differentially expressed, the gene had to be correctly differentially expressed in >= 2 studies; for the dysregulation, it had to be dysregulated >= 3 studies.
#Mouse
ggplot(invivo_genes_melt[!invivo_genes_melt$mouse == "n/a", ], aes(x = as.factor(mouse), y = value)) + geom_violin() + geom_boxplot(width = 0.4, outlier.size = 0)+ geom_jitter(alpha = 0.02) + facet_wrap(~prettylabel, scale = 'free_y', labeller = labeller(sample_type = supp.labs)) + xlab("Gene Group") + ylab("Pseudotime") + theme_bw() + theme(axis.title = element_text(size = 18/fig_factor), axis.text = element_text(size = 10/fig_factor), legend.position = "none", strip.text = element_text(size = 16/fig_factor))
#Human
ggplot(invivo_genes_melt[!invivo_genes_melt$human == "n/a", ], aes(x = as.factor(human), y = value)) + geom_violin() + geom_boxplot(width = 0.4, outlier.size = 0)+ geom_jitter(alpha = 0.02) + facet_wrap(~prettylabel, scale = 'free_y', labeller = labeller(sample_type = supp.labs)) + xlab("Gene Group") + ylab("Pseudotime") + theme_bw() + theme(axis.title = element_text(size = 18/fig_factor), axis.text = element_text(size = 10/fig_factor), legend.position = "none", strip.text = element_text(size = 16/fig_factor))

###Supp 04d
#Percentage of genes in correctly differentially expression and dysregulated gene lists with accessible chromatin (based on ATAC-seq peaks) in the promoter-TSS regions.

ggplot(combined_atac, aes(x= variable, y = value * 100, group = sample, color = study)) + geom_point(size = 2/fig_factor) + geom_line(size = 0.5/fig_factor) + ylim(0, 80) + xlab("Gene Group") + ylab("% Genes with TSS-Promoter Peak") + theme_linedraw() + theme(axis.title = element_text(size = 20/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 18/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1/fig_factor, 'lines')) + scale_color_discrete(name = "Study")

###Supp 04e
#STRING-DB network graphs of the dysregulated maturation TFs. This figure was made by simply inputting the names of the TFs into STRING and then downloading the resulting network graph (thus, there is nothing to plot here).

#####Fig 04
#The following section outlines all of the panels associated with Figure 4, which details the consensus dysregulated gene list and upstream analysis.

###Fig 04a
#Workflow and results figure for identifying the consensus dysregulated gene list. There is nothing to plot in R for this figure. Numbers can be found simply by looking at lengths of corresponding gene lists.

###Fig 04b
#Lineplot showing the % of in vivo genes that are dysregulated in vitro at each pseudotime bin (for human data). As above, we are showing *cumulative* percentages.
lineplot_temp = data.frame(pseudotime = seq(1, 42), kannan = unlist(lapply(subset_diff[1:42], function(x){sum(dysreg%in% x)/length(x)})))
ggplot(lineplot_temp, aes(x = as.numeric(pseudotime), y = kannan * 100)) + geom_point(size = 2/fig_factor) + geom_line(lwd = 1/fig_factor)+ xlab("In Vivo Pseudotime Upper Bound n") + ylab("% Genes Dysregulated in PSC-CMs") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + ylim(0, 22)
rm(lineplot_temp)

###Fig 04c
#Barplot showing % of in vivo genes in each cluster that are dysregulated.
bar_temp = as.data.frame(table(invivo_genes[dysreg, ]$cluster_name)/table(invivo_genes$cluster_name))
ggplot(bar_temp, aes(x = Var1, y = Freq * 100, fill = Var1)) + geom_bar(stat = "identity", alpha = 0.7) + geom_text(aes(label = paste(round(Freq * 100, 2), "%", sep ="")), size = 8/fig_factor) + xlab("In Vivo Gene Cluster") + ylab("% Genes Dysregulated in PSC-CMs") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.position = "none")
rm(bar_temp)

###Fig 04d
#Results of TF target enrichment for the dysregulated gene list. This analysis was done in WebGestalt, using genome as the reference, and selecting all TFs with B-H p-val < 0.05. We then used affinity propagation to cluster and select most interesting TFs.
ggplot(dysreg_tf_good, aes(x = reorder(clean_name, enrichmentRatio), y = enrichmentRatio)) + geom_bar(stat = "identity") + coord_flip() + ylab("Fold Enrichment") + xlab("Enriched TFs") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor)) + scale_y_continuous(limits = c(0,2.5), expand = c(0, 0))

###Fig 04e
#Heatmap of expression values across pseudotime in vivo and in vitro for identified TFS of interest. For the sake of visualization, we limited the n vivo pseudotimes to the perinatal periods.
temp_cols = c("red3", "blue")
names(temp_cols) = c("In Vivo", "In Vitro")
pheatmap(global_data_numeric[c("Srf", "Yy1", "Jun", "Nfe2l2", "Sox9", "Ppara", "Esrra", "Nrf1", "Mef2a"), c(1:9, 13:19)], show_rownames = TRUE, show_colnames = TRUE, cluster_rows = FALSE, cluster_cols= FALSE, scale = "row", annotation = data.frame(row.names = colnames(global_data_numeric)[c(1:9, 13:19)], group = c(rep("In Vivo", 9), rep("In Vitro", 7))), annotation_colors = list(group = temp_cols), color =  rev(colorRampPalette(brewer.pal(9, "RdBu"))(100)), fontsize = 10/fig_factor, cellwidth = 10, cellheight = 10)

###Fig 04f
#IPA Activation Z-scores for TFs of interest. Note that NRF1 didn't have detected activity. We initially performed several comparisons in IPA, but we decided that the most biologically relevant were the in vivo activation across the entire perinatal period (e.g. [0, 40]) and the activation across in vitro using the genes differentially expressed in vitro.
colnames(dysreg_ipa) = c("In Vivo 0-20", "In Vivo 20-40", "In Vivo", "In Vitro", "In Vitro 0-40")
pheatmap(dysreg_ipa[c("SRF", "YY1", "JUN", "NFE2L2", "PPARA", "Esrra", "NRF1", "MEF2A"), c(3,4)], show_rownames = TRUE, show_colnames = TRUE, cluster_rows = FALSE, cluster_cols= FALSE, scale = "none", color =  c(rep("#2166AC", 4), rev(colorRampPalette(brewer.pal(9, "RdBu"))(43)), "#B2182B"), breaks = c(seq(-6, -3), seq(-2, 2, length.out = 42), 3), fontsize = 10/fig_factor, cellwidth = 10, cellheight = 10)

#####Fig 05
#The following section outlines all of the panels associated with Figure 5, which focuses on analysis of maturation pathways in various tissue engineered datasets drawn from the literature.

###Fig 05a
#Workflow figure showing all of the perturbations studied as well as some relevant experimental details. This is a schematic so there is nothing to plot.

###Fig 05b
#Table showing the number of differentially regulated genes from each dataset that are in the perinatal maturation gene list and either differentially regulated in the same or opposite direction. There is nothing to plot here; the numbers are listed in the corresponding section in the code.R file.

###Fig 05c
#Barplot showing, of the genes that are differentially regulated in the same direction as in vivo in each tissue engineered dataset, what percentage can be found in either the consensus good human gene list and what percentage can be found in the consensus dysregulated gene list.
te_list = list(branco_good, zhao_good, kupp_let_good, kupp_age_good, feyen_good, hcts_good, hcas_good, giacomelli_good)
bar_temp = data.frame(study = c("3D Agg.", "ES", "Let7 OE", "LTC", "MM", "EMT", "AIS", "CC"), dysreg = unlist(lapply(te_list, function(x){table(x %in% dysreg)[2]/length(x)})), good = unlist(lapply(te_list, function(x){table(x %in% human_good)[2]/length(x)})))
bar_temp = melt(bar_temp, id.vars = "study")
bar_temp$pretty_label = "Correct"
bar_temp[bar_temp$variable == "dysreg", ]$pretty_label = "Dysregulated"
bar_temp$study = factor(bar_temp$study, levels = c("3D Agg.", "ES", "Let7 OE", "LTC", "MM", "EMT", "AIS", "CC"))
ggplot(bar_temp, aes(x = study, y = value * 100, fill = pretty_label)) + geom_bar(stat = "identity", alpha = 0.7, position = "dodge") + xlab("Study") + ylab("% of Correctly Diff. Exp. Genes") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.justification=c(1,1), legend.position = c(1,1)) + guides(fill = guide_legend(title = "Gene Group", override.aes = list(size = 3/fig_factor)))
rm(bar_temp, te_list)

###Fig 05d
#Lineplot showing the % of in vivo genes that are differentially regulated in the same direction in each tissue engineered dataset at each pseudotime.
lineplot_temp = data.frame(pseudotime = seq(1, 42), branco = unlist(lapply(subset_diff[1:42], function(x){sum(branco_good %in% x)/length(x[x %in% pseudo_ids_invivo_ortho])})), zhao = unlist(lapply(subset_diff[1:42], function(x){sum(zhao_good %in% x)/length(x[x %in% pseudo_ids_invivo_ortho])})), kupp_let = unlist(lapply(subset_diff[1:42], function(x){sum(kupp_let_good %in% x)/length(x[x %in% pseudo_ids_invivo_ortho])})), kupp_age = unlist(lapply(subset_diff[1:42], function(x){sum(kupp_age_good %in% x)/length(x[x %in% pseudo_ids_invivo_ortho])})), feyen = unlist(lapply(subset_diff[1:42], function(x){sum(feyen_good %in% x)/length(x[x %in% pseudo_ids_invivo_ortho])})), lam_mt = unlist(lapply(subset_diff[1:42], function(x){sum(hcts_good %in% x)/length(x[x %in% pseudo_ids_invivo_ortho])})), lam_sheet = unlist(lapply(subset_diff[1:42], function(x){sum(hcas_good %in% x)/length(x[x %in% pseudo_ids_invivo_ortho])})), giacomelli = unlist(lapply(subset_diff[1:42], function(x){sum(giacomelli_good %in% x)/length(x[x %in% pseudo_ids_invivo_ortho])})))
lineplot_temp = melt(lineplot_temp, id.vars = "pseudotime")
ggplot(lineplot_temp, aes(x = as.numeric(pseudotime), y = value*100, color = variable)) + geom_point(size = 2/fig_factor) + geom_line(lwd = 1/fig_factor) + xlab("In Vivo Pseudotime Upper Bound n") + ylab("% Genes Correctly Diff. Exp.") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.justification = c(1,1), legend.position = c(1,1)) + scale_color_discrete(name = "Study", labels = c(c("3D Agg.", "ES", "Let7 OE", "LTC", "MM", "EMT", "AIS", "CC")))
rm(lineplot_temp)

###Fig 05e
#Barplot showing % of in vivo genes in each cluster that are differentially regulated in the same direction in each perturbation study.
bar_temp = rbind(as.data.frame(table(invivo_genes[branco_good, ]$cluster_name)/table(invivo_genes[pseudo_ids_invivo_ortho, ]$cluster_name)), as.data.frame(table(invivo_genes[zhao_good, ]$cluster_name)/table(invivo_genes[pseudo_ids_invivo_ortho, ]$cluster_name)), as.data.frame(table(invivo_genes[kupp_let_good, ]$cluster_name)/table(invivo_genes[pseudo_ids_invivo_ortho, ]$cluster_name)), as.data.frame(table(invivo_genes[kupp_age_good, ]$cluster_name)/table(invivo_genes[pseudo_ids_invivo_ortho, ]$cluster_name)), as.data.frame(table(invivo_genes[feyen_good, ]$cluster_name)/table(invivo_genes[pseudo_ids_invivo_ortho, ]$cluster_name)), as.data.frame(table(invivo_genes[hcts_good, ]$cluster_name)/table(invivo_genes[pseudo_ids_invivo_ortho, ]$cluster_name)), as.data.frame(table(invivo_genes[hcas_good, ]$cluster_name)/table(invivo_genes[pseudo_ids_invivo_ortho, ]$cluster_name)), as.data.frame(table(invivo_genes[giacomelli_good, ]$cluster_name)/table(invivo_genes[pseudo_ids_invivo_ortho, ]$cluster_name)))
bar_temp$study = c(rep("3D Agg.", 5), rep("ES", 5), rep("Let7 OE", 5), rep("LTC", 5), rep("MM", 5), rep("EMT", 5), rep("AIS", 5), rep("CC", 5))
bar_temp$study = factor(bar_temp$study, levels = c("3D Agg.", "ES", "Let7 OE", "LTC", "MM", "EMT", "AIS", "CC"))
ggplot(bar_temp, aes(x = study, y = Freq * 100, fill = Var1)) + geom_bar(stat = "identity", alpha = 0.7, position = "dodge") + xlab("Study") + ylab("% Genes Correctly Diff. Exp.") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor))
rm(bar_temp)

###Fig 05f
#Heatmap showing whether each dysregulated TF is upregulated in each of the perturbation datasets.
te_up_list = list(branco_up_mouse, zhao_up_mouse, kupp_let_up_mouse, kupp_age_up_mouse, feyen_up_mouse, hcts_ids_up_mouse, hcas_ids_up_mouse, giacomelli_ids_up_mouse)
heatmap_temp = as.data.frame(lapply(te_up_list, function(x){c("Srf", "Yy1", "Jun", "Nfe2l2", "Sox9", "Ppara", "Esrra", "Nrf1", "Mef2a") %in% x}))
rownames(heatmap_temp) = c("Srf", "Yy1", "Jun", "Nfe2l2", "Sox9", "Ppara", "Esrra", "Nrf1", "Mef2a")
colnames(heatmap_temp) = c("3D Agg.", "ES", "Let7 OE", "LTC", "MM", "EMT", "AIS", "CC")
pheatmap(data.matrix(heatmap_temp), show_rownames = TRUE, show_colnames = TRUE, cluster_rows = FALSE, cluster_cols= FALSE, scale = "none", color =  rev(colorRampPalette(brewer.pal(9, "RdBu"))(10)), fontsize = 10/fig_factor, cellwidth = 10, cellheight = 10)
rm(te_up_list, heatmap_temp)

###Fig 05g
#Heatmap showing enrichment ratios for downstream targets of each dysregulated TF in each of the perturbation datasets.
pheatmap(te_tf, show_rownames = TRUE, show_colnames = TRUE, cluster_rows = FALSE, cluster_cols= FALSE, scale = "none", color =  c(rev(colorRampPalette(brewer.pal(9, "RdBu"))(100)), rep("#B2182B", 4)), breaks = c(seq(0, 4, length.out = 100), seq(5, 8)), fontsize = 10/fig_factor, cellwidth = 10, cellheight = 10)

#####Supp Figure 05
#The following section outlines all of the panels associated with Suppemental Figure 5, which deals with batch effects and technical issues in trajectory reconstruction.

###Supp 05a
#In vivo trajectory without any batch correction, coloured by batch and timepoint.
cds_test <- new_cell_data_set(data[, pheno_data$group == "in vivo" & pheno_data$top5_norm < 1.8 & pheno_data$depth_norm > -0.7], cell_metadata = pheno_data[pheno_data$group == "in vivo" & pheno_data$top5_norm < 1.8 & pheno_data$depth_norm > -0.7, ], gene_metadata = G_list[rownames(data), ])
cds_test <- preprocess_cds(cds_test, num_dim = 5)
cds_test <- align_cds(cds_test)
cds_test <- reduce_dimension(cds_test)
cds_test@colData$batch = as.factor(cds_test@colData$batch)
plot_cells(cds_test, label_groups_by_cluster = FALSE, label_cell_groups = FALSE, color_cells_by = "batch", cell_size = 1/fig_factor, show_trajectory_graph = FALSE) + xlab("Component 1") + ylab("Component 2") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1/fig_factor, 'lines')) + guides(color = guide_legend(title = "Batch", override.aes = list(size = 4/fig_factor)))
plot_cells(cds_test, label_groups_by_cluster = FALSE, label_cell_groups = FALSE, color_cells_by = "timepoint", cell_size = 1/fig_factor, show_trajectory_graph = FALSE) + xlab("Component 1") + ylab("Component 2") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1/fig_factor, 'lines')) + guides(color = guide_legend(title = "Timepoint", override.aes = list(size = 4/fig_factor)))
rm(cds_test)

###Supp 05b
#Boxplot of pseudotime scores assigned by Monocle to in vivo trajectory by timepoint, but further labeled by batch
ggplot(as.data.frame(cds_invivo@colData), aes(x = timepoint, y = pseudotime, fill = as.factor(batch))) + geom_boxplot(lwd = 1/fig_factor, outlier.size = 0) + geom_jitter(size = 1.5/fig_factor, alpha = 0.3) + coord_flip() + xlab("Timepoint") + ylab("Pseudotime") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + scale_fill_discrete(name = "Batch")

###Supp 05c
#Unbiased reconstruction of the in vitro trajectory
plot_cells(cds_invitro_unbiased, label_groups_by_cluster = FALSE, label_cell_groups = FALSE, color_cells_by = "timepoint", cell_size = 2/fig_factor, show_trajectory_graph = FALSE) + xlab("Component 1") + ylab("Component 2") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1/fig_factor, 'lines')) + guides(color = guide_legend(title = "Timepoint", override.aes = list(size = 4/fig_factor)))

###Supp 05d
#Combined trajectory, with in vivo cells from Batch 3 specifically labeled.
cds_combined@colData$batch3_invivo = (cds_combined@colData$batch == 3 & cds_combined@colData$group == "in vivo")
plot_cells(cds_combined, label_groups_by_cluster = FALSE, label_cell_groups = FALSE, color_cells_by = "batch3_invivo", cell_size = 2/fig_factor, show_trajectory_graph = FALSE) + xlab("Component 1") + ylab("Component 2") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1/fig_factor, 'lines')) + guides(color = guide_legend(title = "In Vivo Batch 3", override.aes = list(size = 4/fig_factor)))

#####Leftover bits and pieces of code

#We spiked our cells with ERCC spike-ins during the RT step of library preparation. Unfortunately, our spike-in percentages were generally too low to be used confidently for normalization, and at any rate it is likely that each batch received slightly different amounts of spike-in. However, they could be used as a broad approximate for RNA quantity. We focused on batch 3 since this batch contained both in vivo and in vitro cells. For visualization, we also filtered out extreme ERCC expression (>500 counts). We plotted ERCC counts as well total depth divided by ERCC counts (a pseudo-normalization).

ggplot(pheno_data[pheno_data$batch == 3 & rownames(pheno_data) %in% colnames(cds_combined) & pheno_data$ercc < 500, ], aes(x = timepoint, y = ercc, fill = group)) + geom_violin(scale = "width", position = position_dodge(width = 1), lwd =  1/fig_factor) + geom_boxplot(position = position_dodge(width = 1), lwd = 1/fig_factor, outlier.size = 1.5/fig_factor) + xlab("Timepoint") + ylab("ERCC UMI Counts") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + scale_fill_discrete(name = "Group")
ggplot(pheno_data[pheno_data$batch == 3 & rownames(pheno_data) %in% colnames(cds_combined) & pheno_data$ercc < 500, ], aes(x = timepoint, y = log10(depth/ercc), fill = group)) + geom_violin(scale = "width", position = position_dodge(width = 1), lwd =  1/fig_factor) + geom_boxplot(position = position_dodge(width = 1), lwd = 1/fig_factor, outlier.size = 1.5/fig_factor) + xlab("Timepoint") + ylab(expression(log[10](Depth / ERCCs))) + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + scale_fill_discrete(name = "Group")

lineplot_temp = data.frame(pseudotime = seq(1, 60), kannan = unlist(lapply(subset_diff2, function(x){sum(pseudo_ids_invitro %in% x)/length(x)})), churko = unlist(lapply(subset_diff2, function(x){sum(pseudo_ids_churko %in% x)/length(x)})), friedman = unlist(lapply(subset_diff2, function(x){sum(pseudo_ids_friedman %in% x)/length(x)})), gerbin = unlist(lapply(subset_diff2, function(x){sum(pseudo_ids_gerbin %in% x)/length(x)})), ruan = unlist(lapply(subset_diff2, function(x){sum(pseudo_ids_ruan %in% x)/length(x)})))
lineplot_temp = melt(lineplot_temp, id.vars = "pseudotime")
ggplot(lineplot_temp, aes(x = pseudotime, y = value, color = variable)) + geom_point() + geom_line()

lineplot_temp = data.frame(pseudotime = seq(5, 60, 5), kannan = unlist(lapply(seq(0,11), function(x){sum(pseudo_ids_invitro %in% unlist(subset_diff2[((5*x)+1):(5*(x+1))]))/length(unlist(subset_diff2[((5*x)+1):(5*(x+1))]))})), churko = unlist(lapply(seq(0,11), function(x){sum(pseudo_ids_churko %in% unlist(subset_diff2[((5*x)+1):(5*(x+1))]))/length(unlist(subset_diff2[((5*x)+1):(5*(x+1))]))})), friedman = unlist(lapply(seq(0,11), function(x){sum(pseudo_ids_friedman %in% unlist(subset_diff2[((5*x)+1):(5*(x+1))]))/length(unlist(subset_diff2[((5*x)+1):(5*(x+1))]))})), gerbin = unlist(lapply(seq(0,11), function(x){sum(pseudo_ids_gerbin %in% unlist(subset_diff2[((5*x)+1):(5*(x+1))]))/length(unlist(subset_diff2[((5*x)+1):(5*(x+1))]))})), ruan = unlist(lapply(seq(0,11), function(x){sum(pseudo_ids_ruan %in% unlist(subset_diff2[((5*x)+1):(5*(x+1))]))/length(unlist(subset_diff2[((5*x)+1):(5*(x+1))]))})))

###Supp 02d
#Aligned principal components for the unbiased in vitro trajectory reconstruction. We use the aligned coordinates here so that there is no issue of potential confounding from the library batches. (We focus on PC1-2)
plot_cells(cds_invitro_unbiased, label_groups_by_cluster = FALSE, label_cell_groups = FALSE, color_cells_by = "timepoint", cell_size = 2/fig_factor, show_trajectory_graph = FALSE, reduction = "Aligned", x = 1, y = 2) + xlab("Component 1") + ylab("Component 2") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1/fig_factor, 'lines')) + guides(color = guide_legend(title = "Timepoint", override.aes = list(size = 4/fig_factor)))

geneList = list(branco_good, zhao_good, kupp_let_good, kupp_age_good, feyen_good, hcas_good, hcts_good, giacomelli_good)
genedf <- stack(setNames(geneList, nm=c("3DAgg", "Stim.", "Let7", "LTC", "MM", "Sheet", "MT", "Co-EC")))
table(genedf[2:1]) %*% t(table(genedf[2:1]))
