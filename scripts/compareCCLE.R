# Loads CCLE data and maps the TCPA to the CCLE data.
# Comparative analusis of TCPA (RPPA) and CCLE (gene expression from microarrays) data.

rm(list=ls())

library(CePa)
library(data.table)
library(grid)
library(dendextend)  # for tanglegram
library(gplots)
library(RColorBrewer)
library(tsne)
library(hexbin)
library(class)
# library(SDMTools)  # for confusion.matrix
library(mda)


setwd("~/Google Drive/projects/tcpa-rppa")
data_dir = "~/DataProjects/tcpa-rppa"

source("lib/annotate.r")  # annotation mapping functions
source("lib/colors.r")  # tissue color

# Load data
# ------------------------------------------------------------
ccle = fread(file.path(data_dir, "data/CCLE_Expression_Entrez_2012-09-29.gct"))
ccle = as.data.frame(ccle)  # convert from data.table to data.frame

# Sample descriptions
ccle_info = read.table(file.path(data_dir, "data/CCLE_sample_info_file_2012-10-18.txt"),
	header=TRUE, sep="\t")

# load TCPA data
tcpa = read.table(file.path(data_dir, "data/MCLP-RBN-v1.0-whole_set.tsv"),
	header=TRUE, sep="\t")

x = tcpa[,4:ncol(tcpa)]
rownames(x) = tcpa$Sample_Name
colnames(x) = colnames(tcpa)[4:ncol(tcpa)]

prot1 = "ECADHERIN"

prot2 = "BETACATENIN"
# prot2 = "AXL"
# prot2 = "NCADHERIN"
# prot2 = "CLAUDIN7"

i = which(colnames(x) == prot1)
j = which(colnames(x) == prot2)

plot(x[,i], x[,j], main=cor(x[,i], x[,j], use="pairwise.complete.obs"), xlab=prot1, ylab=prot2)

# plot(x[,i], x[,j], xlab=prot1, ylab=prot2)



# match the CCLE data names (in the columns) to the sample data.
# -------------------------------------------------------------

# align CCLE sample info to the CCLE main data
ccle_info_id = match(colnames(ccle), ccle_info$CCLE.name)
ccle_info_align = ccle_info[ccle_info_id,]

# get cell line names of all CCLE data.
ccle_cell_lines = ccle_info_align$Cell.line.primary.name[ccle_info_id]

# match TCPA cell line name to CCLE cell line names.
# find corresponding CCLE for each TCPA cell line.
ccle_match_id = match(stemName(tcpa$Sample_Name), stemName(ccle_info_align$Cell.line.primary.name))

# fraction found statistic
sum(!is.na(ccle_match_id)) / length(ccle_match_id)

# Get the matching CCLE data matrix
ccle_match = ccle[,ccle_match_id[!is.na(ccle_match_id)]]
ccle_match_info = ccle_info_align[ccle_match_id[!is.na(ccle_match_id)],]  # sample information

# Get the subset of the TCPA data which matches CCLE
tcpa_match = tcpa[!is.na(ccle_match_id),]

# check match
data.frame(tcpa=tcpa_match$Sample_Name, ccle=ccle_match_info$Cell.line.primary.name)



# COSMIC sample data
cosmic = fread(file.path(data_dir, "data/CosmicSample.tsv"))
table(cosmic$tumour_source)


cosmic_match_id = match(stemName(tcpa_match$Sample_Name), stemName(cosmic$sample_name))

sum(!is.na(cosmic_match_id)) / length(cosmic_match_id)

cosmic_match = cosmic[cosmic_match_id, ]

# cosmic_match$tumour_source
# table(cosmic_match$tumour_source)
# table(cosmic_match$metastatic_site)
# table(cosmic_match$age_at_tumour_recurrence)
# table(cosmic_match$stage)


# TCPA

# Format TCPA data matrix (separate data from sample information).
tcpa_match_mat = tcpa_match[,4:ncol(tcpa_match)]
tcpa_match_mat = as.matrix(tcpa_match_mat)
rownames(tcpa_match_mat) = tcpa_match[,1]  #

# filter by  observations
min_cell_lines = 100  # the minimum number of cell lines measured for an included protein

# count the number of non missing data entries for a particular protein.
protein_counts = apply(tcpa_match_mat, 2, function(col) {
	return(sum(!is.na(col)))
})

tcpa_match_filter_mat  = tcpa_match_mat[,protein_counts > min_cell_lines]


# Analysis of the matching TCPA and CCLE data.
# ------------------------------------------------

# colors = c(brewer.pal(12, "Paired"), brewer.pal(12, "Set3"), brewer.pal(9, "Set1"), brewer.pal(8, "Set2"))
# colors = c(brewer.pal(8, "Dark2"), brewer.pal(12, "Set3"), brewer.pal(9, "Set1"), brewer.pal(8, "Set2"))

# Distance matrix calculation
tcpa_dist = dist(tcpa_match_filter_mat)
ccle_dist = dist(t(ccle_match))

# Combined distance measure.
# weighted distance, with equal weight to each dataset
comb_dist = ccle_dist + sum(ccle_dist) / sum(tcpa_dist) * tcpa_dist



# t-SNE analysis of both TCPA data and matching CCLE data.
# uses the distance matrices tcpa_dist and ccle_dist
# -------------------------------------------------------------------------------

max_iter = 5000

# tcpa_tsne_embed = tsne(x, perplexity=40)
tcpa_tsne_embed = tsne(tcpa_dist, max_iter=max_iter)

# ccle_match_info$Site.Primary
# table(ccle_match_info$Site.Primary)
# unique(ccle_match_info$Site.Primary)

# CCLE 
ccle_tsne_embed = tsne(ccle_dist, max_iter=max_iter)

# Combined and scaled distance
comb_tsne_embed = tsne(comb_dist, max_iter=max_iter)







# Calculate z-scores for proteins
ccle_match_zscore = t(scale(t(ccle_match), center=TRUE, scale=TRUE))
tcpa_match_filter_zscore = scale(tcpa_match_filter_mat)





# Find a gene symbol in the CCLE data
rids = grep("MAP2K1", ccle$Description)
ccle$Description[rids]
# plot(tcpa_match_filter_zscore[,tcpa_n], ccle_match_zscore[ccle_n,], col=colors)



svg("exploratory/figCompareCCLE/tsne_comparison3.svg", width=11, height=8.5)

par(mfrow=(c(3, 4)))

# Primary Site
plot(tcpa_tsne_embed[,1], tcpa_tsne_embed[,2],
	main="RPPA",
	col=colors[match(ccle_match_info$Site.Primary, unique(ccle_match_info$Site.Primary))],
	xlab="", ylab="",
	pch=16)

plot(ccle_tsne_embed[,1], ccle_tsne_embed[,2],
	main="mRNA",
	xlab="", ylab="",
	col=colors[match(ccle_match_info$Site.Primary, unique(ccle_match_info$Site.Primary))], pch=16)

plot(comb_tsne_embed[,1], comb_tsne_embed[,2],
	main="RPPA + mRNA",
	xlab="", ylab="",
	col=colors[match(ccle_match_info$Site.Primary, unique(ccle_match_info$Site.Primary))], pch=16)

plot(1, 1, main="Primary site", type="n",
	xlab="", ylab="")
legend("topleft", legend=unique(ccle_match_info$Site.Primary),
	col=colors,
	pch=16,
	cex=0.75)

protein_name = "ECADHERIN"
ccle_n = 18925  # ECADHERIN feature in CCLE data. Found by manual lookup.
tcpa_n = match(protein_name, colnames(tcpa_match_filter_mat))

pal = rev(brewer.pal(9, "RdYlBu"))
colors_ccle = colorGradient(ccle_match_zscore[ccle_n,], pal)
colors_tcpa = colorGradient(tcpa_match_filter_zscore[,tcpa_n], pal)
colors_comb = colorGradient((tcpa_match_filter_zscore[,tcpa_n] + ccle_match_zscore[ccle_n,]) / 2, pal)

plot(tcpa_tsne_embed[,1], tcpa_tsne_embed[,2],
	main=paste("RPPA", protein_name),
	col= colors_tcpa,
	xlab="", ylab="",
	pch=16)

plot(ccle_tsne_embed[,1], ccle_tsne_embed[,2],
	main="mRNA E-Cadherin",
	col=colors_ccle,
	xlab="", ylab="",
	pch=16
)

plot(comb_tsne_embed[,1], comb_tsne_embed[,2],
	main=paste("RPPA+mRNA", protein_name),
	col=colors_comb,
	xlab="", ylab="",
	pch=16
)

plot(0, 0, type="n")  # void, for layout

site_colors = brewer.pal(9, "Set1")
site_colors = c(site_colors[1], "grey", site_colors[2:5])
site_col = site_colors[as.integer(factor(cosmic_match$tumour_source))]

plot(tcpa_tsne_embed[,1], tcpa_tsne_embed[,2],
	main="RPPA",
	col=site_col,
	pch=16)

plot(ccle_tsne_embed[,1], ccle_tsne_embed[,2],
	main="mRNA",
	col=site_col,
	pch=16)

plot(comb_tsne_embed[,1], comb_tsne_embed[,2],
	main="mRNA",
	col=site_col,
	pch=16)

legend("bottomleft", legend=levels(factor(cosmic_match$tumour_source)),
	col=site_colors,
	pch=16,
	cex=0.75)

# ccle_n = grep("MAP2K1", ccle$Description)  # MEK1
# tcpa_n = match("MEK1_pS217S221", colnames(tcpa_match_filter_mat))

# colors_ccle = colorGradient(ccle_match_zscore[ccle_n,], pal)
# colors_tcpa = colorGradient(tcpa_match_filter_zscore[,tcpa_n], pal)

# plot(tcpa_tsne_embed[,1], tcpa_tsne_embed[,2],
# 	main="RPPA MEK1_pS217S221",
# 	col= colors_tcpa,
# 	xlab="", ylab="",
# 	pch=16)

# plot(ccle_tsne_embed[,1], ccle_tsne_embed[,2],
# 	main="mRNA, MEK1",
# 	col=colors_ccle,
# 	xlab="", ylab="",
# 	pch=16
# )

dev.off()



# E-Cadherin expression and metastasis?
# t-tests
# -------------------------------------------
ccle_n = 18925  # ECADHERIN
tcpa_n = match("ECADHERIN", colnames(tcpa_match_filter_mat))
tcpa_n = match("SNAIL", colnames(tcpa_match_filter_mat))

metastatic = cosmic_match$tumour_source == "metastasis"
primary = cosmic_match$tumour_source == "primary"

# ccle_match_zscore[ccle_n,]

# boxplot(
# 	list(
# 		metastatic=tcpa_match_filter_zscore[metastatic, tcpa_n],
# 		primary=tcpa_match_filter_zscore[primary, tcpa_n]
# 	)
# )

# t.test(
# 	tcpa_match_filter_zscore[metastatic, tcpa_n],
# 	tcpa_match_filter_zscore[primary, tcpa_n]
# )


# t.test(
# 	ccle_match_zscore[metastatic, tcpa_n],
# 	ccle_match_zscore[primary, tcpa_n]
# )


t_tests = lapply(1:ncol(tcpa_match_filter_zscore), function(i) {
	t.test(
		tcpa_match_filter_zscore[metastatic, i],
		tcpa_match_filter_zscore[primary, i]
	)
})
names(t_tests) = colnames(tcpa_match_filter_zscore)

p_vals = sapply(t_tests, function(x) x$p.value)

estimates = sapply(t_tests, function(x) x$estimate)

barplot(estimates[, p_vals < 0.05])

tcpa_match_filter_zscore[, p_vals < 0.05]

sort(p_vals)




# Bootstrap tSNE with different perplexity parameters
# ------------------------------------------------------------------------
perp_vals = seq(from=5, to=50, by=5)

tcpa_tsne_bootstrap = lapply(perp_vals, function(perplexity) {
	message(perplexity)
	embed = tsne(tcpa_dist, perplexity=perplexity, max_iter=1000)
	return(embed)
})

ccle_tsne_bootstrap = lapply(perp_vals, function(perplexity) {
	message(perplexity)
	embed = tsne(ccle_dist, perplexity=perplexity, max_iter=1000)
	return(embed)
})

# par(mfrow=c(length(perp_vals), 2))

svg("exploratory/figCompareCCLE/perplexity_bootstrap.svg", width=6, height=13)
par(mfrow=c(10, 4), mar=c(1, 1, 1, 1))
protein_name = "ECADHERIN"
ccle_n = 18925  # ECADHERIN feature in CCLE data. Found by manual lookup.
tcpa_n = match(protein_name, colnames(tcpa_match_filter_mat))

for (k in 1:length(perp_vals)) {
	plot(tcpa_tsne_bootstrap[[k]][,1], tcpa_tsne_bootstrap[[k]][,2],
		main=paste0("RPPA p=", perp_vals[k]),
		col=colors[match(ccle_match_info$Site.Primary, unique(ccle_match_info$Site.Primary))],
		xlab="", ylab="",
		cex=0.8, cex.main=0.7,
		xaxt="n", yaxt="n",
		pch=16)

	pal = rev(brewer.pal(9, "RdYlBu"))
	colors_ccle = colorGradient(ccle_match_zscore[ccle_n,], pal)
	colors_tcpa = colorGradient(tcpa_match_filter_zscore[,tcpa_n], pal)
	# colors_comb = colorGradient((tcpa_match_filter_zscore[,tcpa_n] + ccle_match_zscore[ccle_n,]) / 2, pal)

	plot(tcpa_tsne_bootstrap[[k]][,1], tcpa_tsne_bootstrap[[k]][,2],
		main=paste0("RPPA ", protein_name, " p=", perp_vals[k]),
		col= colors_tcpa,
		xlab="", ylab="",
		cex=0.8, cex.main=0.7, 
		xaxt="n", yaxt="n",
		pch=16)


	plot(ccle_tsne_bootstrap[[k]][,1], ccle_tsne_bootstrap[[k]][,2],
		main=paste0("CCLE", " p=", perp_vals[k]),
		col=colors[match(ccle_match_info$Site.Primary, unique(ccle_match_info$Site.Primary))],
		xlab="", ylab="",
		cex=0.8, cex.main=0.7,
		xaxt="n", yaxt="n",
		pch=16)

	plot(ccle_tsne_bootstrap[[k]][,1], ccle_tsne_bootstrap[[k]][,2],
		main=paste0("CCLE ", protein_name, " p=", perp_vals[k]),
		col= colors_ccle,
		xlab="", ylab="",
		cex=0.8, cex.main=0.7,
		xaxt="n", yaxt="n",
		pch=16)
}
dev.off()


# Principal component analysis (PCA)
# --------------------------------------

# library(mice)
library(pcaMethods)


tcpa_pca = pca(tcpa_match_filter_mat, method="svdImpute", center=TRUE, scale="uv")
ccle_pca = pca(t(ccle_match), method="svdImpute", center=TRUE, scale="uv")


protein_name = "ECADHERIN"
ccle_n = 18925  # ECADHERIN feature in CCLE data. Found by manual lookup.
tcpa_n = match(protein_name, colnames(tcpa_match_filter_mat))

pal = rev(brewer.pal(9, "RdYlBu"))
colors_ccle = colorGradient(ccle_match_zscore[ccle_n,], pal)
colors_tcpa = colorGradient(tcpa_match_filter_zscore[,tcpa_n], pal)


svg("exploratory/figCompareCCLE/PCA_RPPA_mRNA.svg")
par(mfrow=c(2, 2))
plot(scores(tcpa_pca)[, 1], scores(tcpa_pca)[, 2],
	col=colors[match(ccle_match_info$Site.Primary, unique(ccle_match_info$Site.Primary))],
	xlab="PC1", ylab="PC2",
	main="TCPA, tissue of origin",
	pch=16)

plot(scores(tcpa_pca)[, 1], scores(tcpa_pca)[, 2],
	col=colors_tcpa,
	xlab="PC1", ylab="PC2",
	main="TCPA, E-cadherin expresison",
	pch=16)

plot(scores(ccle_pca)[, 1], scores(ccle_pca)[, 2],
	col=colors[match(ccle_match_info$Site.Primary, unique(ccle_match_info$Site.Primary))],
	xlab="PC1", ylab="PC2",
	main="mRNA, tissue of origin",
	pch=16)

plot(scores(ccle_pca)[, 1], scores(ccle_pca)[, 2],
	col=colors_ccle,
	xlab="PC1", ylab="PC2",
	main="mRNA, E-cadherin expression",
	pch=16)
dev.off()


# tcpa_pca = prcomp(tcpa_match_filter_mat)
# tcpa_pca = prcomp(na.omit(tcpa_match_filter_mat), na.action=na.omit)
# tcpa_pca = prcomp(t(tcpa_match_filter_mat), na.action=na.omit)

# t(ccle_match)




protein_name = "ECADHERIN"
# protein_name = "CAVEOLIN1"
protein_name = "MEK1"
tcpa_n = match(protein_name, colnames(tcpa_match_filter_mat))
colors_tcpa = colorGradient(tcpa_match_filter_zscore[,tcpa_n], pal)

par(mfrow=c(2, 1))
plot(tcpa_tsne_embed[,1], tcpa_tsne_embed[,2],
	main="RPPA",
	col=colors[match(ccle_match_info$Site.Primary, unique(ccle_match_info$Site.Primary))],
	pch=16)
plot(tcpa_tsne_embed[,1], tcpa_tsne_embed[,2],
	main=paste("RPPA", protein_name),
	col= colors_tcpa,
	pch=16)



# Old, comparison plot
pdf("exploratory/tsne_comparison.pdf", width=16, height=1)
par(mfrow=(c(3, 4)))

# Primary Site
plot(tcpa_tsne_embed[,1], tcpa_tsne_embed[,2],
	main="RPPA",
	col=colors[match(ccle_match_info$Site.Primary, unique(ccle_match_info$Site.Primary))],
	pch=16)

plot(ccle_tsne_embed[,1], ccle_tsne_embed[,2],
	main="mRNA",
	col=colors[match(ccle_match_info$Site.Primary, unique(ccle_match_info$Site.Primary))], pch=16)

plot(comb_tsne_embed[,1], comb_tsne_embed[,2],
	main="RPPA + mRNA",
	col=colors[match(ccle_match_info$Site.Primary, unique(ccle_match_info$Site.Primary))], pch=16)

plot(1, 1, main="Primary site", type="n")
legend("topleft", legend=unique(ccle_match_info$Site.Primary),
	col=colors,
	pch=16,
	cex=0.75)

# Histology
plot(tcpa_tsne_embed[,1], tcpa_tsne_embed[,2],
	main="RPPA",
	col=colors[match(ccle_match_info$Histology, unique(ccle_match_info$Histology))], pch=16)

plot(ccle_tsne_embed[,1], ccle_tsne_embed[,2],
	main="mRNA",
	col=colors[match(ccle_match_info$Histology, unique(ccle_match_info$Histology))], pch=16)

plot(comb_tsne_embed[,1], comb_tsne_embed[,2],
	main="RPPA + mRNA",
	col=colors[match(ccle_match_info$Histology, unique(ccle_match_info$Histology))], pch=16)

plot(1, 1, main="Histology", type="n")
legend("topleft", legend=unique(ccle_match_info$Histology),
	col=colors,
	pch=16,
	cex=0.75)

# Hist subtype

common_hist_subtypes = names(sort(table(ccle_match_info$Hist.Subtype1), decreasing=TRUE)[1:length(colors)])
hist_col = colors[match(ccle_match_info$Hist.Subtype1, common_hist_subtypes)]
hist_col[is.na(hist_col)] = "grey"

plot(tcpa_tsne_embed[,1], tcpa_tsne_embed[,2],
	main="RPPA",
	# col=colors[match(ccle_match_info$Hist.Subtype1, unique(ccle_match_info$Hist.Subtype1))],
	col=hist_col,
	pch=16)

plot(ccle_tsne_embed[,1], ccle_tsne_embed[,2],
	main="mRNA",
	# col=colors[match(ccle_match_info$Hist.Subtype1, unique(ccle_match_info$Hist.Subtype1))],
	col=hist_col,
	pch=16)

plot(comb_tsne_embed[,1], comb_tsne_embed[,2],
	main="RPPA + mRNA",
	# col=colors[match(ccle_match_info$Hist.Subtype1, unique(ccle_match_info$Hist.Subtype1))],
	col=hist_col,
	pch=16)

plot(1, 1, main="26 most common hist. subtypes", type="n")
legend("topleft", legend=common_hist_subtypes,
	col=colors,
	pch=16,
	cex=0.75)





hist_col = colors[match(ccle_match_info$Hist.Subtype1, common_hist_subtypes)]

# # Gender
# plot(tcpa_tsne_embed[,1], tcpa_tsne_embed[,2],
# 	main="MCLP",
# 	col=colors[match(ccle_match_info$Gender, unique(ccle_match_info$Gender))], pch=16)

# plot(ccle_tsne_embed[,1], ccle_tsne_embed[,2],
# 	main="CCLE",
# 	col=colors[match(ccle_match_info$Gender, unique(ccle_match_info$Gender))], pch=16)

# plot(1, 1, main="Gender", type="n")
# legend("topleft", legend=unique(ccle_match_info$Gender),
# 	col=colors,
# 	pch=16,
# 	cex=0.75)
dev.off()


# Similarity analysis
# -------------------------------------------------------------------------------------
# x = tcpa_match_mat


tcpa_match_norm = apply(tcpa_match_filter_mat, 1, scale, scale=TRUE, center=TRUE)
ccle_match_norm = t(apply(ccle_match, 1, scale, scale=TRUE, center=TRUE))

# tcpa_match_norm = apply(tcpa_match_mat, 1, scale, scale=TRUE, center=TRUE)
# ccle_match_norm = t(apply(ccle_match, 1, scale, scale=TRUE, center=TRUE))


# check
# ccle_match_norm = as.matrix(ccle_match_norm)

# tcpa_cor = cor(tcpa_match_norm, method="spearman", use="pairwise.complete.obs")
# ccle_cor = cor(ccle_match_norm, method="spearman")

# same analysis based on distance measures instead of correlation
# tcpa_cor = dist(tcpa_match)
# ccle_cor = dist(ccle_match_norm)


# Format dist object to distance matrices
tcpa_dist_mat = as.matrix(tcpa_dist)
ccle_dist_mat = as.matrix(ccle_dist)

# # Overwrite CCLE cell line names with TCPA names (as matched previously)
# warning("overwrite CCLE cell line names:")
# # check overwritten names
# overwrite_id = which(rownames(ccle_dist_mat) != rownames(tcpa_dist_mat))
# data.frame(from=rownames(ccle_dist_mat)[overwrite_id], to=rownames(tcpa_dist_mat)[overwrite_id])

# Overwrite names
rownames(ccle_dist_mat) = rownames(tcpa_dist_mat)
colnames(ccle_dist_mat) = colnames(tcpa_dist_mat)


# Regression model of distances
dist_regr = lm(as.vector(tcpa_dist_mat) ~ as.vector(ccle_dist_mat))

# Get matrix (cell x cell) of residuals of linear regression
residuals = matrix(resid(dist_regr), nrow=nrow(tcpa_dist_mat), ncol=ncol(tcpa_dist_mat))
rownames(residuals) = rownames(tcpa_dist_mat)
colnames(residuals) = colnames(tcpa_dist_mat)

# Flatten residuals next to pairwise distances for plotting
x = data.frame(tcpa=as.vector(as.matrix(tcpa_dist_mat)), ccle=as.vector(as.matrix(ccle_dist_mat)), res=as.vector(residuals))
x = x[order(x$res, decreasing=T),]

# pdf("exploratory/top_distance_residuals.pdf")
# k = 500
# plot(c(50, 250), c(5, 35), type="n",
# 	main=paste("Top", k, "residuals"),
# 	xlab="mRNA distance", ylab="RPPA distance")
# points(x$ccle[1:k], x$tcpa[1:k], col=color_scale[length(color_scale)])
# points(x$ccle[(nrow(x)-k+1):nrow(x)], x$tcpa[(nrow(x)-k+1):nrow(x)], col=color_scale[1])
# abline(dist_regr)
# dev.off()

# hexbin plot of all pairwise distances with the top and bottom residuals highligted as points.
# -------------------------------------------------------------------------
pdf("exploratory/pairwise_dist.pdf")
k = 200  # top and bottom k pairwise distances
p = plot(hexbin(ccle_dist_mat[lower.tri(ccle_dist_mat)], tcpa_dist_mat[lower.tri(tcpa_dist_mat)]),
	xlab="mRNA distance",
	ylab="RPPA distance",
	main=paste0("Pairwise cell line distances, cor=", round(cor(tcpa_dist_mat[lower.tri(tcpa_dist_mat)], ccle_dist_mat[lower.tri(ccle_dist_mat)], method="spearman"), 3))
	)

# push plot viewport
pushHexport(p$plot.vp)

# Plot top and bottom
grid.points(x$ccle[1:k], x$tcpa[1:k], gp=gpar(cex=0.5, col=color_scale[length(color_scale)]))
grid.points(x$ccle[(nrow(x)-k+1):nrow(x)], x$tcpa[(nrow(x)-k+1):nrow(x)], gp=gpar(cex=0.5, col=color_scale[1]))

# Plot regression line
grid.abline(dist_regr$coefficients[1], dist_regr$coefficients[2], gp=gpar(col=rgb(44, 162, 95, maxColorValue=255)))
dev.off()


# Heatmaps of pairwise distance residuals
# -------------------------------------------------------------------------

# Heatmap of all pairwise distance residuals
pdf("exploratory/pairwise_dist_residuals_all.pdf")
n_colors = 100  # number of color intervals in scale
color_scale = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(n_colors))
heatmap.2(residuals,
	main="Pairwise cell line distance residuals",
	trace="none",
	col=color_scale,
	RowSideColors=colors[match(ccle_match_info$Site.Primary, unique(ccle_match_info$Site.Primary))],
	ColSideColors=colors[match(ccle_match_info$Site.Primary, unique(ccle_match_info$Site.Primary))],
	cexCol=0.1,
	cexRow=0.1,
	breaks=seq(-10, 10, length.out=n_colors+1)
)
dev.off()

# Heatmaps of tissue-specific pairwise distance residuals
for (tissue in unique(ccle_match_info$Site.Primary)) {
	print(tissue)
	# Heatmap of tissue-specific residuals
	pdf(paste0("exploratory/figCompareCCLE/pairwise_dist_residuals_", tissue, ".pdf"))
	heatmap.2(residuals[ccle_match_info$Site.Primary == tissue, ccle_match_info$Site.Primary == tissue],
		main=paste("Pairwise dist residuals", tissue),
		trace="none",
		col=color_scale,
		cexCol=0.5,
		cexRow=0.5,
		breaks=seq(-10, 10, length.out=n_colors+1)
	)
	dev.off()
}

# Pairwise residual cluster with low values. Identified from the primary heatmap.
# ---------------------------------------------------
cell_lines = c("MCAS", "JVM2", "KMS11", "KURAMOCHI", "OPM2", "L363", "OVMANA", "FU97", "OV90", "HUH7", "H1385", "H1395", "H2126", "KLE", "IM95", "COLO684", "H2106", "DMS79", "AGS", "DU4475", "SHP77", "MINO", "REC1", "SKMEL28", "SKMEL2", "RPMI8226", "KMS20", "NCIH929", "AMO1", "MOLP8", "LP1", "H69", "H211", "H810", "H838", "JHH7", "A204", "H446", "H524")

id = match(cell_lines, rownames(residuals))

n_colors = 100  # number of color intervals in scale
color_scale = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(n_colors))
pdf("exploratory/figCompareCCLE/pairwise_dist_residuals_low.pdf")
heatmap.2(residuals[id, id],
	main="Pairwise cell line distance residuals low cluster",
	trace="none",
	col=color_scale,
	RowSideColors=colors[match(ccle_match_info$Site.Primary[id], unique(ccle_match_info$Site.Primary))],
	ColSideColors=colors[match(ccle_match_info$Site.Primary[id], unique(ccle_match_info$Site.Primary))],
	cexCol=0.4,
	cexRow=0.4,
	breaks=seq(-10, 10, length.out=n_colors+1)
)
dev.off()

# Slice of RPPA data associated with the residual cluster with low values.
pdf("exploratory/figCompareCCLE/RPPAslice_residuals_low_cluster.pdf")
x = tcpa_match_filter_mat
x[is.na(x)] = 0.0
heatmap.2(t(x[id,]),
	main="RPPA low residuals",
	trace="none",
	col=color_scale,
	# RowSideColors=colors[match(ccle_match_info$Site.Primary[id], unique(ccle_match_info$Site.Primary))],
	ColSideColors=colors[match(ccle_match_info$Site.Primary[id], unique(ccle_match_info$Site.Primary))],
	cexCol=0.4,
	cexRow=0.2,
	breaks=seq(-5, 5, length.out=n_colors+1)
)
dev.off()


# Slice of RPPA data from breast tissue.
tissue = "breast"
id = ccle_match_info$Site.Primary == tissue

pdf("exploratory/figCompareCCLE/RPPAslice_breast.pdf")
x = tcpa_match_filter_mat
x[is.na(x)] = 0.0
heatmap.2(t(x[id,]),
	main=tissue,
	trace="none",
	col=color_scale,
	# RowSideColors=colors[match(ccle_match_info$Site.Primary[id], unique(ccle_match_info$Site.Primary))],
	# ColSideColors=colors[match(ccle_match_info$Site.Primary[id], unique(ccle_match_info$Site.Primary))],
	cexCol=0.4,
	cexRow=0.2,
	breaks=seq(-5, 5, length.out=n_colors+1)
)
dev.off()


# Tanglegram tree comparison of cell lines
# -------------------------------------------------------------------------

# tissue = "lung"
# tissue = "stomach"
tissue = "skin"
# tissue = "breast"

id = which(ccle_match_info$Site.Primary == tissue)

# par(mfrow=c(2, 1))
# for a guide on how to use "dendextend" see:
# http://www.r-statistics.com/2014/07/the-dendextend-package-for-visualizing-and-comparing-trees-of-hierarchical-clusterings-slides-from-user2014/
# for another intro see
# https://cran.r-project.org/web/packages/dendextend/vignettes/introduction.html

# nodePar = list(lab.cex=0.6, pch=c(NA, 19), 
# 	cex=0.7, col=c("blue", "red", "green"))

# Examples of how to use dendextend
# dend_tcpa = as.dendrogram(hc_tcpa)
# dend_tcpa = color_branches(dend_tcpa, k=5)
# dend_tcpa = color_labels(dend_tcpa,
# 	col=colors[match(ccle_match_info$Hist.Subtype1[id], unique(ccle_match_info$Hist.Subtype1))])

# plot(dend_tcpa)

# Hierarchical clustering of selection
hc_tcpa = hclust(as.dist(tcpa_dist_mat[id, id]))
hc_ccle = hclust(as.dist(ccle_dist_mat[id, id]))

# Convert to dendogram
dend_tcpa = as.dendrogram(hc_tcpa)
dend_ccle = as.dendrogram(hc_ccle)



# Untagle the two dendrograms. First random search then 2 side rotation.
dend_untangle_rand = untangle_random_search(dend_tcpa, dend_ccle)  # returns list in order: RPPA, CCLE
dend_untangle_rotate = untangle_step_rotate_2side(dend_untangle_rand[[1]], dend_untangle_rand[[2]])

pdf(paste0("exploratory/figCompareCCLE/tanglegram_", tissue, ".pdf"))
tanglegram(dend_untangle_rotate[[1]], dend_untangle_rotate[[2]],
	main=paste(tissue, "cancer cell lines"),
	main_left="RPPA", main_right="mRNA",
	lwd=2.0,  # connecting label line width
	color_lines=rgb(0,0,0,0.5),  # connectin line color
	cex_main=2,  
	margin_inner=6,
	columns_width=c(5,2,5),  # relative sizes of 3 plots: left dendrogram, connectin lines, right dendrogram
	k_branches=7,
	# rank_branches=TRUE
)
dev.off()


# RPPA heatmap of matching data only.
# ---------------------------------------------
x = tcpa_match_filter_mat
x[is.na(x)] = 0.0
n_colors = 100  # number of color intervals in scale
color_scale = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(n_colors))

# x = tcpa_match_mat
heatmap.2(t(x),
	trace="none",
	col=color_scale,
	cexCol=0.2,
	cexRow=0.2,
	breaks=seq(-2, 2, length.out=n_colors+1)
)




# Ordered residual boxplots, per cell line.
pdf("exploratory/figCompareCCLE/residual_ranking.pdf")
par(mfrow=c(2, 1))
id = order(apply(residuals, 2, median), decreasing=T)[1:30]
boxplot(residuals[,id], las=2,
	ylab="RPPA-mRNA residuals",
	col=colors[match(ccle_match_info$Site.Primary[id], unique(ccle_match_info$Site.Primary))]
	)

abline(0,0, lty=3, col="grey")

id = order(apply(residuals, 2, median), decreasing=FALSE)[1:30]
boxplot(residuals[,id], las=2,
	ylab="RPPA-mRNA residuals",
	col=colors[match(ccle_match_info$Site.Primary[id], unique(ccle_match_info$Site.Primary))]
	)
abline(0,0, lty=3, col="grey")
dev.off()


# Centroid distance calculation
# -------------------------------------------------

# Calculates the euclidean centroid distances from the columns in mat to a provided centroid
centroidDist = function(mat, centroid) {
	if (nrow(mat) != length(centroid)) {
		stop("Dimensions of matrix and centroid disagree.")
	}

	distances = sqrt(
		apply(
			(mat - centroid)^2,
			2,
			sum, na.rm=TRUE)
	)
	return(distances)
}

# Combine every other element of two lists
braidLists = function(list1, list2) {
	if (length(list1) != length(list2)) {
		stop("List dimension disagree.")
	}

	out = list()
	for (i in 1:length(list1)) {

		out[[i * 2 - 1]] = list1[[i]]
		out[[i * 2]] = list2[[i]]
	}

	return(out)
}

# mat is a features x samples matrix
# cl is a vector of cluster identities 
centroidMat = function(mat, group) {

	if (ncol(mat) != length(group)) {
		stop("Matrix and group dimension mismatch.")
	}

	# matrix of centroids
	centroids = matrix(NA, ncol=length(group), nrow(mat))

	colnames(centroids) = group

	for (i in 1:length(unique(group))) {
		
		# find sample columns of group
		cols = which(group == unique(group)[i])
	}
}

# ccle_dist_scaled = ccle_dist
# tcpa_dist_scaled = sum(ccle_dist) / sum(tcpa_dist) * tcpa_dist

# Calculate tissue centroids

# Init centroid matrices
tcpa_tissue_centroid = matrix(NA,
	nrow=ncol(tcpa_match_filter_mat),
	ncol=length(unique(ccle_match_info$Site.Primary)))
colnames(tcpa_tissue_centroid) = unique(ccle_match_info$Site.Primary)

ccle_tissue_centroid = matrix(NA,
	nrow=nrow(ccle_match),
	ncol=length(unique(ccle_match_info$Site.Primary)))
colnames(ccle_tissue_centroid) = unique(ccle_match_info$Site.Primary)

for (i in 1:length(unique(ccle_match_info$Site.Primary))) {
	tissue = unique(ccle_match_info$Site.Primary)[i]
	id = which(ccle_match_info$Site.Primary == tissue)

	tcpa_tissue_centroid[,i] = apply(tcpa_match_filter_mat[id,], 2, mean, na.rm=TRUE)
	ccle_tissue_centroid[,i] = apply(ccle_match[,id], 1, mean, na.rm=TRUE)
}


# For sample belonging to a tissue type, calculate the distance to the tissue centroid.
tcpa_dist_centroid = list()
ccle_dist_centroid = list()
for (i in 1:length(unique(ccle_match_info$Site.Primary))) {
	tissue = as.character(unique(ccle_match_info$Site.Primary)[i])
	id = which(ccle_match_info$Site.Primary == tissue)

	tcpa_tissue = t(tcpa_match_filter_mat[id,])
	ccle_tissue = ccle_match[,id]

	tcpa_dist_centroid[[tissue]] = centroidDist(tcpa_tissue, tcpa_tissue_centroid[,i])
	ccle_dist_centroid[[tissue]] = centroidDist(ccle_tissue, ccle_tissue_centroid[,i])
}


tcpa_dist_scale = mean(unlist(tcpa_dist_centroid))
ccle_dist_scale = mean(unlist(ccle_dist_centroid))
tcpa_dist_centroid_scaled = lapply(tcpa_dist_centroid, function(distances) {
	return(distances / tcpa_dist_scale)
})

ccle_dist_centroid_scaled = lapply(ccle_dist_centroid, function(distances) {
	return(distances / ccle_dist_scale)
})


# par(mfrow=c(2, 1))
# boxplot(tcpa_dist_centroid_scaled, las=2)
# boxplot(ccle_dist_centroid_scaled, las=2)

# barplot(
# 	as.vector(rbind(
# 		sapply(tcpa_dist_centroid_scaled, mean),
# 		sapply(ccle_dist_centroid_scaled, mean)
# 	)),
# 	las=2)

# par(mfrow=c(2, 1))
# boxplot(tcpa_dist_centroid, las=2)
# boxplot(ccle_dist_centroid, las=2)

# Combine lists of dist centroids
comb_dist_centroid_scaled = braidLists(tcpa_dist_centroid_scaled, ccle_dist_centroid_scaled)

# barplot(sapply(comb_dist_centroid_scaled, mean))

# pairwise t-test, filter the comparisons shown.
ttest = list()
ftest = list()
for (i in 1:length(tcpa_dist_centroid_scaled)) {
	try({
		ttest[[i]] = t.test(tcpa_dist_centroid_scaled[[i]], ccle_dist_centroid_scaled[[i]])
	})

	try({
		ftest[[i]] = var.test(tcpa_dist_centroid_scaled[[i]], ccle_dist_centroid_scaled[[i]])
	})
}

ttest_pvalues = sapply(ttest, function(test) {
	if (!is.null(test$p.value)) {
		return(test$p.value)
	} else {
		return(1.0)
	}
})

ftest_pvalues = sapply(ftest, function(test) {
	if (!is.null(test$p.value)) {
		return(test$p.value)
	} else {
		return(1.0)
	}
})

# Calculate tissue order based on t-test and F-test status groups: t-test, both, F-test
tissue_order = c(
	# significantly different means
	which(ttest_pvalues < 0.05 & ftest_pvalues >= 0.05),
	# significant means and variances
	which(ttest_pvalues < 0.05 & ftest_pvalues < 0.05),
	# sig variances
	which(ttest_pvalues >= 0.05 & ftest_pvalues < 0.05)
)

n_mean_sig = length(which(ttest_pvalues < 0.05 & ftest_pvalues >= 0.05))
n_both = length(which(ttest_pvalues < 0.05 & ftest_pvalues < 0.05))
n_var_sig = length(which(ttest_pvalues >= 0.05 & ftest_pvalues < 0.05))

separators = c(2*n_mean_sig, 2*n_mean_sig + 2*n_both)  # for drawing vertical lines indicating signifiance groups

# Get list for plotting figure
# comb_dist_centroid_scaled_sig = braidLists(tcpa_dist_centroid_scaled[ttest_pvalues < 0.05], ccle_dist_centroid_scaled[ttest_pvalues < 0.05])
comb_dist_centroid_scaled_sig = braidLists(tcpa_dist_centroid_scaled[tissue_order], ccle_dist_centroid_scaled[tissue_order])

# mRNA-RPPA colors associated with 
col = unlist(braidLists(
	sapply(tcpa_dist_centroid_scaled[tissue_order], function(nums) {
		return(rep(rgb(31, 120, 180, maxColorValue=255), length(nums)))
	}),
	sapply(ccle_dist_centroid_scaled[tissue_order], function(nums) {
		return(rep(rgb(51, 160, 44, maxColorValue=255), length(nums)))
	})
))

pdf("exploratory/figCompareCCLE/rel_dist_tissue_centroid.pdf", height=4, width=7)
pair_offset = -0.3  # how much RPPA and mRNA pairs are closer horizontally
plot(
	# x coordinates of comparative jitter plot
	jitter(rep(
		1:length(comb_dist_centroid_scaled_sig) + rep(c(0, pair_offset), length(comb_dist_centroid_scaled_sig)/2),
		sapply(comb_dist_centroid_scaled_sig, length)
	)),
	# y coordinate (distance to centroid)
	unlist(comb_dist_centroid_scaled_sig),
	main="Significant mean (t-test) and variance (F-test). blue: RPPA, green: mRNA",
	ylab="Relative distance to tissue centroid",
	xlab=paste(names(tcpa_dist_centroid_scaled[tissue_order]), collapse=", "),
	xaxt="n",  # no x-axis labels
	bty="n",  # no bounding box
	col=col
	# cex.main=0.6
	)

abline(1,0, lty=3, col="grey")  # baseline
abline(v=separators[1] + 0.3, col="grey")
abline(v=separators[2] + 0.3, col="grey")

# Mean indication
means = sapply(comb_dist_centroid_scaled_sig, function(nums) {
	return(mean(nums))
})

line_width = 0.3
for (i in 1:length(means)) {
	if (i %% 2 == 0) {
		lines(c(i - line_width + pair_offset, i + line_width + pair_offset), c(means[i], means[i]),
			col="black",
			lwd=1.5)
	} else {
		lines(c(i - line_width, i + line_width), c(means[i], means[i]),
			col="black",
			lwd=1.5)
	}
}
dev.off()


# Between-tissue distances, comparative hierarchical clustering.
# ----------------------------------------------------
tcpa_centroid_dist = dist(t(tcpa_tissue_centroid))
ccle_centroid_dist = dist(t(ccle_tissue_centroid))


hc_tcpa_centroid = hclust(tcpa_centroid_dist)
hc_ccle_centroid = hclust(ccle_centroid_dist)

# Convert to dendogram
dend_tcpa = as.dendrogram(hc_tcpa_centroid)
dend_ccle = as.dendrogram(hc_ccle_centroid)

dend_tcpa %>% set("nodes_col", c(3,4))

# Untagle the two dendrograms. First random search then 2 side rotation.
dend_untangle_rand = untangle_random_search(dend_tcpa, dend_ccle)  # returns list in order: RPPA, CCLE
dend_untangle_rotate = untangle_step_rotate_2side(dend_untangle_rand[[1]], dend_untangle_rand[[2]])

pdf(paste0("exploratory/figCompareCCLE/tanglegram_tissue_centroid.pdf"))
tanglegram(dend_untangle_rotate[[1]], dend_untangle_rotate[[2]],
	main=paste("Tissue centroids"),
	main_left="RPPA", main_right="mRNA",
	lwd=2.0,  # connecting label line width
	# color_lines=rgb(0,0,0,0.5),  # connectin line color
	# color_lines=c("red", "blue"),

	color_lines=colors[order.dendrogram(dend_untangle_rotate[[1]])],
	cex_main=2,  
	margin_inner=6,
	columns_width=c(5,2,5)  # relative sizes of 3 plots: left dendrogram, connectin lines, right dendrogram
	# k_branches=7,
	# rank_branches=TRUE
)
dev.off()


# K-nearest neighbor cross-validation estimate of confusion matrix
# ---------------------------------------------------------------

# tcpa_match_mat_filter
# ccle_match


# ccle_match_info$Site.Primary

x = tcpa_match_filter_mat
x[is.na(x)] = 0.0

tcpa_knn_pred = knn.cv(x, factor(ccle_match_info$Site.Primary), k=3)  # leave-one-out cross validation
tcpa_conf_mat = confusion(tcpa_knn_pred, factor(ccle_match_info$Site.Primary))


ccle_knn_pred = knn.cv(t(ccle_match), factor(ccle_match_info$Site.Primary), k=3)  # leave-one-out cross validation
ccle_conf_mat = confusion(ccle_knn_pred, factor(ccle_match_info$Site.Primary))

pdf("exploratory/figCompareCCLE/cv_tissue_errors_knn.pdf", width=8, height=12)
par(mfrow=c(5, 4))

for (tissue in unique(ccle_match_info$Site.Primary)) {
	# Tissue index
	index = colnames(ccle_conf_mat) == tissue


	n = sum(tcpa_conf_mat[,index])

	if (n != sum(ccle_conf_mat[,index])) {
		stop("Sample number disagrees.")
	}

	tcpa_mistakes = tcpa_conf_mat[!index, index]
	tcpa_mistakes = tcpa_mistakes[tcpa_mistakes > 0]

	ccle_mistakes = ccle_conf_mat[!index, index]  # error counts for tissue type
	ccle_mistakes = ccle_mistakes[ccle_mistakes > 0]

	ccle_mistakes_sorted = sort(ccle_mistakes, decreasing=TRUE)
	tcpa_mistakes_sorted = sort(tcpa_mistakes, decreasing=TRUE)


	plot(c(0, max(8, length(tcpa_mistakes), length(ccle_mistakes)) + 2), range(c(0, 1, -tcpa_mistakes, tcpa_mistakes, -ccle_mistakes, ccle_mistakes)), type = "n",
		main=paste0(tissue, " , n=", n),
		ylab="CV misclassifications",
		xlab="",
		xaxt="n",  # no x-axis labels
		bty="n",  # no bounding box
		cex.main=0.7
	)

	abline(h=0, col="grey")

	barplot(tcpa_mistakes_sorted,
		col=colors[match(names(tcpa_mistakes_sorted), unique(ccle_match_info$Site.Primary))],
		names.arg=rep("", length(tcpa_mistakes)),  # no labels
		# border=rgb(0.9, 0.9, 0.9),
		border=NA,
		add=TRUE,
		axes=FALSE)

	barplot(-ccle_mistakes_sorted,
		col=colors[match(names(ccle_mistakes_sorted), unique(ccle_match_info$Site.Primary))],
		names.arg=rep("", length(ccle_mistakes)),  # no labels
		# border=rgb(0.9, 0.9, 0.9),
		border=NA,
		add=TRUE,
		axes=FALSE)
}
dev.off()



# boxplot(comb_dist_centroid_scaled_sig)
# boxplot(comb_dist_centroid_scaled)

# plot(
# 	# x coordinates of comparative jitter plot
# 	jitter(rep(
# 		1:length(comb_dist_centroid_scaled),
# 		sapply(comb_dist_centroid_scaled, length)
# 	)),
# 	# y coordinate (distance to centroid)
# 	unlist(comb_dist_centroid_scaled)
# 	)

# i = 16
# plot(
# 	c(
# 		jitter(rep(1, length(tcpa_dist_centroid_scaled[[i]]))),
# 		jitter(rep(2, length(ccle_dist_centroid_scaled[[i]])))
# 	),
# 	c(
# 		tcpa_dist_centroid_scaled[[i]],
# 		ccle_dist_centroid_scaled[[i]]
# 	)
# )



# within tissue slice
for (tissue in unique(ccle_match_info$Site.Primary)) {
	print(tissue)

	within_residuals = residuals[ccle_match_info$Site.Primary == tissue, ccle_match_info$Site.Primary == tissue]
	between_residuals = residuals[ccle_match_info$Site.Primary == tissue, ccle_match_info$Site.Primary != tissue]

	print(t.test(as.vector(within_residuals), as.vector(between_residuals)))
}

within_residuals = residuals[ccle_match_info$Site.Primary == tissue, ccle_match_info$Site.Primary == tissue]
between_residuals = residuals[ccle_match_info$Site.Primary == tissue, ccle_match_info$Site.Primary != tissue]

t.test(as.vector(within_residuals), as.vector(between_residuals))

# image(residuals[ind[1,], ind[2,]])


i = which(rownames(residuals) == "AU565")
# plot(tcpa_cor[i,], ccle_dist[i,])

# residuals[i,]


id = ccle_match_info$Site.Primary == tissue

# x = data.frame(tcpa=tcpa_cor[i,], ccle=ccle_dist[i,], res=residuals[i,])
# x = x[ccle_match_info$Site.Primary == tissue,]

# x = data.frame(tcpa=as.vector(tcpa_cor[id,id]), ccle=as.vector(ccle_dist[id,id]), res=as.vector(residuals[id,id]))

x = data.frame(tcpa=as.vector(tcpa_cor), ccle=as.vector(ccle_dist), res=as.vector(residuals))


x = x[!x$tcpa > 0.99999,]

x = x[order(x$res, decreasing=T),]

plot(x$tcpa, x$ccle, main=paste(rownames(residuals)[i]))
points(x$tcpa[1:20], x$ccle[1:20], col="red")
points(x$tcpa[(nrow(x)-20+1):nrow(x)], x$ccle[(nrow(x)-20+1):nrow(x)], col="blue")


# Averages of averages
x = matrix(rnorm(5*30), nrow=5, ncol=30)

apply(x, 1, sd)
apply(x, 1, mean)
