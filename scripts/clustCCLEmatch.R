rm(list=ls())

library(data.table)
library(dendextend)

setwd("/Users/sk/Desktop/tcpa-rppa")

source("lib/annotate.r")  # annotation mapping functions
source("lib/clustVal.r")
source("lib/basePlots.r")
source("lib/colors.r")  # Defines colors


# Load CCLE and MCLP (TCPA) data
# ------------------------------------------------------------
ccle = fread("data/CCLE_Expression_Entrez_2012-09-29.gct")
ccle = as.data.frame(ccle)  # convert from data.table to data.frame
# ccle$Description

# Sample descriptions
ccle_info = read.table("data/CCLE_sample_info_file_2012-10-18.txt", header=TRUE, sep="\t")

# load TCPA data
tcpa = read.table("data/MCLP-RBN-v1.0-whole_set.tsv", header=TRUE, sep="\t")

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


# TCPA

# Format TCPA data matrix (separate data from sample information).
tcpa_match_mat = tcpa_match[,4:ncol(tcpa_match)]
tcpa_match_mat = as.matrix(tcpa_match_mat)
rownames(tcpa_match_mat) = tcpa_match[,1]  #

# filter by number of observations
min_cell_lines = 100  # the minimum number of cell lines measured for an included protein

# count the number of non missing data entries for a particular protein.
protein_counts = apply(tcpa_match_mat, 2, function(col) {
	return(sum(!is.na(col)))
})

tcpa_match_filter_mat  = tcpa_match_mat[,protein_counts > min_cell_lines]


# Calculate distance matrices of matching CCLE and TCPA data
# --------------------------------------------------------

# tcpa_dist = dist(tcpa_match_filter_mat)
# ccle_dist = dist(t(ccle_match))


# Filter matching CCLE gene expression data based on probe variance
# n_probes = 500  # the number of probes to include
n_probes = ncol(tcpa_match_filter_mat)

ccle_gene_var = apply(ccle_match, 1, var)
high_var_probes = sort(ccle_gene_var, index.return=TRUE, decreasing=TRUE)

ccle_match_filter = ccle_match[high_var_probes$ix[1:n_probes],]

# Filter CCLE genes (array probes) based on protein symbols from the TCPA
include_probes = which(ccle$Description %in% colnames(tcpa_match_filter_mat))
ccle_match_same_genes = ccle_match[include_probes,]


# Gap statistic calculation
# ----------------------------------------
max_clusters = 20  # maximum number of clusters considered
n_null = 10  # the number of reference (null) datasets. Notation: B in (Tibshirani, 2001).


gap_tcpa = gapStatHc(mat=tcpa_match_filter_mat, max_clusters=max_clusters, n_null=n_null)
gap_ccle = gapStatHc(mat=t(ccle_match), max_clusters=max_clusters, n_null=n_null)

gap_ccle_filter = gapStatHc(mat=t(ccle_match_filter), max_clusters=max_clusters, n_null=n_null)

gap_ccle_same_genes = gapStatHc(mat=t(ccle_match_same_genes), max_clusters=max_clusters, n_null=n_null)


green = rgb(51, 160, 44, maxColorValue=255)
blue = rgb(31, 120 , 180, maxColorValue=255)
lightblue = rgb(166, 206 , 227, maxColorValue=255)
lightgreen = rgb(178, 223 , 138, maxColorValue=255)

pdf("exploratory/figClustCCLEmatch/gap_statistic_CCLE_match.pdf", width=7, height=4)
par(mfcol=c(2, 4))
gapPlot(gap_tcpa, col=blue, main=paste(ncol(tcpa_match_filter_mat), "RPPA"))
gapPlot(gap_ccle, col=green, main=paste(nrow(ccle_match), "mRNA"))
gapPlot(gap_ccle_filter, col=lightgreen, main=paste(nrow(ccle_match_filter), "high var mRNA"))
gapPlot(gap_ccle_same_genes, col=lightgreen, main=paste(nrow(ccle_match_same_genes), "matching mRNA"))
dev.off()



dmat = dist(t(ccle_match))
hc = hclust(dmat, method="average")

# clusters = cutree(hc, 3)
clusters = cutree(hc, 14)
sort(clusters)

hc$labels = ccle_match_info$Cell.line.primary.name

dend = as.dendrogram(hc)
dend = hang.dendrogram(dend)

labels_cex(dend) = 0.28

# dend = color_branches(dend, 5, colors)
dend = color_labels(dend,
	col=colors[match(ccle_match_info$Site.Primary, unique(ccle_match_info$Site.Primary))][hc$order]
)



# RPPA
# -----------------------------------
dmat2 = dist(tcpa_match_filter_mat)
hc2 = hclust(dmat2, method="average")
# clusters2 = cutree(hc2, 6)
clusters2 = cutree(hc2, 13)
names(clusters2) = ccle_match_info$CCLE.name
sort(clusters2)

dend2 = as.dendrogram(hc2)
dend2 = hang.dendrogram(dend2)

labels_cex(dend2) = 0.28

dend2 = color_labels(dend2,
	col=colors[match(ccle_match_info$Site.Primary, unique(ccle_match_info$Site.Primary))][hc2$order]
)

pdf("exploratory/figClustCCLEmatch/cell_line_dendrogram.pdf", height=16)
par(mfrow=c(1, 2))
plot(dend, main="mRNA", horiz=TRUE)

plot(dend2, main="RPPA", horiz=TRUE)
dev.off()


# plot(hc, cex=0.35, col=c("red", "blue"))

# names(clusters2)
# colnames(ccle)

ccle_match_info
data.frame(n1=names(clusters2), n2=ccle_match_info$CCLE.name)
data.frame(n1=names(clusters2), n2=ccle_match_info$CCLE.name, )

data.frame(name=names(clusters), cl1=clusters, cl2=clusters2)

# Cophonetic analysis
# --------------------------------------------------------------

tcpa_dist = dist(tcpa_match_filter_mat)
hc_tcpa = hclust(tcpa_dist, method="average")
tcpa_cophen_dist = cophenetic(hc_tcpa)

ccle_dist = dist(t(ccle_match))
hc_ccle = hclust(ccle_dist, method="average")
ccle_cophen_dist = cophenetic(hc_ccle)

# high variance filter
ccle_filter_dist = dist(t(ccle_match_filter))
hc_ccle_filter = hclust(ccle_filter_dist, method="average")
ccle_filter_cophen_dist = cophenetic(hc_ccle_filter)

# random sample of mRNA features
ccle_match_sample = ccle_match[sample(1:nrow(ccle_match), ncol(tcpa_match_filter_mat)),]
ccle_sample_dist = dist(t(ccle_match_sample))
hc_ccle_sample = hclust(ccle_sample_dist, method="average")
ccle_sample_cophen_dist = cophenetic(hc_ccle_sample)

# matched mRNA probes
ccle_same_genes_dist = dist(t(ccle_match_same_genes))
hc_ccle_same_genes = hclust(ccle_same_genes_dist, method="average")
ccle_same_genes_cophen_dist = cophenetic(hc_ccle_same_genes)

cor_test = list()
cor_test[["RPPA"]] = cor.test(tcpa_dist, tcpa_cophen_dist)
cor_test[["mRNA"]] = cor.test(ccle_dist, ccle_cophen_dist)
cor_test[["mRNA var. filter"]] = cor.test(ccle_filter_dist, ccle_filter_cophen_dist)
cor_test[["mRNA random"]] = cor.test(ccle_sample_dist, ccle_sample_cophen_dist)
cor_test[["mRNA same genes"]] = cor.test(ccle_same_genes_dist, ccle_same_genes_cophen_dist)


green = rgb(77, 175, 74, 200, maxColorValue=255)
lightgreen = rgb(77, 175, 74, 150, maxColorValue=255)
blue = rgb(55, 126 , 184, 200, maxColorValue=255)

# lightblue = rgb(166, 206 , 227, maxColorValue=255)
# lightgreen = rgb(178, 223 , 138, maxColorValue=255)

pdf("exploratory/figClustCCLEmatch/cophenetic_clust.pdf", width=1.8, height=4)
err_width = 0.02
bar_centers = barplot(
	sapply(cor_test, function(test) test$estimate),
	names.arg=names(cor_test),
	col=c(
		blue,
		green,
		lightgreen,
		lightgreen,
		lightgreen,
		lightgreen
	),
	cex.names=0.5,
	las=2,
	beside=TRUE,
	ylab="Cophenetic correlation",
	ylim=c(0.5, 1.0),
	xpd=FALSE,  # don't show outside
	border=NA
	)

plotErrBars(
	bar_centers,
	lapply(cor_test, function(test) test$conf.int),
	err_width)

dev.off()

# CCLE
# hc_ccle = hclust(ccle_dist, method="average")  # CCLE

# Plots

par(mfrow=c(2, 1))
# plot(log(within_tcpa, 2), type="l")
plot(within_tcpa, type="l")

plot(gap, type="l", col=green)
# points(gap, pch=16, cex=0.7, col=green)
segments(c(1:max_clusters), gap - sd_log_within_ref, c(1:max_clusters), gap + sd_log_within_ref, col=green)
# top and bottom bar
segments(c(1:max_clusters) - seg_width, gap + sd_log_within_ref, c(1:max_clusters) + seg_width, gap + sd_log_within_ref, col=green)
segments(c(1:max_clusters) - seg_width, gap - sd_log_within_ref, c(1:max_clusters) + seg_width, gap - sd_log_within_ref, col=green)



# lines(gap + sd_log_within_ref, type="l", col="grey")
# lines(gap - sd_log_within_ref, type="l", col="grey")
# polygon(c(1:max_clusters, max_clusters:1), c(gap + sd_log_within_ref, rev(gap - sd_log_within_ref)),
# 	col="grey",
# 	border=NA)

# lines(gap)


# apply(log(within_ref, 2), 2, mean)
# apply(log(within_ref, 2), 2, sd)




plot(within_ref[[1]], type="l")
lines(within_ref[[2]], type="l")
lines(within_ref[[3]], type="l")
lines(within_ref[[4]], type="l")
lines(within_ref[[5]], type="l")







hc = hclust(dist(mat), method="average")
plot(hc)

plot(mat[1,], mat[2,])

i = 16
j = 5
par(mfrow=c(2, 1))
plot(tcpa_match_filter_mat[,i], tcpa_match_filter_mat[,j])
plot(mat[,i], mat[,j])



cuts = list()
for (k in 1:max_clusters) {
	cuts[[k]] = cutree(hc, k=k)
}

ws = sapply(cuts, function(groups) {
	return(withinClustSquares(dist(mat), groups))
})
plot(log(ws, 2), type="l")





