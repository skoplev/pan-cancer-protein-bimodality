# Bimocality analysis of TCPA data.

# Layered deconvolution of TCPA/MCLP (RPPA) data from MD Anderson.

rm(list=ls())

library(class)
library(data.table)
library(RCurl)
library(XML)
library(RColorBrewer)

setwd("~/Google Drive/projects/tcpa-rppa")
data_dir = "~/DataProjects/tcpa-rppa"

source("lib/annotate.r")  # annotation mapping functions
source("lib/clustVal.r")
source("lib/bimodal.r")


# Load MCLP protein expression data
tcpa = read.table(file.path(data_dir, "data/MCLP-RBN-v1.0-whole_set.tsv"), header=TRUE, sep="\t")

# Format data matrix, cell lines x protein measurements
x = tcpa[,4:ncol(tcpa)]
x = as.matrix(x)
rownames(x) = tcpa[,1]  #


# Complete fraction
sum(!is.na(x)) / length(x)

# cell_line = "K562"
# cell_line = "HEPG2"
# cell_line = "A549"
# cell_line = "MCF7"
# cell_line = "PC3"
# cell_line = "OCILY7"
# cell_line = "HCT116"
# protein = "ECADHERIN"
# protein = "VIMENTIN"

# x[rownames(x) == cell_line, colnames(x) == protein]
# hist(x[, colnames(x) == protein])
# x[, colnames(x) == protein]

# Load matching CCLE data
# ------------------------------------------------------------
ccle_table = fread(file.path(data_dir, "data/CCLE_Expression_Entrez_2012-09-29.gct"))
ccle_table = as.data.frame(ccle_table)  # convert from data.table to data.frame

# Sample descriptions
ccle_info = read.table(file.path(data_dir, "data/CCLE_sample_info_file_2012-10-18.txt"), header=TRUE, sep="\t")

# Format CCLE data
ccle_mat = ccle_table[,3:ncol(ccle_table)]
rownames(ccle_mat) = make.names(ccle_table$Description, unique=TRUE)

# Match metadata to ccle_table
ccle_info_id = match(colnames(ccle_table), ccle_info$CCLE.name)
ccle_info_id = ccle_info_id[!is.na(ccle_info_id)]
ccle_info = ccle_info[ccle_info_id,]

# get cell line names of all CCLE data.
# ccle_cell_lines = ccle_info_align$Cell.line.primary.name[ccle_info_id]

# match TCPA cell line name to CCLE cell line names.
# find corresponding CCLE for each TCPA cell line.
# ccle_match_id = match(stemName(tcpa$Sample_Name), stemName(ccle_info_align$Cell.line.primary.name))

# filter based on protein counts
# ----------------------------------------------------------
# count the number of non missing data entries for a particular protein.
protein_counts = apply(x, 2, function(col) {
	return(sum(!is.na(col)))
})

min_cell_lines = 40  # the minimum number of cell lines measured for an included protein
x = x[,protein_counts > min_cell_lines]  # remove protein without sufficient measurements

# Complete fraction after filtering
sum(!is.na(x)) / length(x)


# Cell line names, matching the data matrix
cell_line_names = tcpa$Sample_Name

# load CCLE annotations
annot = list()
annot$ccle = read.table(file.path(data_dir, "data/CCLE_sample_info_file_2012-10-18.txt"), header=TRUE, sep="\t")

# CCLE annotations
matches = match(stemName(cell_line_names), stemName(annot$ccle$Cell.line.primary.name))
ccle = annot$ccle[matches,]  # data frame of ccle annotations matchine the MCLP data

# Annotate proteins
# ------------------------------

# Predict missing annotations from CCLE based on the MCLP data.
# Nearest neighbor prediction.
# -----------------------------------------------------------------
missing = is.na(ccle$Cell.line.primary.name) 

ccle$inferred = missing  # keep track of which tissue types are inferred directly annotated.

# copy data matrix with missing values replaced with zero
x_zeroth = x
x_zeroth[is.na(x_zeroth)] = 0.0

# K-nearest neighbor prediction
knn_prediction = class::knn(
	train=x_zeroth[!missing,],
	test=x_zeroth[missing,],
	cl=factor(ccle$Site.Primary[!missing]),
	k=3,
	prob=TRUE
)

# Fill in infered tissue 
ccle$Site.Primary[missing] = knn_prediction

proteins = makeProteinTable(colnames(x))  # downloads additional annotations for proteins in RPPA dataset

# Univariate analysis of distributions
# ------------------------------------------
require(moments)
require(RColorBrewer)
require(mixtools)
require(gplots)
require(igraph)
require(beeswarm)
require(qgraph)


# Calculate protein-specific moments, means, variances, ...
# And other basic calculations, the outcomes of which are shared elsewhere in the script.
# -----------------------------------------------------------
prot_means = apply(x, 2, mean, na.rm=TRUE)
prot_means_sort = sort(prot_means, decreasing=TRUE, index.return=TRUE)

prot_var = apply(x, 2, var, na.rm=TRUE)
prot_var_sort = sort(prot_var, decreasing=TRUE, index.return=TRUE)

prot_sd = apply(x, 2, sd, na.rm=TRUE)

prot_skew = apply(x, 2, skewness, na.rm=TRUE)
prot_skew_sort = sort(prot_skew, decreasing=TRUE, index.return=TRUE)

prot_kurt = apply(x, 2, kurtosis, na.rm=TRUE)
prot_kurt_sort = sort(prot_kurt, decreasing=TRUE, index.return=TRUE)


# Fit Gaussian models univariately to the data. One model fitted for each protein observation.
# Includes calculation of the log-likelihood of the two-component Gaussian mixture model and a unimodal Gaussian model.
# ---------------------------------------------------------------

# Fit Gaussian mixture models
em = fitGMM(x)

em_ccle = fitGMM(t(ccle_mat))

# Evaluate probabilities of the low bin assignment
prob_low = probLowModelMat(x, em)
prob_low_ccle = probLowModelMat(t(ccle_mat), em_ccle)  # WARNING: takes 5 min to execute. Could be optimized.


# Entropy and purity calculation based on classification from the mixture model
# ---------------------------------------------------------

# alpha = 0.5  # significance level for the classification based on a univarite Gaussian mixture model.
alpha = 0.5

class_stat = calcTissueEntropy(prob_low, ccle$Site.Primary, alpha)
class_stat_ccle = calcTissueEntropy(prob_low_ccle, ccle_info$Site.Primary, alpha)

# Likelihood comparisons and tests;
# ------------------------------------------------------------------------------
# Likelihood ratio

# norm_loglik = normLoglik(x)

# diff_loglik = -(em$loglik - norm_loglik)  # higher is in favor of EM fit

# diff_loglik_sort = sort(diff_loglik, decreasing=TRUE, index.return=TRUE)
# diff_loglik_order = match(names(diff_loglik_sort$x), colnames(x))

# Bayesian information criteria (BIC)
# ------------------------------------------------------

bic_cut = 2  # BIC cutoff

delta_bic = deltaBIC(x, em)  # RPPA BIC

delta_bic_ccle = deltaBIC(t(ccle_mat), em_ccle)

bimodal_prot = names(delta_bic)[delta_bic > bic_cut]
bimodal_prot = bimodal_prot[!is.na(bimodal_prot)]

write.table(bimodal_prot, "exploratory/bimodal/bimodal_prot.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)


write.csv(delta_bic, "exploratory/bimodal/delta_bic_RPPA.csv", quote=FALSE, col.names=FALSE)

# Compare BIC between mRNA, protein, and PTM.
# ------------------------------------

# Match BIC based on gene names

# Get RPPA ids for protein-based measurements that has a uniprot mapping
# and is not associated with a PTM.
rppa_ids = which(!is.na(proteins$uniprot_name) & is.na(proteins$mclp_ptm))


# Get matching CCLE transcript data ids
ccle_ids_match = match(proteins$uniprot_name[rppa_ids], rownames(ccle_mat))

# Gene matching statistics
list(
	RPPA=sum(is.na(proteins$mclp_ptm)),
	CCLE=nrow(ccle_mat),
	overlap=sum(!is.na(ccle_ids_match)),
	RPPA_only=sum(is.na(proteins$mclp_ptm)) - sum(!is.na(ccle_ids_match)),
	CCLE_only=nrow(ccle_mat) - sum(!is.na(ccle_ids_match))
)


# data.frame(n1=proteins$uniprot_name[rppa_ids], n2=rownames(ccle_mat)[ccle_ids_match])

# Match cell lines names from RPPA to CCLE
ccle_match_id = match(
	stemName(rownames(prob_low)),
	stemName(ccle_info$Cell.line.primary.name)
)
# Cell line matching stats
list(
	RPPA=nrow(x),
	CCLE=nrow(ccle_info),
	overlap=sum(!is.na(ccle_match_id)),
	RPPA_only=nrow(x) - sum(!is.na(ccle_match_id)),
	CCLE_only=nrow(ccle_info) - sum(!is.na(ccle_match_id))
)

# Align posterior probabilities matching both cell lines and genes
prob_low_align = prob_low[
	which(!is.na(ccle_match_id)),  # the RPPA cell lines found
	rppa_ids  # the matching protein ids
]

prob_low_ccle_align = prob_low_ccle[
	ccle_match_id[!is.na(ccle_match_id)],  # existing ids, removes NA
	ccle_ids_match
]

# Phosphorylation statistics: number of proteins that are also measured by their phosphorylation state
n_protein_ptm_overlap = sum(
	!is.na(
		match(
			proteins$mclp_name[!is.na(proteins$mclp_ptm)],  # protein names of PTM measurements
			proteins$mclp_id  # pure protein id, assumes base id specifies protein measurement (no PTM)
		)
	)
)

# Phosphorylation only measurements, no base protein
phos_only = sum(
	is.na(
		match(
			proteins$mclp_name[!is.na(proteins$mclp_ptm)],
			proteins$mclp_id
		)
	)
)

list(
	proteins=sum(is.na(proteins$mclp_ptm)),  # no PTM antibody measurements
	phos=sum(!is.na(proteins$mclp_ptm)),  # assumes that all PTM "_" annotations are phoshosites
	overlap=n_protein_ptm_overlap,
	proteins_only=sum(is.na(proteins$mclp_ptm)) - n_protein_ptm_overlap,
	phos_only=phos_only
)

# Check that names of cell lines align
# data.frame(n1=rownames(prob_low_align), n2=rownames(prob_low_ccle_align))


# Calculate bin coupling for protein and mRNA measurements
prob_low_cor = cor(prob_low_align, prob_low_ccle_align,
	use="pairwise.complete.obs",
	method="spearman")



# plot(density(diag(prob_low_cor), na.rm=TRUE))
cross_coupling = diag(prob_low_cor)  # mRNA protein cross coupling

# mRNA-protein, BIC comparison scatter plot
# ----------------------------------------------

delta_bic_mrna = delta_bic_ccle[ccle_ids_match]
delta_bic_prot = delta_bic[rppa_ids]

data.frame(n1=names(delta_bic_mrna), n2=names(delta_bic_prot))

svg("exploratory/figCouplingCCLE/delta_bic_mRNA_prot.svg", width=6, height=6)
plotDeltaBIC(delta_bic_mrna, delta_bic_prot, cross_coupling,
	xlab=expression(paste("mRNA bimodality, sign ", log["10"], "(", Delta, "BIC + 1", ")")),
	ylab=expression(paste("Protein bimodality, sign ", log["10"], "(", Delta, "BIC + 1", ")"))
)
dev.off()

log10(delta_bic_mrna + 1)

data.frame(delta_bic_mrna, delta_bic_prot)

idx = signPseudoLog(delta_bic_mrna) > 1 & signPseudoLog(delta_bic_prot) > 1
mean(idx, na.rm=TRUE)

# Transcriptional bimodality
idx = signPseudoLog(delta_bic_mrna) > 1 & signPseudoLog(delta_bic_prot) > 1 & cross_coupling > 0.5
mean(idx, na.rm=TRUE)

# idx = signPseudoLog(delta_bic_mrna) <= 1 & signPseudoLog(delta_bic_prot) > 1 & cross_coupling <= 0.5
# Translational/post-translational bimodality
idx = signPseudoLog(delta_bic_mrna) <= 1 & signPseudoLog(delta_bic_prot) > 1
mean(idx, na.rm=TRUE)

idx = signPseudoLog(delta_bic_mrna) > 1 & signPseudoLog(delta_bic_prot) <= 1
mean(idx, na.rm=TRUE)

idx = signPseudoLog(delta_bic_mrna) <= 1 & signPseudoLog(delta_bic_prot) <= 1
mean(idx, na.rm=TRUE)





# Bin assignment counts
# ---------------------------------------------

alpha = 0.1  # significance level of posterior bin assignment

low_bin = prob_low_align > 1 - alpha

high_bin = prob_low_align < alpha

low_bin_ccle = prob_low_ccle_align > 1 - alpha

high_bin_ccle = prob_low_ccle_align < alpha

coupling_counts = list()
for (i in 1:ncol(low_bin)) {
	counts = list()
	counts[["--"]] = sum(low_bin_ccle[,i] & low_bin[,i], na.rm=TRUE)
	counts[["++"]] = sum(high_bin_ccle[,i] & high_bin[,i], na.rm=TRUE)
	counts[["-+"]] = sum(low_bin_ccle[,i] & high_bin[,i], na.rm=TRUE)
	counts[["+-"]] = sum(high_bin_ccle[,i] & low_bin[,i], na.rm=TRUE)

	coupling_counts[[colnames(low_bin)[i]]] = unlist(counts)
}

# Write coupling counts to file
coupling_counts_mat = Reduce(rbind, coupling_counts)
rownames(coupling_counts_mat) = names(coupling_counts)
colnames(coupling_counts_mat) = c("--", "++", "-+", "+-")

write.csv(coupling_counts_mat, "data/couplingCCLE/coupling_counts.csv", quote=FALSE)


# E-Cadherin translationally regulated cell lines.
# Find and write CCLE info to table.

# Feature id in RPPA and CCLE aligned data
ecad_id = which(colnames(low_bin) == "ECADHERIN")

# Get CCLE id names for 
transl_regulated_cell_lines_id = which(high_bin_ccle[,ecad_id] & low_bin[,ecad_id])
transl_regulated_cell_lines = names(transl_regulated_cell_lines_id)

# Map to CCLE info table
ccle_id = match(transl_regulated_cell_lines, ccle_info[["CCLE.name"]])

ccle_info_ecad_transl = ccle_info[ccle_id,]

# Reorder based on tissue of origin
ccle_info_ecad_transl = ccle_info_ecad_transl[
	order(ccle_info_ecad_transl$Site.Primary),
]

write.csv(ccle_info_ecad_transl, "exploratory/figCouplingCCLE/ECADHERIN_translational_regulated.csv", quote=FALSE)


ccle_ecad = ccle_mat[
	which(rownames(ccle_mat) == "CDH1"),
	ccle_match_id[!is.na(ccle_match_id)]
]
ccle_ecad = unlist(ccle_ecad)

rppa_ecad = x[
	which(!is.na(ccle_match_id)),
	which(colnames(x) == "ECADHERIN")
]
rppa_ecad = unlist(rppa_ecad)

coupling_color = rep("black", length(ccle_ecad))
coupling_color[transl_regulated_cell_lines_id] = "red"

svg("exploratory/figCouplingCCLE/ECADHERIN_mRNA_prot.svg", width=3.5, height=4)
plot(ccle_ecad, rppa_ecad, col=coupling_color,
	xlab="mRNA", ylab="Protein", main="E-Cadherin")
dev.off()

cor.test(ccle_ecad, rppa_ecad)

prob_low_cor[rownames(prob_low_cor) == "ECADHERIN", which(colnames(prob_low_cor) == "CDH1")]


coupling_counts[["ECADHERIN"]]

pdf("exploratory/figCouplingCCLE/EMT_bin_coupling_mrna_protein.pdf", width=1.8, height=10)
orange = rgb(253, 174, 97, maxColorValue=255)
green = rgb(171, 221, 164, maxColorValue=255)
prot_sel = c("ECADHERIN", "CLAUDIN7", "RAB25", "PAI1", "AXL", "HEREGULIN")
# prot_sel = c("ECADHERIN", "CLAUDIN7", "RAB25")
par(mfrow=c(length(prot_sel), 1))
for (i in 1:length(prot_sel)) {
	barplot(coupling_counts[[prot_sel[i]]],
		main=prot_sel[i],
		ylab="Cell lines",
		col=c(orange, orange, green, green),
		border=rgb(0.1, 0.1, 0.1)
	)
}
dev.off()

# Protein-PTM bimodality
# ------------------------------------------------

# Construct matching ids
rppa_phos_ids = which(!is.na(proteins$uniprot_name) & !is.na(proteins$mclp_ptm))

# Match to redundant base protein ids
rppa_phos_prot_ids = match(proteins$mclp_name[rppa_phos_ids], proteins$mclp_id)


delta_bic_phos = delta_bic[rppa_phos_ids]
delta_bic_base_prot = delta_bic[rppa_phos_prot_ids]

# Calculate correlations (coupling) between posterior assignment.
prob_low_phos_cor = cor(
	prob_low[,rppa_phos_ids],
	prob_low[,rppa_phos_prot_ids],
	use="pairwise.complete.obs",
	method="spearman")

svg("exploratory/figCouplingCCLE/delta_bic_prot_phos.svg", width=6, height=6)
plotDeltaBIC(delta_bic_base_prot, delta_bic_phos, diag(prob_low_phos_cor),
	xlab=expression(paste("Protein bimodality, sign ", log["10"], "(", Delta, "BIC + 1", ")")),
	ylab=expression(paste("Phosphosite bimodality, sign ", log["10"], "(", Delta, "BIC + 1", ")"))
)
dev.off()

prob_low_phos_cor[rownames(prob_low_phos_cor) == "HER2_pY1248", colnames(prob_low_phos_cor) == "HER2"]

idx = signPseudoLog(delta_bic_base_prot) > 1 & signPseudoLog(delta_bic_phos) > 1
mean(idx, na.rm=TRUE)

idx = signPseudoLog(delta_bic_base_prot) > 1 & signPseudoLog(delta_bic_phos) > 1 & diag(prob_low_phos_cor) > 0.5
mean(idx, na.rm=TRUE)

idx = signPseudoLog(delta_bic_base_prot) <= 1 & signPseudoLog(delta_bic_phos) > 1
mean(idx, na.rm=TRUE)

idx = signPseudoLog(delta_bic_base_prot) > 1 & signPseudoLog(delta_bic_phos) <= 1
mean(idx, na.rm=TRUE)







data.frame(n1=names(delta_bic_phos), n2=names(delta_bic_base_prot))


# protein-mRNA cross coupling
# ------------------------------------------------------------------------

# beta = 0.2   # fraction defining normal and rare
beta = 0.25   # fraction defining normal and rare
entropy_quant = 1/3  # minimum entropy quantile cutoff


plot(density(em$lambda, na.rm=TRUE))
hist(em$lambda, breaks=50)

prob_low_matched_lines = prob_low[which(!is.na(ccle_match_id)),]
prob_low_ccle_matched_lines = prob_low_ccle[ccle_match_id[!is.na(ccle_match_id)],]

# The complete cross correlation takes a long time to execute.
prob_low_cor_all = cor(
	prob_low_matched_lines,
	prob_low_ccle_matched_lines,
	use="pairwise.complete.obs",
	method="spearman"
)



protein_names = c("ECADHERIN", "VEGFR2", "PAR")  # Selected proteins
protein_names = c("ECADHERIN", "PAI1", "AXL")  # Selected proteins
protein_names = c("ECADHERIN", "PAI1", "NCADHERIN", "PCADHERIN")  # Selected proteins
protein_names = c("ECADHERIN", "PAI1")  # Selected proteins
# protein_names = c("ECADHERIN", "PCADHERIN")  # Selected proteins


# Index of protein
indices = match(protein_names, rownames(prob_low_cor_all))

pdf("exploratory/figCouplingCCLE/protein_mRNA_cross_coupling.pdf", height=8, width=2.5)

# Colors
colors = brewer.pal(9, "Set1")
colors = colors[3:length(colors)]


pdf("exploratory/figCouplingCCLE/protein_mRNA_cross_coupling_density.pdf", height=2.5, width=3.5)
# par(mfrow=c(4, 1))
plot(density(prob_low_cor_all, na.rm=TRUE), col="black", lwd=1.5,
	main="",
	xlab="Protein-mRNA coupling")
for (k in 1:length(indices)) {
	lines(density(prob_low_cor_all[indices[k],], na.rm=TRUE), col=colors[k], lwd=1.5)
}
legend("topright",
	legend=c("All", protein_names),
	col=c("black", colors),
	pch="-", lwd=2, cex=0.5)
dev.off()


# Protein-mRNA trans coupling barplot
# -------------------------------------------------
prot_name = "ECADHERIN"
prot_name = "NCADHERIN"
prot_name = "PAI1"
prot_name = "VEGFR2"
prot_name = "FOXM1"
prot_name = "DUSP4"
prot_name = "CAVEOLIN1"

# Barplots of top n high and low cross correlations
n = 25

pdf(
	paste0("exploratory/figCouplingCCLE/prot-mRNA-barplot/", prot_name, "_protein_mRNA_cross_coupling.pdf"),
		width=2.5)
par(mar=c(5.1, 6, 4.1, 2.1))
# red = rgb(253, 219, 199, maxColorValue=255)
# blue = rgb(209, 229, 240, maxColorValue=255)
purple = rgb(231, 212, 232, maxColorValue=255)
green = rgb(217, 240, 211, maxColorValue=255)
# orange = rgb(254, 224, 139, maxColorValue=255)
# green2 = rgb(230, 245, 152, maxColorValue=255)

cross_cor_sort = sort(
	prob_low_cor_all[
		rownames(prob_low_cor_all) == prot_name,
	]
)

m = length(cross_cor_sort)

barplot(
	abs(rev(c(cross_cor_sort[m:(m - n + 1)], cross_cor_sort[n:1]))),
	main=prot_name,
	cex.names=0.6,
	xlab="Abs. protein-mRNA coupling",
	border=rgb(0.6, 0.6, 0.6),
	col=c(rep(purple, n), rep(green, n)),
	las=2,
	horiz=TRUE
)

abline(v=0.2, col="white", lwd=1)
abline(v=0.4, col="white", lwd=1)
abline(v=0.6, col="white", lwd=1)

dev.off()


# Lists of E-Cadherin associated transcripts
# -----------------------------------------------------
prot_name = "ECADHERIN"
alpha = 0.5  # coupling cutoff for signatures
i = which(rownames(prob_low_cor_all) == prot_name)

pos_genes = colnames(prob_low_cor_all)[prob_low_cor_all[i,] > alpha]
pos_genes = pos_genes[!is.na(pos_genes)]  # Remove no-name genes, from NA posterior probabilities

write.table(pos_genes, paste0("exploratory/figCouplingCCLE/signatures/", prot_name, "_pos.txt"),
	quote=FALSE, row.names=FALSE, col.names=FALSE)

neg_genes = colnames(prob_low_cor_all)[prob_low_cor_all[i,] < -alpha]
neg_genes = neg_genes[!is.na(neg_genes)]

write.table(neg_genes, paste0("exploratory/figCouplingCCLE/signatures/", prot_name, "_inverse.txt"),
	quote=FALSE, row.names=FALSE, col.names=FALSE)


idx = which(abs(prob_low_cor_all[i,]) > alpha)

bimodal_cor = prob_low_cor_all[i, idx]

bimodal_cor_tab = data.frame(gene_symbol=names(bimodal_cor), bimodal_cor=bimodal_cor)

bimodal_cor_tab = bimodal_cor_tab[order(bimodal_cor_tab$bimodal_cor, decreasing=TRUE), ]

write.table(bimodal_cor_tab,
	paste0("exploratory/figCouplingCCLE/signatures/", prot_name, "_signature.csv"),
	sep=",",
	row.names=FALSE,
	quote=FALSE)



# Enrichment analysis of positively coupled transcripts using Enrichr.
# -------------------------------------------------------

n = 20  # top enrichment to consider

# Coloring scheme
purple = rgb(222, 203, 228, maxColorValue=255)
orange = rgb(254, 217, 166, maxColorValue=255)
yellow = rgb(255, 255, 204, maxColorValue=255)
pink = rgb(253, 218, 236, maxColorValue=255)
brown = rgb(229, 216, 189, maxColorValue=255)


# GO Biological process
pos_go_process = fread("exploratory/figCouplingCCLE/emt_signature/up_enrichment/GO_Biological_Process_2015_table.txt")
pos_go_process = as.data.frame(pos_go_process)

# Order by adjusted P-values
pos_go_process = pos_go_process[order(pos_go_process[["Adjusted P-value"]]),]

# GO cellular component
pos_go_comp = fread("exploratory/figCouplingCCLE/emt_signature/up_enrichment/GO_Cellular_Component_2015_table.txt")
pos_go_comp = as.data.frame(pos_go_comp)

# Order by adjusted P-values
pos_go_comp = pos_go_comp[order(pos_go_comp[["Adjusted P-value"]]),]


# GO cellular component
pos_chea = fread("exploratory/figCouplingCCLE/emt_signature/up_enrichment/ChEA_2015_table.txt")
pos_chea = as.data.frame(pos_chea)

# Order by adjusted P-values
pos_chea = pos_chea[order(pos_chea[["Adjusted P-value"]]),]

# Fraction of pos genes explained by association to top-n terms
go_process_genes = unique(
	str_split(
		paste(pos_go_process$Genes[1:n], collapse=";")  # Combined top-n gene associations
		, ";"
	)[[1]]
)

go_process_genes = getUniqueGenes(pos_go_process, n)
go_process_coverage = pos_genes %in% go_process_genes

go_comp = getUniqueGenes(pos_go_comp, n)
go_comp_coverage = pos_genes %in% go_comp

go_chea = getUniqueGenes(pos_chea, n)
go_chea_coverage = pos_genes %in% go_chea


pdf("exploratory/figCouplingCCLE/emt_signature/up_enrichment/bar_plot_pies.pdf", height=1.8)

par(mfrow=c(1, 3))
pie(
	c(sum(go_process_coverage), sum(!go_process_coverage)), 
	labels=c(sum(go_process_coverage), sum(!go_process_coverage)), 
	main="GO biological process",
	col=c(orange, "white"))

pie(
	c(sum(go_comp_coverage), sum(!go_comp_coverage)), 
	labels=c(sum(go_comp_coverage), sum(!go_comp_coverage)), 
	main="GO biological component",
	col=c(pink, "white"))

pie(
	c(sum(go_chea_coverage), sum(!go_chea_coverage)), 
	labels=c(sum(go_chea_coverage), sum(!go_chea_coverage)), 
	main="ChEA",
	col=c(brown, "white"))

dev.off()


# Enrichment barplots
# -----------------------------------------------------
pdf("exploratory/figCouplingCCLE/emt_signature/up_enrichment/bar_plot.pdf", height=2.5)
par(mar=c(5.1, 10, 4.1, 2.1), mfrow=c(1, 3))

barplot(
	rev(-log10(pos_go_process[1:n, "Adjusted P-value"])),
	main="GO biological process",
	names.arg=rev(pos_go_process$Term[1:n]),
	cex.names=0.5,
	space=0.7,
	border=NA,
	col=orange,
	las=2,
	horiz=TRUE,
	xlab=expression(paste(-log["10"], "p (BH)"))
)
abline(v=-log10(0.05), col="grey")



barplot(
	rev(-log10(pos_go_comp[1:n, "Adjusted P-value"])),
	main="GO cellular component",
	names.arg=rev(pos_go_comp$Term[1:n]),
	cex.names=0.5,
	space=0.7,
	border=NA,
	col=pink,
	las=2,
	horiz=TRUE,
	xlab=expression(paste(-log["10"], "p (BH)"))

)
abline(v=-log10(0.05), col="grey")


barplot(
	rev(-log10(pos_chea[1:n, "Adjusted P-value"])),
	main="Upstream TFs",
	names.arg=rev(pos_chea$Term[1:n]),
	cex.names=0.5,
	space=0.7,
	border=NA,
	col=brown,
	las=2,
	horiz=TRUE,
	xlab=expression(paste(-log["10"], "p (BH)"))

)
abline(v=-log10(0.05), col="grey")

dev.off()


# -------------------------------------------------------------

n_colors = 100
color_scale = rev(colorRampPalette(brewer.pal(9, "RdBu"))(n_colors))

row_quant = apply(prob_low_cor_all, 1, function(row) {
	q = quantile(abs(row), 0.95, na.rm=TRUE)
	return(q)
})

col_quant = apply(prob_low_cor_all, 2, function(col) {
	q = quantile(abs(col), 0.95, na.rm=TRUE)
	return(q)
})

# Bimodal and high entropy proteins
high_entropy = which(
	delta_bic > bic_cut  &
	class_stat$min_entropy > quantile(class_stat$min_entropy[delta_bic > bic_cut], entropy_quant, na.rm=TRUE)
)

high_entropy_mrna = which(
	delta_bic_mrna > bic_cut  &
	class_stat_ccle$min_entropy > quantile(class_stat_ccle$min_entropy[delta_bic_mrna > bic_cut], entropy_quant, na.rm=TRUE)
)


# Subset of the cross coupling matrix
# cor_mat = prob_low_cor_all[
# 	high_entropy,
# 	delta_bic_mrna > 2
# ]

cor_mat = prob_low_cor_all[
	high_entropy,
	high_entropy_mrna
]

cor_mat[is.na(cor_mat)] = 0.0  # for weighting missing values

row_mean = apply(abs(cor_mat), 1, median, na.rm=TRUE)
col_mean = apply(abs(cor_mat), 2, median, na.rm=TRUE)

mat = cor_mat[
	row_mean > quantile(row_mean, 0.75),
	col_mean > quantile(col_mean, 0.99)
]

pdf("exploratory/figCouplingCCLE/protein_mRNA_cross_coupling_filtered_heatmap.pdf", height=5)
heatmap.2(
	# prob_low_cor,
	mat,
	distfun=function(x) {
		# x[is.na(x)] = 0.0  # missing values assumed to indicate no 
		# dmat = dist(x, method="minkowski", p=1.5)
		# dmat = dist(x, method="minkowski", p=2.0)
		dmat = dist(x)
		dmat[is.na(dmat)] = max(dmat, na.rm=TRUE) + 1
		return(dmat)
	},
	breaks=seq(-1, 1, length=n_colors + 1),  # full range color scheme
	trace="none",
	col=color_scale,
	cexCol=0.3,
	cexRow=0.3,
	# cexRow=0.25,
	# ColSideColors=col_row,
	main="Row-column median quantile filter",
	key.title="",
	key.xlab="Coupling",
	key.ylab="",
	xlab="mRNA", ylab="Protein",
)
dev.off()

# Construct aligned 

# delta_bic_sort = sort(delta_bic, decreasing=FALSE, index.return=TRUE)
delta_bic_order = match(names(delta_bic_sort$x), colnames(x))

norm_loglik_ccle = normLoglik(t(ccle_mat))

# par(mfrow=c(2, 5))
# plotHist(x, delta_bic_order[1:10])


# sum(delta_bic > 2, na.rm=TRUE)

# bic_threshold = 2
# bic_sig = which(delta_bic > 2)
# bic_sig_equal = which(delta_bic > 2 & em_lambda >= 1/3 & em_lambda <= 2/3)


# Classes of genes, separated by Delta BIC, tissue entropy, and prior mixture.
non_bimodal = which(
	delta_bic <= bic_cut |
	is.na(delta_bic)
)

# hist(class_stat$min_entropy[delta_bic > bic_cut], breaks=50, border=NA, col="black")
# hist(class_stat$total_entropy[delta_bic > bic_cut], breaks=50, border=NA, col="black")

# Bimodal subclasses
# Low entropy => tissue-specific
low_entropy = which(
	delta_bic > bic_cut &
	class_stat$min_entropy <= quantile(class_stat$min_entropy[delta_bic > bic_cut], entropy_quant, na.rm=TRUE)
)

high_entropy = which(
	delta_bic > bic_cut  &
	class_stat$min_entropy > quantile(class_stat$min_entropy[delta_bic > bic_cut], entropy_quant, na.rm=TRUE)
)

# hist(class_stat)

high_entropy_mid = which(
	delta_bic > bic_cut & 
	class_stat$min_entropy > quantile(class_stat$min_entropy[delta_bic > bic_cut], entropy_quant, na.rm=TRUE) &
	em$lambda >= beta & em$lambda <= (1 - beta)
)

high_entropy_rare_down = which(
	delta_bic > bic_cut & 
	class_stat$min_entropy > quantile(class_stat$min_entropy[delta_bic > bic_cut], entropy_quant, na.rm=TRUE) &
	em$lambda < beta
)

high_entropy_rare_up = which(
	delta_bic > bic_cut & 
	class_stat$min_entropy > quantile(class_stat$min_entropy[delta_bic > bic_cut], entropy_quant, na.rm=TRUE) &
	em$lambda > (1 - beta)
)




# Collect table of statistics 
bimodal_prot_tab = data.frame(
	RPPA=names(delta_bic),
	delta_BIC=delta_bic,
	mixture_coef_low=em$lambda,
	min_tissue_entropy=class_stat$min_entropy,
	annot=NA
)

bimodal_prot_tab$annot[non_bimodal] = "non-bimodal"
bimodal_prot_tab$annot[low_entropy] = "bimodal_low-entropy"
bimodal_prot_tab$annot[high_entropy_mid] = "bimodal_high-entropy_common"
bimodal_prot_tab$annot[high_entropy_rare_down] = "bimodal_high-entropy_rare-down"
bimodal_prot_tab$annot[high_entropy_rare_up] = "bimodal_high-entropy_rare-up"

# Counts of annotations
table(bimodal_prot_tab$annot)

bimodal_prot_tab = bimodal_prot_tab[order(bimodal_prot_tab$delta_BIC, decreasing=TRUE), ]

prot_annot = makeProteinTable(bimodal_prot_tab$RPPA)

bimodal_prot_tab$gene_symbol = as.character(prot_annot$uniprot_name)


write.table(bimodal_prot_tab,
	file="data/bimodal/bimodal_prot2.csv",
	sep=",",
	quote=FALSE,
	row.names=FALSE
)


# all(rownames(class_stat) == colnames(x))

# Write lists of Uniprot gene symbols to files.
write.table(
	t(data.frame(proteins$uniprot_name[non_bimodal])),
	file="exploratory/figLayerDeconv/geneLists/non_bimodal.tsv",
	row.names=FALSE, col.names=FALSE,
	quote=FALSE
)


write.table(
	t(data.frame(proteins$uniprot_name[low_entropy])),
	file="exploratory/figLayerDeconv/geneLists/low_entropy.tsv",
	row.names=FALSE, col.names=FALSE,
	quote=FALSE
)

write.table(
	t(data.frame(proteins$uniprot_name[high_entropy_mid])),
	file="exploratory/figLayerDeconv/geneLists/high_entropy_mid.tsv",
	row.names=FALSE, col.names=FALSE,
	quote=FALSE
)

write.table(
	t(data.frame(proteins$uniprot_name[high_entropy_rare_down])),
	file="exploratory/figLayerDeconv/geneLists/high_entropy_rare_down.tsv",
	row.names=FALSE, col.names=FALSE,
	quote=FALSE
)

write.table(
	t(data.frame(proteins$uniprot_name[high_entropy_rare_up])),
	file="exploratory/figLayerDeconv/geneLists/high_entropy_rare_up.tsv",
	row.names=FALSE, col.names=FALSE,
	quote=FALSE
)


# non_bimodal

# delta_bic[non_bimodal]

i = which(colnames(x) == "JAB1")
i = which(colnames(x) == "P27")
# class_stat$total_entropy

hist(x[,low_entropy[4]], breaks=100, col="black", border=NA)

hist(x[,i], breaks=100, col="black", border=NA)

# data.frame(val=x[,i], p=prob_low[,i])


colors = c(rgb(240,163,255,maxColorValue=255),rgb(0,117,220,maxColorValue=255),rgb(153,63,0,maxColorValue=255),rgb(76,0,92,maxColorValue=255),rgb(25,25,25,maxColorValue=255),rgb(0,92,49,maxColorValue=255),rgb(43,206,72,maxColorValue=255),rgb(255,204,153,maxColorValue=255),rgb(128,128,128,maxColorValue=255),rgb(148,255,181,maxColorValue=255),rgb(143,124,0,maxColorValue=255),rgb(157,204,0,maxColorValue=255),rgb(194,0,136,maxColorValue=255),rgb(0,51,128,maxColorValue=255),rgb(255,164,5,maxColorValue=255),rgb(255,168,187,maxColorValue=255),rgb(66,102,0,maxColorValue=255),rgb(255,0,16,maxColorValue=255),rgb(94,241,242,maxColorValue=255),rgb(0,153,143,maxColorValue=255),rgb(224,255,102,maxColorValue=255),rgb(116,10,255,maxColorValue=255),rgb(153,0,0,maxColorValue=255),rgb(255,255,128,maxColorValue=255),rgb(255,255,0,maxColorValue=255),rgb(255,80,5,maxColorValue=255))

n_colors=100
# color_scale = colorRampPalette(brewer.pal(9, "PRGn"))(n_colors)
color_scale = rev(colorRampPalette(brewer.pal(9, "RdBu"))(n_colors))
col_row = colors[match(ccle$Site.Primary, unique(ccle$Site.Primary[!ccle$inferred]))]
# col_row[ccle$inferred] = "white"


pdf("exploratory/figLayerDeconv/post_prob_high_entropy_rare.pdf")
a = t(1 - prob_low[,high_entropy_rare_up])
# a[a < 0.9] = NA  # Clamp uncertain values
a[a < 0.9] = 0.0  # Clamp uncertain values

b = t(prob_low[,high_entropy_rare_down])
# b[b < 0.9] = NA
b[b < 0.9] = 0.0
b = -b  # negative encoding of "down" probabilities

c = rbind(a, b)
hv = heatmap.2(
	c,
	main=paste0(nrow(c), " bimodal, high tissue entropy, rare."),
	# the distance function minimally replaces the pairwise distances that 
	# cannot be computed with the maximum distance.
	distfun=function(x) {
		# x[is.na(x)] = 0.0  # missing values assumed to indicate no 
		# dmat = dist(x, method="minkowski", p=1.5)
		dmat = dist(abs(x), method="minkowski", p=2.0)
		dmat[is.na(dmat)] = max(dmat, na.rm=TRUE) + 1
		return(dmat)
	},
	trace="none",
	col=color_scale,
	cexCol=0.07,
	cexRow=0.5,
	# cexRow=0.25,
	ColSideColors=col_row,
	key.title="Posterior probability",
	key.xlab="",
	key.ylab=""
)
dev.off()


pdf("exploratory/figLayerDeconv/post_prob_high_entropy_mid.pdf")
a = t(1 - prob_low[,high_entropy_mid])
a = 2 * a - 1  # probability -> [-1, 1]
a[a > -0.8 & a < 0.8] = 0.0  # clamping, not probability!
# a[a > -0.8 & a < 0.8] = NA  # clamping, not probability!

hv = heatmap.2(
	a,
	main=paste0(nrow(a), " bimodal, high tissue entropy, common."),
	# the distance function minimally replaces the pairwise distances that 
	# cannot be computed with the maximum distance.
	distfun=function(x) {
		dmat = dist(x, method="minkowski", p=2.0)
		dmat[is.na(dmat)] = max(dmat, na.rm=TRUE) + 1
		return(dmat)
	},
	trace="none",
	col=color_scale,
	cexCol=0.07,
	cexRow=0.5,
	# cexRow=0.25,
	ColSideColors=col_row,
	key.title="Posterior probability",
	key.xlab="",
	key.ylab=""
)
dev.off()


# # Low tissue entropy
# a = t(1 - prob_low[,low_entropy])
# a = 2 * a - 1  # probability -> [-1, 1]
# # a[a > -0.8 & a < 0.8] = 0.0  # clamping, not probability!
# # a[a > -0.8 & a < 0.8] = NA  # clamping, not probability!

# hv = heatmap.2(
# 	a,
# 	main=paste0(nrow(a), " bimodal, low tissue entropy."),
# 	# the distance function minimally replaces the pairwise distances that 
# 	# cannot be computed with the maximum distance.
# 	distfun=function(x) {
# 		dmat = dist(x, method="minkowski", p=1.5)
# 		dmat[is.na(dmat)] = max(dmat, na.rm=TRUE) + 1
# 		return(dmat)
# 	},
# 	trace="none",
# 	col=color_scale,
# 	cexCol=0.07,
# 	cexRow=0.5,
# 	# cexRow=0.25,
# 	ColSideColors=col_row,
# 	key.title="Posterior probability",
# 	key.xlab="",
# 	key.ylab=""
# )

# Distribution of average, variance, ... 
# -------------------------------------
pdf("exploratory/figLayerDeconv/univar_moments_distr.pdf", width=3, height=8)
colors = brewer.pal(12, "Paired")
nbreaks = 100
par(mfrow=c(5, 1))
hist(prot_means, main="Average", xlab="Averages", border=NA, col=colors[2], breaks=nbreaks)
hist(prot_var, main="Variance", xlab="Variances", border=NA, col=colors[4], breaks=nbreaks)
hist(prot_skew, main="Skewness", xlab="Skewness", border=NA, col=colors[6], breaks=nbreaks)
hist(prot_kurt, main="Kurtosis", xlab="Kurtosis", border=NA, col=colors[8], breaks=nbreaks)
hist(-em_loglik, main="Log-likelihood of two-component Gaussian mixture model", xlab="log-likelihood", border=NA, col=colors[10], breaks=nbreaks, cex.main=0.7)
dev.off()


pdf("exploratory/figLayerDeconv/univar_moments_ranked.pdf", width=16, height=8)
n = 10
m = 10
par(
	mfcol=c(n, m),
	mar=c(2, 2, 1, 1)
)
# n lowest means
plotHist(x, prot_means_sort$ix[ncol(x):((ncol(x) - n + 1))], col=colors[1])
# n highest means
plotHist(x, prot_means_sort$ix[1:n], col=colors[2])

# Lowest var
plotHist(x, prot_var_sort$ix[ncol(x):((ncol(x) - n + 1))], col=colors[3])
# highest var
plotHist(x, prot_var_sort$ix[1:n], col=colors[4])

# Lowest skew
plotHist(x, prot_skew_sort$ix[ncol(x):((ncol(x) - n + 1))], col=colors[5])
# Highest skew
plotHist(x, prot_skew_sort$ix[1:n], col=colors[6])

# Lowest kurt
plotHist(x, prot_kurt_sort$ix[ncol(x):((ncol(x) - n + 1))], col=colors[7])
# Highest kurt
plotHist(x, prot_kurt_sort$ix[1:n], col=colors[8])

# Unfavorable lok likelihood, highest
# plotHist(x, em_loglik_sort$ix[1:n], col=colors[9])
plotHist(x, em_loglik_order[1:n], col=colors[9])
# Favorable log likelihood, lowest
plotHist(x, em_loglik_order[length(em_loglik_sort$ix):((length(em_loglik_sort$ix) - n + 1))], col=colors[10])
dev.off()


# Analysis of the most bimodal proteins identified
# --------------------------------------------------------------------

# Get ids of significantly bimodal proteins
# n = 100

# Based on EM loglikelihood (not theoretically valid)
# ids = em_loglik_order[length(em_loglik_sort$ix):((length(em_loglik_sort$ix) - n + 1))]

# Based on log likelihood difference.
# ids = diff_loglik_order[length(diff_loglik_sort$ix):((length(diff_loglik_sort$ix) - n + 1))]
# ids = diff_loglik_order[1:n]

# ids = delta_bic_order[1:n]
# ids = delta_bic_order[1:50]
# Protein bimodality order based on BIC. Cutoff at 10 indicating proteins which are strongy bimodal.
# ids = delta_bic_order[1:sum(delta_bic > 5, na.rm=TRUE)]
# ids = delta_bic_order[1:sum(delta_bic > 5, na.rm=TRUE)]
ids = delta_bic_order[1:sum(delta_bic > 2, na.rm=TRUE)]

# ids = ids[1:26]


# Histograms of all top bimodal proteins
# par(mfrow=c(5, 10))
par(mfrow=c(5, 5))

plotHist(x, ids[1:20], col="black")

# pdf("exploratory/figLayerDeconv/bimodal_log_likelihood.pdf", height=3, width=4)
# plot(-em_loglik_sort$x, type="n", bty="n",
# 	main="Drop-off in bimodal score",
# 	xlab="Order",
# 	ylab="log-likelihood"
# 	)
# abline(h=0, col="grey")
# lines(-rev(em_loglik_sort$x), type="l", lwd=2)
# dev.off()

pdf("exploratory/figLayerDeconv/bimodal_delta_bic.pdf", height=3, width=7)
par(mfrow=c(1, 2))

# Plots of the E-Cadherin mixture model fit
# ------------------------------------------------
i = match("ECADHERIN", colnames(x))
em = em_results[[i]]

# plotHist(x, i, col=colors[9], xlim=range(x[,i], na.rm=TRUE), freq=FALSE)

# pdf("exploratory/figLayerDeconv/hist_em_fit_ecadherin.pdf", height=2.5, width=3)
hist(x[,i], border=NA, col="grey", breaks=50, freq=FALSE, ylim=c(0, 0.65),
	xlab="Protein expression",
	main=colnames(x)[i])

blue = rgb(55, 126, 184, maxColorValue=255)
red = rgb(228, 26, 28, maxColorValue=255)

x_view = seq(-10, 10, length.out=300)
lines(x_view, dnorm(x_view, em$mu[1], em$sigma[1]), col=blue, lwd=2)
lines(x_view, dnorm(x_view, em$mu[2], em$sigma[2]), col=red, lwd=2)


# Rank plot of the BIC metric
# ---------------------------------
# k = 400
plot(delta_bic_sort$x, type="n", bty="n",
	main="Support for 2 components",
	xlab="Protein order",
	# ylab="Delta BIC"
	ylab=expression(paste(Delta, "BIC"))
	)
# abline(h=0, col="grey")
# lines(-rev(delta_bic_sort$x[1:k]), type="l", lwd=2)
order_cutoff = max(which(delta_bic_sort$x > 10))
abline(v=order_cutoff, col="grey")

lines(delta_bic_sort$x, type="l", lwd=2)

dev.off()

# CCLE delta BIC 
pdf("exploratory/figLayerDeconv/bimodal_delta_bic_CCLE.pdf", height=3, width=3.5)
delta_bic_sort_ccle = sort(delta_bic_ccle, decreasing=TRUE, index.return=TRUE)
plot(delta_bic_sort_ccle$x,
	type="n", bty="n",
	main="Support for 2 components",
	xlab="Transcript order",
	# ylab="Delta BIC"
	ylab=expression(paste(Delta, "BIC"))
	)
order_cutoff = max(which(delta_bic_sort_ccle$x > 10))
abline(v=order_cutoff, col="grey")

lines(delta_bic_sort_ccle$x, type="l", lwd=2)
dev.off()



plot(density(delta_bic, na.rm=TRUE))

# E-cadherin example
# Separate histograms by tissue
# --------------------------------------------------
# Colors of maximum separability adapted from the colour alphabet.
# https://eleanormaclure.files.wordpress.com/2011/03/colour-coding.pdf
colors = c(rgb(240,163,255,maxColorValue=255),rgb(0,117,220,maxColorValue=255),rgb(153,63,0,maxColorValue=255),rgb(76,0,92,maxColorValue=255),rgb(25,25,25,maxColorValue=255),rgb(0,92,49,maxColorValue=255),rgb(43,206,72,maxColorValue=255),rgb(255,204,153,maxColorValue=255),rgb(128,128,128,maxColorValue=255),rgb(148,255,181,maxColorValue=255),rgb(143,124,0,maxColorValue=255),rgb(157,204,0,maxColorValue=255),rgb(194,0,136,maxColorValue=255),rgb(0,51,128,maxColorValue=255),rgb(255,164,5,maxColorValue=255),rgb(255,168,187,maxColorValue=255),rgb(66,102,0,maxColorValue=255),rgb(255,0,16,maxColorValue=255),rgb(94,241,242,maxColorValue=255),rgb(0,153,143,maxColorValue=255),rgb(224,255,102,maxColorValue=255),rgb(116,10,255,maxColorValue=255),rgb(153,0,0,maxColorValue=255),rgb(255,255,128,maxColorValue=255),rgb(255,255,0,maxColorValue=255),rgb(255,80,5,maxColorValue=255))

n = 10  # plot nth most bimodal histograms by tissue type

# Catalogue of proteins outside top n, to be plotted "manually"
i = which(colnames(x) == "VEGFR2")
i = which(colnames(x) == "BIM")

for (k in 1:n) {
	i = ids[k]  # column id

	protein = colnames(x)[i]

	pdf(paste0("exploratory/figLayerDeconv/bimodal_hist", protein, ".pdf"), width=7, height=2.5)

	par(mfrow=c(4, 5),
		mar=c(2, 2, 1, 1))
	plotHist(x, i, col="black")  # histogram of all RPPA data for protein

	# Break protein measurements by tissue type
	for (tissue in unique(ccle$Site.Primary)) {
		rows = which(ccle$Site.Primary == tissue)  # including inferred tissue type
		# rows = which(ccle$Site.Primary == tissue & !ccle$inferred)  # only directly mapped tissue types.

		plotHist(x[rows,], i, main=tissue,
			xlim=range(x, na.rm=TRUE), 
			breaks = seq(range(x, na.rm=TRUE)[1], range(x, na.rm=TRUE)[2], length.out=100),
			col=colors[which(tissue == unique(ccle$Site.Primary[!ccle$inferred]))])
	}

	dev.off()
}


# Mixture model assignments and posterior correlations
# ------------------------------------

# ids = delta_bic_order[1:sum(delta_bic > 2, na.rm=TRUE)]
# ids = which(delta_bic > 2)
ids = high_entropy
subx = x[,ids]

# Associated mixture models for the selected proteins
sub_em_results = em_results[ids] 

prob_low = matrix(NA, nrow=nrow(subx), ncol=ncol(subx))
colnames(prob_low) = colnames(subx)
rownames(prob_low) = rownames(subx)
for (i in 1:nrow(subx)) {
	for (j in 1:ncol(subx)) {
		prob_low[i, j] = probLowModel(
			subx[i, j],
			lambda=sub_em_results[[j]]$lambda,
			mu=sub_em_results[[j]]$mu, 
			sigma=sub_em_results[[j]]$sigma
		)
	}
}

# Calculate bimodal protein coupling (correlation between posterior probabilities).
# prob_cor_low = cor(prob_low, use="pairwise.complete.obs")
prob_cor_low = cor(prob_low, use="pairwise.complete.obs", method="spearman")

# Test that the posterior calculations are identical to the R mixtools package.
# i = 10
# plot(em_results[ids][[i]]$posterior[,1], prob_low[!is.na(prob_low[,i]), i])

# Filter the posterior probabilities based on missing values.
missing_tolerance = 0.8  # the fraction of observations required to not be missing for inclusion.
nprob = apply(prob_low, 2, function(row) {
	return(sum(!is.na(row)))
})
prob_low_filter = prob_low[, nprob > nrow(prob_low) * missing_tolerance]

prob_cor_low_filter = cor(prob_low_filter, use="pairwise.complete.obs")





# Heatmap plot of posterior probabilities of the Gaussian mixture models.
pdf(paste0("exploratory/figLayerDeconv/post_prob_heatmap", length(ids), ".pdf"), height=5, width=8)
n_colors=100
# color_scale = colorRampPalette(brewer.pal(9, "PRGn"))(n_colors)
color_scale = rev(colorRampPalette(brewer.pal(9, "RdBu"))(n_colors))
col_row = colors[match(ccle$Site.Primary, unique(ccle$Site.Primary[!ccle$inferred]))]
# col_row[ccle$inferred] = "white"

hv = heatmap.2(t(1 - prob_low),
	main=paste0(ncol(prob_low), " bimodal proteins. tol=", missing_tolerance),
	# the distance function minimally replaces the pairwise distances that 
	# cannot be computed with the maximum distance.
	distfun=function(x) {
		dmat = dist(x, method="minkowski", p=1.5)
		dmat[is.na(dmat)] = max(dmat, na.rm=TRUE) + 1
		return(dmat)
	},
	trace="none",
	col=color_scale,
	cexCol=0.07,
	cexRow=0.5,
	# cexRow=0.25,
	ColSideColors=col_row,
	key.title="Posterior probability",
	key.xlab="",
	key.ylab=""
)
dev.off()

# Get the plotting protein order
prot_order = colnames(prob_low_filter)[hv$rowInd]

prot_ids = match(prot_order, rownames(class_stat))

class_stat[prot_ids,]
# barplot(c(class_stat$high_entropy[prot_ids], class_stat$low_entropy[prot_ids]))

# Image of tissue entropy matching the posterior probability heatmap.
# White indicates high entropy, black indicates low entropy.
n_colors = 100
pdf(paste0("exploratory/figLayerDeconv/post_prob_heatmap", length(ids), "_entropy.pdf"), height=5, width=1.5)
image(t(as.matrix(class_stat[c("low_entropy", "high_entropy")][prot_ids,])),
	axes=FALSE,
	box=FALSE,
	col=gray.colors(n_colors, start=0.1, end=0.95)
	# col=colorRampPalette(brewer.pal(9, "Greens"))(n_colors)
	# col=rev(colorRampPalette(brewer.pal(9, "PRGn"))(n_colors))
)
dev.off()



# color_scale = colorRampPalette(c("grey100", "grey0"))(n_colors)


pdf(paste0("exploratory/figLayerDeconv/bimodal_filter", length(ids), "_cor_heatmap.pdf"))
color_scale = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(n_colors))
heatmap.2(
	prob_cor_low_filter,
	main="Posterior probability correlation",
	trace="none",
	col=color_scale,
	distfun=function(x) {
		dmat = dist(abs(x))
		dmat[is.na(dmat)] = max(dmat, na.rm=TRUE) + 1
		# dist(abs(x))
		return(dmat)
	},
	# RowSideColors=colors[match(ccle_match_info$Site.Primary, unique(ccle_match_info$Site.Primary))],
	# ColSideColors=colors[match(ccle_match_info$Site.Primary, unique(ccle_match_info$Site.Primary))],
	cexCol=0.3,
	cexRow=0.3,
	breaks=seq(-1, 1, length.out=n_colors+1)
)
dev.off()

pdf(paste0("exploratory/figLayerDeconv/bimodal_all", length(ids), "_cor_heatmap.pdf"))
color_scale = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(n_colors))
heatmap.2(
	prob_cor_low,
	main="Posterior probability correlation",
	trace="none",
	col=color_scale,
	distfun=function(x) {
		dmat = dist(abs(x))
		dmat[is.na(dmat)] = max(dmat, na.rm=TRUE) + 1
		# dist(abs(x))
		return(dmat)
	},
	# RowSideColors=colors[match(ccle_match_info$Site.Primary, unique(ccle_match_info$Site.Primary))],
	# ColSideColors=colors[match(ccle_match_info$Site.Primary, unique(ccle_match_info$Site.Primary))],
	cexCol=0.3,
	cexRow=0.3,
	breaks=seq(-1, 1, length.out=n_colors+1)
)
dev.off()




# Graph analysis
# ------------------------------------------------------------
# see
# https://rpubs.com/kateto/netviz
# for different network vizualizations

# Merge function for weighted graphs
mergeGraph = function(g1, g2) {
  e1 = get.data.frame(g1, what="edges")
  e2 = get.data.frame(g2, what="edges")
  e = merge(e1, e2, by=c("from", "to"), all=TRUE)
  newe = data.frame(e[,c("from", "to"), drop=FALSE],
          weight=rowSums(e[, c("weight.x", "weight.y")], na.rm=TRUE))
  graph.data.frame(newe, directed=is.directed(g1))
}


# Filter based on bimodal and high entropy
ids = high_entropy
subx = x[,ids]

# Associated mixture models for the selected proteins
sub_em_results = em$results[ids] 

# Calculate posterior probability of low bin assignment
prob_low = matrix(NA, nrow=nrow(subx), ncol=ncol(subx))
colnames(prob_low) = colnames(subx)
rownames(prob_low) = rownames(subx)
for (i in 1:nrow(subx)) {
	for (j in 1:ncol(subx)) {
		prob_low[i, j] = probLowModel(
			subx[i, j],
			lambda=sub_em_results[[j]]$lambda,
			mu=sub_em_results[[j]]$mu, 
			sigma=sub_em_results[[j]]$sigma
		)
	}
}

# Calculate bimodal protein coupling (correlation between posterior probabilities).
# prob_cor_low = cor(prob_low, use="pairwise.complete.obs")
prob_cor_low = cor(prob_low, use="pairwise.complete.obs", method="spearman")

sort(rownames(prob_cor_low))
prob_cor_low[which(rownames(prob_cor_low) == "ECADHERIN"), which(colnames(prob_cor_low) == "BETACATENIN")]
prob_cor_low[which(rownames(prob_cor_low) == "ECADHERIN"), which(colnames(prob_cor_low) == "NCADHERIN")]

prob_cor_low[which(rownames(prob_cor_low) == "ECADHERIN"), which(colnames(prob_cor_low) == "CLAUDIN7")]

# prob_cor_low[which(rownames(prob_cor_low) == "ECADHERIN"), ]

prob_cor_low_test = corAndPvalue(prob_low, use="pairwise.complete.obs", method="spearman")

# p.adjust(prob_cor_low_test$p)

pvals = prob_cor_low_test$p[lower.tri(prob_cor_low_test$p)]
r = prob_cor_low_test$cor[lower.tri(prob_cor_low_test$cor)]

sum(abs(r) > 0.3, na.rm=TRUE)
sum(abs(r) > 0.3 & p.adjust(pvals) > 0.05, na.rm=TRUE)



sum(p.adjust(pvals) < 0.05, na.rm=TRUE)
# sum((pvals) < 0.05, na.rm=TRUE)
sum(abs(r) > 0.3, na.rm=TRUE)
sum(abs(r) > 0.3 & p.adjust(pvals) < 0.05, na.rm=TRUE)

bonferoni_p = 0.05 / (nrow(prob_cor_low)^2 / 2 - nrow(prob_cor_low))

sum(pvals < bonferoni_p, na.rm=TRUE)

# sum(p.adjust(pvals, method="bonferroni") < 0.05, na.rm=TRUE)

# sum(abs(r) > 0.3 & pvals > bonferoni_p, na.rm=TRUE)


# sum(abs(r) > 0.3 & pvals < 0.0001, na.rm=TRUE)

# sum(abs(r) > 0.3 & pvals > 0.05, na.rm=TRUE)


# sum(abs(r) > 0.3 & pvals > 0.05, na.rm=TRUE)
# sum(abs(prob_cor_low) > 0.3, na.rm=TRUE)


pdf("exploratory/figLayerDeconv/bimodal_top_cor_network.pdf")
qgraph(prob_cor_low,
	layout="spring",
	minimum=0.15,  # minimum weights
	negCol=rgb(55, 126, 184, maxColorValue=255),
	posCol=rgb(228, 26, 28, maxColorValue=255),
	cut=0.5,  # over which sizes of edges change
	vsize=3,
	labels=colnames(prob_cor_low)
)
dev.off()

# Construct adjecency matrix
# beta = 2  # soft threshold exponent

# Cutoff
w = prob_cor_low

# w[abs(w) < 0.25] = 0.0  # hard threshold for correlations
w[abs(w) < 0.3] = 0.0  # hard threshold for correlations

# FDR threshol
bonferoni_p = 0.05 / (nrow(prob_cor_low)^2 / 2 - nrow(prob_cor_low))
w[abs(w) < 0.3 | prob_cor_low_test$p > bonferoni_p] = 0.0  # hard threshold for correlations

sum(abs(w) > 0.2, na.rm=TRUE)/2


# w[abs(w) < 0.4] = 0.0  # hard threshold for correlations
# w[abs(w) < 0.5] = 0.0  # hard threshold for correlations

# w = w^beta  # soft threshold
w[is.na(w)] = 0.0  # encoding for non-existing edges
diag(w) = 0.0  # no self edges

w_sign = w  # signed adjacency matrix
w = abs(w)

# Construct graph object from adjacency matrix
graph = graph.adjacency(w, mode="lower", weighted=TRUE)  # unsigned graph for calculating communities
sgraph = graph.adjacency(w_sign, mode="lower", weighted=TRUE)  # signed graph

# Set graph attribute; EM mixture coefficients from fit.
# V(graph)$lambda = em_lambda[match(names(V(graph)), colnames(x))]
# V(sgraph)$lambda = em_lambda[match(names(V(sgraph)), colnames(x))]

# plot(log(degree.distribution(graph)), log="x")
# 
# plot(degree(graph)

# Find community clusters using the Girvan-Newman algorithm.
# clust = cluster_fast_greedy(graph)

# Girvan-Newman, slower algorithm. Results in ~117 cluster.
# clust = cluster_edge_betweenness(graph)  

# Random walk model
# clust = cluster_walktrap(graph)

# Spin-class model and simulated annealing.
# clust = cluster_spinglass(graph)

# Leading eigenvector.
clust = cluster_leading_eigen(graph)

# write.graph(sgraph, file="data/layerDeconv/networks/prot_adj2.gml", format="gml")
write.graph(sgraph, file="data/layerDeconv/networks/prot_adj_FDR05.gml", format="gml")



# Members of each cluster
clust_sizes = rep(NA, length(clust))
for (i in 1:length(clust)) {
	clust_sizes[i] = length(clust[[i]])
}

clust_include = which(clust_sizes >= 2)

# Concatenate communities into single graph
com_graph = induced_subgraph(sgraph, clust[[clust_include[1]]])

for (i in 2:length(clust_include)) {
	com_graph = mergeGraph(com_graph,
		induced_subgraph(sgraph, clust[[clust_include[i]]])
	)
}

# Annotate combined graph with EM mixtures (labmda values)
V(com_graph)$lambda = em$lambda[match(names(V(com_graph)), colnames(x))]

# write.graph(com_graph, file="data/layerDeconv/networks/prot_adj_community2.gml", format="gml")
write.graph(com_graph, file="data/layerDeconv/networks/prot_adj_community_FDR05.gml", format="gml")


# sub_graphs = list()
# for (i in clust_include) {

# 	clust_ids = match(clust[[i]], colnames(w_sign))

# 	sub_graphs[[i]] = graph.adjacency(w_sign[clust_ids, clust_ids], mode="lower", weighted=TRUE)

# 	# write sub graph
# 	write.graph(sub_graphs[[i]], file=paste0("exploratory/figLayerDeconv/networks/prot_adj_", i, ".gml"), format="gml")
# }

# Calculate average between cluster interaction
w_between = matrix(NA, nrow=length(clust_include), ncol=length(clust_include))
for (i in clust_include) {
	clust_ids1 = match(clust[[i]], colnames(w_sign))

	for (j in clust_include) {
		clust_ids2 = match(clust[[j]], colnames(w_sign))
		w_between[i, j] = mean(w[clust_ids1, clust_ids2])
		# w_between[i, j] = max(w[clust_ids1, clust_ids2])
	}
}

diag(w_between) = 0.0

dmat = as.dist(1/w_between)
dmat[is.infinite(dmat)] = 1000

hc = hclust(dmat)
plot(hc)



i = 1
prot_ids = colnames(prob_cor_low) %in% clust[[i]]
# cor_mat = prob_cor_low[prot_ids, prot_ids]
# diag(cor_mat) = 0.0
qgraph(w_sign[prot_ids, prot_ids],
	layout="spring",
	minimum=0.1,  # minimum weights
	negCol=rgb(55, 126, 184, maxColorValue=255),
	posCol=rgb(228, 26, 28, maxColorValue=255),
	cut=0.5,  # over which sizes of edges change
	vsize=3,
	labels=colnames(cor_mat)
)


heatmap.2(
	cor_mat,
	main="Posterior probability correlation",
	trace="none",
	col=color_scale,
	distfun=function(x) {
		dmat = dist(abs(x))
		dmat[is.na(dmat)] = max(dmat, na.rm=TRUE) + 1
		# dist(abs(x))
		return(dmat)
	},
	# RowSideColors=colors[match(ccle_match_info$Site.Primary, unique(ccle_match_info$Site.Primary))],
	# ColSideColors=colors[match(ccle_match_info$Site.Primary, unique(ccle_match_info$Site.Primary))],
	cexCol=0.3,
	cexRow=0.3,
	breaks=seq(-1, 1, length.out=n_colors+1)
)


# proteins$uniprot_name[match(clust[[i]], colnames(x))]

# clust_ids = match(clust[[i]], colnames(x))
# write.table(
# 	t(data.frame(proteins$uniprot_name[clust_ids])),
# 	file="exploratory/figLayerDeconv/geneLists/clust.tsv",
# 	row.names=FALSE, col.names=FALSE,
# 	quote=FALSE
# )


# hrg_fit = hrg.fit(graph)
# dend = hrg.dendrogram(hrg_fit)
# dendPlot(hrg_fit)

# qgraph(dend)



# qgraph(w,
# 	layout="spring",
# 	minimum=0.1,  # minimum weights
# 	negCol=rgb(55, 126, 184, maxColorValue=255),
# 	posCol=rgb(228, 26, 28, maxColorValue=255),
# 	cut=0.5,  # over which sizes of edges change
# 	vsize=3,
# 	labels=colnames(cor_mat)
# )



# Calculate network centrality.
cen = subgraph_centrality(graph)

# Weights are not taken into account.
between = betweenness(graph, directed=FALSE)

# centrality = evcent(graph)
# sort(centrality$vector)


pdf("exploratory/figLayerDeconv/network_measures.pdf", width=10, height=4)
par(mfrow=c(1, 3))

dotchart(rev(sort(between, decreasing=TRUE)[1:20]), xlim=c(0, max(between)),
	# lcolor="white",
	pch=16,
	xlab="Betweenness")

dotchart(rev(sort(cen, decreasing=TRUE)[1:20]), xlim=c(0, max(cen)),
	pch=16,
	xlab="Centrality")

dotchart(rev(sort(degree(graph), decreasing=TRUE)[1:20]), xlim=c(0, max(degree(graph))),
	pch=16,
	xlab="Degree")
dev.off()



# plot(graph)

cliq = cliques(graph)
cliq = maximal.cliques(graph)



# Hive plot test
# ---------------------------
# hive = adj2HPD(w, type="2D")
hive = adj2HPD(w)

hive2 = mineHPD(hive, option = "rad <- tot.edge.count")

hive3 = mineHPD(hive2, option = "axis <- source.man.sink")

# hive = mineHPD(hive)
hive3$edges$weight <- sqrt(hive3$edges$weight)*0.5
hive3$nodes$size <- 0.5

plotHive(hive3)


hive4 <- mineHPD(hive3, option = "remove zero edge")

plotHive(graph)



qgraph(prob_cor_low,
	layout="spring",
	minimum=0.25,  # minimum weights
	negCol=rgb(55, 126, 184, maxColorValue=255),
	posCol=rgb(228, 26, 28, maxColorValue=255),
	cut=0.5,  # over which sizes of edges change
	vsize=3,
	labels=colnames(prob_cor_low)
)


sort(prob_cor_low[rownames(prob_cor_low) == "ECADHERIN", ])
sort(prob_cor_low[rownames(prob_cor_low) == "CLAUDIN7", ])
sort(prob_cor_low[rownames(prob_cor_low) == "PCADHERIN", ])

# Scatterplots of particual protein measurements colored by tissue of origin.
# ---------------------------------------------
protein1 = "ECADHERIN"
protein2 = "BETACATENIN"

protein1 = "VEGFR2"
protein2 = "BETACATENIN"

protein1 = "ECADHERIN"
protein2 = "CLAUDIN7"

protein1 = "ECADHERIN"
protein2 = "PCADHERIN"

protein1 = "ECADHERIN"
protein2 = "SYK"

protein1 = "ECADHERIN"
protein2 = "CAVEOLIN1"

protein1 = "ECADHERIN"
protein2 = "AXL"

protein1 = "ECADHERIN"
protein2 = "PLK1"

protein1 = "ECADHERIN"
protein2 = "MDM2_pS166"


protein1 = "ECADHERIN"
protein2 = "INPP4B"

protein1 = "ECADHERIN"
protein2 = "CD49B"

protein1 = "ECADHERIN"
protein2 = "EGFR"

protein1 = "ECADHERIN"
protein2 = "BETACATENIN"

protein1 = "ECADHERIN"
protein2 = "EGFR"

protein1 = "ECADHERIN"
protein2 = "JAGGED1"

protein1 = "ECADHERIN"
protein2 = "EGFR_pY1068"

protein1 = "ECADHERIN"
protein2 = "TUBERIN_pT1462"


protein1 = "CLAUDIN7"
protein2 = "SYK"

protein1 = "YAP_pS127"
protein2 = "BETACATENIN"

protein1 = "BIM"
protein2 = "BETACATENIN"

protein1 = "PAI1"
protein2 = "CAVEOLIN1"

protein1 = "YAP_pS127"
protein2 = "BIM"

protein1 = "ANNEXIN1"
protein2 = "CAVEOLIN1"

protein1 = "P53"
protein2 = "P16INK4A"

# pdf(paste0("exploratory/figLayerDeconv/scatter_", protein1, "_", protein2, ".pdf"), width=8, height=4.25)
svg(paste0("exploratory/figLayerDeconv/scatter_", protein1, "_", protein2, ".svg"), width=8, height=4.25)

# setEPS()
# postscript(paste0("exploratory/figLayerDeconv/scatter_", protein1, "_", protein2, ".eps"), width=8, height=4.25)

i = which(colnames(x) == protein1)
j = which(colnames(x) == protein2)

tissues = factor(ccle$Site.Primary)
tissues[ccle$inferred] = NA  # remove inferred tissues

par(mfrow=c(1, 2))

tissue_map = match(tissues, unique(ccle$Site.Primary[!ccle$inferred]))

# col = colors[as.numeric(tissues)]
col = colors[tissue_map]
col[is.na(col)] = rgb(0,0,0,0.05)

# plot(x[,i], x[,j], col=col, xlab=protein1, ylab=protein2, type="n")
plot(x[,i], x[,j], col=col, xlab=protein1, ylab=protein2, pch=16)
# points(x[,i], x[,j], col=col, xlab=protein1, ylab=protein2, pch=16)

plot(1, 1, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
legend("topleft", legend=unique(ccle$Site.Primary[!ccle$inferred]), col=colors, pch=16, bty="n")
dev.off()



# Figures explaining the EM algorithm for a general audience
# ------------------------------------------------------------------


i = ids[1]
j = 205


bee_data = list()
bee_data[[colnames(x)[i]]] = x[,i]
bee_data[[colnames(x)[[j]]]] = x[,j]

# par()  # adjust to size
beeswarm(bee_data, cex=0.5, bty="n", ylab="Protein expression")




# Clustering analysis of TCPA data
# ------------------------------------------------------------



# dmat = dist(x)

# hc = hclust(dmat, method="average")
# plot(hc)

# mat = x
max_clusters = 20
n_null = 20

source("lib/clustVal.r")

# Hierarchical clustering
gap_hc = gapStatHc(x, max_clusters=20, n_null)

# k-means clustering
gap_kmeans = gapStatKmeans(x, max_clusters, n_null)


pdf("exploratory/figLayerDeconv/kmeans_hc_comparison.pdf", width=5, height=5)
par(mfrow=c(2, 2))
gapPlot(gap_hc, colors[6], main="UPGMA")
gapPlot(gap_kmeans, colors[10], main="k-means")
dev.off()

# plot(gap_hc$within, type="l")
# lines(gap_kmeans$within)


# Spearman's correlations tests
spearmant = function(r, n) {
	r * sqrt((n - 2) / (1 - r^2))
}

mean(apply(x, 2, function(col) sum(!is.na(col))))

median(apply(x, 2, function(col) sum(!is.na(col))))

min(apply(x, 2, function(col) sum(!is.na(col))))



sum(apply(x, 2, function(col) sum(!is.na(col))) <= 46)

# Average
test = spearmant(0.3, 327)
# p = 1 - pt(test, 736 - 2)
# p = 1 - pt(test, 327 - 2)
p = 2 * pt(test, 327 - 2, lower.tail=FALSE)  # two-tailed

0.05 / 450


test = spearmant(0.3, 46)
# p = 1 - pt(test, 736 - 2)
# p = 1 - pt(test, 46 - 2)
p = 2*pt(test, 46 - 2, lower.tail=FALSE)


