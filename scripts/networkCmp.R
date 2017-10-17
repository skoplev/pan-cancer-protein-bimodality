# Comparisions of networks, both bimodal coupling and directed acyclic graphs
# from Bayesian network inference.
# Compares TCPA cell line with patient data.

rm(list=ls())

library(bdmerge)  # for loading TCPA data

dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_51.jdk/Contents/Home/jre/lib/server/libjvm.dylib')
library(rcausal)
library(RColorBrewer)

setwd("~/Google Drive/projects/tcpa-rppa")

source("lib/parseData.r")
source("lib/bimodal.r")
source("lib/annotate.r")
source("lib/basePlots.r")


# Load TCPA data. Rows are samples, cols are proteins.
tcpa_pat = readData(file_path="data/tcpa_tcga.h5")
rownames(tcpa_pat$data) = tcpa_pat$meta_row$pr_symbol  # ensure valid rowname for data matrix
colnames(tcpa_pat$data) = tcpa_pat$meta_col$TCGA_patient_barcode
tcpa_pat = invert(tcpa_pat)
tcpa_line = parseTCPA(file_path="data/MCLP-RBN-v1.0-whole_set.tsv")


# Bimodal coupling coefficient networks
# ---------------------------------------------------------

# Fit two-component Gaussian distributions
em_pat = fitGMM(tcpa_pat$data)

# Calculate delta BIC
delta_bic_pat = deltaBIC(tcpa_pat$data, em_pat)

# calculate posterior probabilities
prob_low_pat = probLowModelMat(tcpa_pat$data, em_pat)

# Calculate bimodal coupling coefficients
bimodal_coupling_coeff_pat = cor(prob_low_pat, use="pairwise.complete.obs", method="spearman")

# Correlation network, for comparison
cor_pat = cor(tcpa_pat$data, use="pairwise.complete.obs", method="pearson")

# Same for the cell line data
em_line = fitGMM(tcpa_line$data)
delta_bic_line = deltaBIC(tcpa_line$data, em_line)
prob_low_line = probLowModelMat(tcpa_line$data, em_line)
bimodal_coupling_coeff_line = cor(prob_low_line, use="pairwise.complete.obs", method="spearman")
cor_line = cor(tcpa_line$data, use="pairwise.complete.obs", method="pearson")


# Bayesian network inference using fast greedy search algorithm
# -----------------------------------------------------------------
n_bootstrap = 200  # 100 => 4h runtime

causal_net_line = list()
causal_net_pat = list()
for (i in 1:n_bootstrap) {
	cat("Bootstrap iteration ", i, " out of ", n_bootstrap, "\n")

	# run fast greedy search algorithm and store results
	causal_net_line[[i]] = fgs(
		tcpa_line$data[
			sample(nrow(tcpa_line$data), replace=TRUE),  # sample with replacement, same number of samples (but repeated)
		]
	)

	causal_net_pat[[i]] = fgs(
		tcpa_pat$data[
			sample(nrow(tcpa_pat$data), replace=TRUE),  # sample with replacement, same number of samples (but repeated)
		]
	)
}

# No bootstrap
# causal_net_line = fgs(tcpa_line$data)  # cell line data
# causal_net_pat = fgs(tcpa_pat$data)  # patient data


# Calculates the average adjacency matrix of list of Bayesian networks.
# Used for bootstrap calculations
meanCausalAdjMat = function(causal_nets) {
	# Get all adjacency matrices
	adj_mats = lapply(causal_nets, function(net) {
		adjacency_mat = as(net$graphNEL, "matrix")
	})

	# Calculate mean matrix
	mean_adj_mat = Reduce("+", adj_mats) / length(adj_mats)

	return(mean_adj_mat)
}


# Calculate average adjacency matrix
causal_net_line_adj = meanCausalAdjMat(causal_net_line)
causal_net_pat_adj = meanCausalAdjMat(causal_net_pat)


# causal_net_pat_adj = as(causal_net_pat$graphNEL, "matrix")
# causal_net_line_adj = as(causal_net_line$graphNEL, "matrix")

# Match protein ids from cell line and patient data
line2pat_id = match(
	stemName(rownames(bimodal_coupling_coeff_pat)),
	stemName(rownames(bimodal_coupling_coeff_line))
)

network_simil = list()

# Correlation matrix comparison
network_simil[["Correlation"]] = cor.test(
	cor_pat[lower.tri(cor_pat)],
	cor_line[line2pat_id, line2pat_id][lower.tri(cor_pat)],
	use="pairwise.complete.obs"
)

# Bimodal coupling
network_simil[["Bimodal coupling"]] = cor.test(
	bimodal_coupling_coeff_pat[lower.tri(bimodal_coupling_coeff_pat)],
	bimodal_coupling_coeff_line[line2pat_id, line2pat_id][lower.tri(bimodal_coupling_coeff_pat)],
	use="pairwise.complete.obs"
)

# Causal network
network_simil[["Causal bootstrap"]] = cor.test(
	as.vector(causal_net_pat_adj),
	as.vector(causal_net_line_adj[line2pat_id, line2pat_id]),
	use="pairwise.complete.obs"
)


# # Causal network-vs bimodal
# network_simil[["Causal-bimodal cell line"]] = cor.test(
# 	as.vector(causal_net_line_adj),
# 	as.vector(bimodal_coupling_coeff_line),
# 	use="pairwise.complete.obs"
# )



# From Claudin-7
sort(causal_net_line_adj[rownames(causal_net_line_adj) == "CLAUDIN7",])
sort(causal_net_pat_adj[rownames(causal_net_pat_adj) == "Claudin-7",])

# Targeting E-cadherin
sort(causal_net_pat_adj[, colnames(causal_net_pat_adj) == "E-Cadherin"])
sort(causal_net_line_adj[, colnames(causal_net_line_adj) == "ECADHERIN"])

# Targeting Rab-25
sort(causal_net_pat_adj[, colnames(causal_net_pat_adj) == "Rab-25"])
sort(causal_net_line_adj[, colnames(causal_net_line_adj) == "RAB25"])

# From E-cadherin
sort(causal_net_pat_adj[rownames(causal_net_pat_adj) == "E-Cadherin", ])
sort(causal_net_line_adj[rownames(causal_net_line_adj) == "ECADHERIN", ])



library(reshape2)


netw_tab = melt(causal_net_line_adj)
colnames(netw_tab) = c("from", "to", "causal_weight")

bimodal_coupling_coeff_line_edges = melt(bimodal_coupling_coeff_line)

# Test that all edges are aligned
if (!all(netw_tab[, 1] == bimodal_coupling_coeff_line_edges[, 1]) | 
	!all(netw_tab[, 2] == bimodal_coupling_coeff_line_edges[, 2]))
{
		stop("edge mismatch")
}

netw_tab$bimodal_cor = bimodal_coupling_coeff_line_edges$value


netw_tab = netw_tab[order(netw_tab$causal_weight, decreasing=TRUE), ]



write.csv(netw_tab, "figures/tables/causal_netw_cell_line.csv", row.names=FALSE)


netw_tab[netw_tab$from == "CLAUDIN7" & netw_tab$to == "RAB25", ]
netw_tab[netw_tab$from == "MACC1" & netw_tab$to == "ECADHERIN", ]
netw_tab[netw_tab$from == "HSP27_pS82" & netw_tab$to == "ECADHERIN", ]
netw_tab[netw_tab$from == "LKB1" & netw_tab$to == "CLAUDIN7", ]
netw_tab[netw_tab$from == "STATHMIN" & netw_tab$to == "EIF2A_pS51", ]
netw_tab[netw_tab$from == "STATHMIN" & netw_tab$to == "ECADHERIN", ]
netw_tab[netw_tab$from == "CHK1" & netw_tab$to == "ECADHERIN", ]
netw_tab[netw_tab$from == "CHK1" & netw_tab$to == "RAB25", ]
netw_tab[netw_tab$from == "CLAUDIN7" & netw_tab$to == "MACC1", ]
# netw_tab[netw_tab$from == "ECADHERIN" & netw_tab$to == "HSP27_pS82", ]





causal_net_pat_adj_edges = melt(causal_net_pat_adj)
causal_net_pat_adj_edges = causal_net_pat_adj_edges[order(causal_net_pat_adj_edges$value, decreasing=TRUE), ]

causal_net_line_adj_edges = melt(causal_net_line_adj)
causal_net_line_adj_edges = causal_net_line_adj_edges[order(causal_net_line_adj_edges$value, decreasing=TRUE), ]


head(causal_net_line_adj_edges, 50)

# bimodal_coupling_coeff_pat_edges = melt(bimodal_coupling_coeff_pat)




pat_proteins = c("E-Cadherin", "Claudin-7", "Rab-25")
pat_proteins = c("E-Cadherin", "Claudin-7", "Rab-25")
idx = causal_net_pat_adj_edges[, 1] %in% pat_proteins & causal_net_pat_adj_edges[, 2] %in% pat_proteins
causal_net_pat_adj_edges[idx, ]


sort(colnames(causal_net_pat_adj))
sort(colnames(causal_net_line_adj))

# pat_proteins = c("ECADHERIN", "CLAUDIN7", "RAB25")
# pat_proteins = c("ECADHERIN", "CLAUDIN7", "RAB25", "MACC1")
pat_proteins = c("ECADHERIN", "CLAUDIN7", "RAB25", "BETACATENIN")
# pat_proteins = c("ECADHERIN", "CLAUDIN7", "RAB25", "MACC1", "ALPHACATENIN")
# pat_proteins = c("ECADHERIN", "CLAUDIN7", "RAB25", "MACC1", "ALPHACATENIN", "LKB1")
# pat_proteins = c("ECADHERIN", "CLAUDIN7", "RAB25", "MACC1", "ALPHACATENIN", "LKB1", "GATA3")
idx = causal_net_line_adj_edges[, 1] %in% pat_proteins & causal_net_line_adj_edges[, 2] %in% pat_proteins
causal_net_line_adj_edges[idx, ]



# Plots
# ----------------------------------------------------
pdf("exploratory/figNetworkCmp/patient_cell_line2.pdf", height=2.7)
par(mfrow=c(1, 3))

# Barplots
bar_centers = barplot(
	sapply(network_simil, function(test) test$estimate),
	names.arg=c("Correlation", "Bimodal coupling", "Causal bootstrap"),
	ylim=c(0, 0.35),
	ylab="Network similarity",
	main="Tumor samples vs cell lines",
	# col=brewer.pal(9, "Pastel2"),
	col=rgb(0.9, 0.9, 0.9),
	border=rgb(0.5, 0.5, 0.5),
	# border=NA,
	cex.names=0.6,
	las=2)
plotErrBars(
	x=bar_centers,
	confidence=lapply(network_simil, function(test) test$conf.int),
	err_width=0.1,
	lwd=1.5
	)


# Interaction rank plot from bootstrap simulation
# ----------------------------------------------------

line_cols = brewer.pal(2, "Set1")
m = 1000  # top protein-protein causal interactions to show in rank plot


# Proteins to highlight
prot1 = "Claudin-7"
prot2 = "Rab-25"

i = which(rownames(bimodal_coupling_coeff_pat) == prot1)
j = which(rownames(bimodal_coupling_coeff_pat) == prot2)

plot(sort(causal_net_pat_adj, decreasing=TRUE)[1:m],
	type="l",
	ylim=c(0, 1),
	main="Causal bootstrap",
	ylab="Weight",
	xlab="Rank",
	col=line_cols[1]
	)
lines(sort(causal_net_line_adj, decreasing=TRUE)[1:m],
	col=line_cols[2]
	)

pat_weight = causal_net_pat_adj[i, j]
line_weight = causal_net_line_adj[line2pat_id, line2pat_id][i, j]

pat_rank = sum(sort(causal_net_pat_adj) >= pat_weight)
line_rank = sum(sort(causal_net_line_adj, ) >= line_weight)


points(pat_rank, pat_weight,
	pch=16, col=line_cols[1])
text(pat_rank, pat_weight, labels=paste(prot1, "->", prot2), 
	col=line_cols[1], pos=1,
	cex=0.7
	)
points(line_rank, line_weight,
	pch=16, col=line_cols[2])
text(line_rank, line_weight, labels=paste(prot1, "->", prot2),
	col=line_cols[2], pos=1,
	cex=0.7
	)


plot(
	sort(abs(bimodal_coupling_coeff_pat[lower.tri(bimodal_coupling_coeff_pat)]), decreasing=TRUE)[1:m],
	# sort(abs(bimodal_coupling_coeff_pat[lower.tri(bimodal_coupling_coeff_pat))]), decreasing=TRUE)[1:m],
	type="l",
	ylim=c(0, 1),
	main="Bimodal coupling",
	ylab="Weight",
	xlab="Rank",
	col=line_cols[1]
	)
lines(
	sort(abs(bimodal_coupling_coeff_line[lower.tri(bimodal_coupling_coeff_line)]), decreasing=TRUE)[1:m],
	type="l",
	ylim=c(0, 1),
	col=line_cols[2]
	)

pat_coupling = bimodal_coupling_coeff_pat[i, j]
line_coupling = bimodal_coupling_coeff_line[line2pat_id, line2pat_id][i, j]

pat_rank = sum(sort(abs(bimodal_coupling_coeff_pat[lower.tri(bimodal_coupling_coeff_pat)])) >= pat_coupling)
line_rank = sum(sort(abs(bimodal_coupling_coeff_line[lower.tri(bimodal_coupling_coeff_line)])) >= line_coupling)

points(pat_rank, pat_coupling,
	pch=16, col=line_cols[1])
text(pat_rank, pat_coupling, labels=paste(prot1, ":", prot2), 
	col=line_cols[1], pos=1,
	cex=0.7
	)
points(line_rank, line_coupling,
	pch=16, col=line_cols[2])
text(line_rank, line_coupling, labels=paste(prot1, ":", prot2),
	col=line_cols[2], pos=1,
	cex=0.7
	)

legend("topright", legend=c("Tumor samples", "Cell lines"), col=line_cols, pch=15)
dev.off()


prot1 = "E-Cadherin"
prot2 = "Claudin-7"
# prot2 = "Rab-25"


prot1 = "Claudin-7"
prot2 = "Rab-25"

i = which(rownames(bimodal_coupling_coeff_pat) == prot1)
j = which(rownames(bimodal_coupling_coeff_pat) == prot2)
# j = 63

# rownames(bimodal_coupling_coeff_pat)
bimodal_coupling_coeff_pat[i, j]
bimodal_coupling_coeff_line[line2pat_id, line2pat_id][i, j]


x = melt(bimodal_coupling_coeff_pat)
x = x[order(x$value, decreasing=TRUE), ]

x = x[x[, 1] != x[, 2], ]
head(x, 100)


cor_pat[i, j]
cor_line[line2pat_id, line2pat_id][i, j]


plot(tcpa_pat$data[,i], tcpa_pat$data[,j], xlab=prot1, ylab=prot2)
cor(tcpa_pat$data[,i], tcpa_pat$data[,j])


plot(delta_bic_pat, delta_bic_line[line2pat_id])
plot(signPseudoLog(delta_bic_pat), signPseudoLog(delta_bic_line[line2pat_id]))
abline(coef=c(0, 1), col="grey")
