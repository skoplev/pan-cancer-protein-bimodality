# L1000CDS2 prediction visualizations of EMT and MET signatures.
# Based on mRNA bimodal coupling to E-Cadherin protein.

rm(list=ls())

# Collect scores of same small molecules based ranking from L1000CDS2 analysis
# Input: L1000CDS2 output table.
collectL1000Scores = function(l1000cds2) {
	l1000_scores = list()
	for (drug in unique(l1000cds2$Perturbation)) {
		l1000_scores[[drug]] = l1000cds2$score[l1000cds2$Perturbation == drug]
	}
	return(l1000_scores)
}

setwd("/Users/sk/Desktop/tcpa-rppa")

# Reverse signature indicicating EMT.
# ----------------------------------------------------
l1000_rev = fread("data/couplingCCLE/signatures/L1000/L1000CDS2_reverse.tsv")
l1000_rev = as.data.frame(l1000_rev)

l1000_rev_scores = collectL1000Scores(l1000_rev)

pdf("exploratory/figCouplingCCLE/emt_signature/L1000/ecad_rev_EMT.pdf", width=4, height=4)
par(
	# bty="n",
	mar=c(11, 5, 4, 2), cex.axis=0.7)
boxplot(
	l1000_rev_scores,
	las=2,
	ylab="EMT score [0, 1]",
	col=rgb(0.9, 0.9, 0.9),
	xaxt="n"
)
axis(side=1, at=1:length(l1000_rev_scores), names(l1000_rev_scores), tick=FALSE, las=2)
dev.off()

# Mimic signature indicating MET
# --------------------------------------------------------
l1000_mimic = fread("data/couplingCCLE/signatures/L1000/L1000CDS2_mimic.tsv")
l1000_mimic = as.data.frame(l1000_mimic)

# Collect scores of same small molecules based ranking from L1000CDS2 analysis
l1000_mimic_scores = collectL1000Scores(l1000_mimic)


pdf("exploratory/figCouplingCCLE/emt_signature/L1000/ecad_mimic_MET.pdf", width=4, height=4)
par(
	# bty="n",
	mar=c(11, 5, 4, 2), cex.axis=0.7)
boxplot(
	l1000_mimic_scores,
	las=2,
	ylab="MET score [0, 1]",
	col=rgb(0.9, 0.9, 0.9),
	xaxt="n"
)
axis(side=1, at=1:length(l1000_mimic_scores), names(l1000_mimic_scores), tick=FALSE, las=2)
dev.off()
