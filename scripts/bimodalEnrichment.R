
library(data.table)

setwd("~/Google Drive/projects/tcpa-rppa")

go_bimodal = fread("data/bimodal/enrichment/GO_BP_bimodal_common.txt")
go_unimodal = fread("data/bimodal/enrichment/GO_BP_nonbimodal.txt")


go_unimodal_sig = go_unimodal[go_unimodal[["upload_1 (P-value)"]] < 0.05, ]
go_bimodal_sig = go_bimodal[go_bimodal[["upload_1 (P-value)"]] < 0.05, ]


idx = !go_bimodal_sig[["GO biological process complete"]] %in% go_unimodal_sig[["GO biological process complete"]]

go_bimodal_sig[idx, ]

head(go_bimodal_sig[idx, ], 30)

pvals = unlist(go_bimodal_sig[idx, "upload_1 (P-value)"])
names(pvals) = unlist(go_bimodal_sig[idx, "GO biological process complete"])

# Remove GO IDs from names
names(pvals) = sapply(strsplit(names(pvals), " \\("), function(x) x[1])


pdf("exploratory/bimodal/enrichment/bimodal_common_not_nonbimodal.pdf", width=6, height=5)
par(mar=c(16, 4, 2, 2))
barplot(-log10(pvals[1:30]), las=2, ylab=expression("-log"[10] * " p (BH)"),
	main="",
	cex.names=0.6)
dev.off()

barplot(-log10(pvals), las=2, ylab=expression("-log"[10] * " p"),
	cex.names=0.6)


