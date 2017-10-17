setwd("~/Google Drive/projects/tcpa-rppa")

library(RColorBrewer)

data_dir = "~/DataProjects/tcpa-rppa"


# Load signatures

emt =list()
emt$epi = as.character(read.table("exploratory/figCouplingCCLE/signatures/ecad_signature_up.txt")[, 1])
emt$mes = as.character(read.table("exploratory/figCouplingCCLE/signatures/ecad_signature_down.txt")[, 1])

# emt$mes

lit_emt = list()

lit_emt$byers_2012 = list()
lit_emt$byers_2012$epi = as.character(read.table("exploratory/EMT_signatures/byers_2012/pos_ecad.txt")[, 1])
lit_emt$byers_2012$mes = as.character(read.table("exploratory/EMT_signatures/byers_2012/neg_ecad.txt")[, 1])

lit_emt$huang_2013 = list()
lit_emt$huang_2013$epi = as.character(read.table("exploratory/EMT_signatures/huang_2013/emt_epi.txt")[, 1])
lit_emt$huang_2013$mes = as.character(read.table("exploratory/EMT_signatures/huang_2013/emt_mes.txt")[, 1])

lit_emt$mak_2013 = list()
lit_emt$mak_2013$epi = as.character(read.table("exploratory/EMT_signatures/mak_2013/emt_epi.txt")[, 1])
lit_emt$mak_2013$mes = as.character(read.table("exploratory/EMT_signatures/mak_2013/emt_mes.txt")[, 1])

lit_emt$tan_2014 = list()
lit_emt$tan_2014$epi = as.character(read.table("exploratory/EMT_signatures/tan_2014/cell_line_epi.txt")[, 1])
lit_emt$tan_2014$mes = as.character(read.table("exploratory/EMT_signatures/tan_2014/cell_line_mes.txt")[, 1])

lit_emt$rokavec_2017 = list()
lit_emt$rokavec_2017$epi = as.character(read.table("exploratory/EMT_signatures/rokavec_2017/emt_epi.txt")[, 1])
lit_emt$rokavec_2017$mes = as.character(read.table("exploratory/EMT_signatures/rokavec_2017/emt_mes.txt")[, 1])



# lit_emt$byers_2012$mes = as.character(read.table("exploratory/figCouplingCCLE/signatures/ecad_signature_down.txt")[, 1])


# table(emt$epi %in% lit_emt$byers_2012$epi)


pdf("exploratory/EMT_signatures/counts_cmp.pdf", height=3)
par(mfrow=c(1, 2), mar=c(8, 4, 2, 2))
epi_counts = sapply(lit_emt, function(emt_cmp) {
	table(emt$epi %in% emt_cmp$epi)
})
epi_counts = epi_counts[2:1, ]

total = sapply(lit_emt, function(x) length(unique(x$epi)))

barplot(epi_counts,
	main="Epithelial",
	ylab="Genes",
	names.arg=paste0(colnames(epi_counts), " (", total, ")"),
	las=2
)


mes_counts = sapply(lit_emt, function(emt_cmp) {
	table(emt$mes %in% emt_cmp$mes)
})
mes_counts = mes_counts[2:1, ]

total = sapply(lit_emt, function(x) length(unique(x$mes)))

barplot(mes_counts,
	main="Mesenchymal",
	ylab="Genes",
	names.arg=paste0(colnames(epi_counts), " (", total, ")"),
	las=2
)
dev.off()

# sapply(lit_emt, function(x) length(x$mes))

