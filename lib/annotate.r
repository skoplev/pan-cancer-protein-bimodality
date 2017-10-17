# Common helper functions for annotating cell lines

# Stem cell line names for standardized string comparison
stemName = function(name) {
	require(stringr)

	stripped_name = str_replace_all(name, "^NCI-", "")  # NCI cell line prefix
	stripped_name = str_replace_all(stripped_name, "^NCI", "")

	# Remove all non-alhph-numercial characters and convert to upper case
	alpha_num = str_replace_all(stripped_name, "[^[:alnum:]]", "")
	return(toupper(alpha_num))
}


# Returns url query for searching UniProt database
uniprotQueryUrl = function(query) {
	url = paste0(
		"http://www.uniprot.org/uniprot/?query=",
		query,
		"+AND+",
		"organism:9606",  # homo sapiens
		"&sort=score&columns=id,entry%20name,reviewed,protein%20names,genes,organism,length&format=tab"
	)
	return(url)
	# return(URLencode(url))  # return sanitized URL
}

# Uniprot url for getting data of an entry
uniprotUrl = function(id) {
	url = paste0(
		"http://www.uniprot.org/uniprot/",
		id,
		".xml"
	)
	return(url)
}

# Name list: from->to
# Manual override of wrong names due to bad Uniport search results 
std_rename = list(
	X1433EPSILON="YWHAE",
	X4EBP1="EIF4EBP1",
	X4EBP1_pS65="EIF4EBP1",
	X4EBP1_pT37T46="EIF4EBP1",
	X53BP1="TP53BP1",
	AMPKALPHA="PRKAA2",
	AMPK="PRKAA2",
	BCLXL="BCL2L1",
	CMYC="MYC",
	CASPASE7CLEAVEDD198="CASP7",
	COLLAGENVI="COL6A1",  # first in complex
	CYCLINB1="CCNB1",
	CYCLIND1="CCND1",
	CYCLINE1="CCNE1",
	FIBRONECTIN="FN1",
	GAB2="GAB2",
	GSK3ALPHABETA="GSK3A",
	GSK3ALPHABETA_pS21S9="GSK3A",
	NCADHERIN="CDH2",
	NFKBP65_pS536="RELA",
	PCADHERIN="CDH3",
	P38MAPK="MAPK14",
	PI3KP110ALPHA="PIK3CA",
	PR="PGR",
	STAT5ALPHA="STAT5A",
	X4EBP1_pT70="EIF4EBP1",
	BAP1C4="BAP1",
	CYCLINE2="CCNE2",
	MYOSINIIA_pS1943="MYH9",
	P21="CDKN1A",
	PI3KP85="PIK3R1",
	PKCPANBETAII_pS660="PRKCB",
	RICTOR_pT1135="RICTOR",
	RICTOR="RICTOR",
	SCD="SCD",
	TAZ="TAZ",
	ACETYLATUBULINLYS40="TUBA1A",
	P62LCKLIGAND="SQSTM1",
	X1433BETA="YWHAB",
	X1433ZETA="YWHAZ",
	DIRAS3="DIRAS3",
	PARPCLEAVED="PARP1",
	SETD2="SETD2",
	SNAIL="SNAI1",
	MYOSINIIA="MYH9",
	GYS_pS641="GYS1",
	GYS="GYS1",
	LCN2A="LCN2",
	NRF2="NFE2L2",
	PARPAB3="PARP1",
	THYMIDILATESYNTHASE="TYMS",
	TTF1="TTF1",
	CHROMOGRANINANTERM="CHGA",
	ALPHACATENIN="CTNNA1",
	ATR="ATR",
	BACTIN="ACTB",
	BCATENIN_pT41S45="CTNNB1",
	CASPASE9CLEAVED="CASP9",
	COXIV="COX4I1",
	CYCLOPHILINF="PPIF",
	DETYROSINATEDALPHATUBULIN="TUBA1A",
	DMHISTONEH3="H3F3A",
	DMK9HISTONEH3="H3F3A",
	ER="ESR1",
	FAK="PTK2",
	FAK_pY397="PTK2",
	GPBB="PYGM",
	HEXOKINASEII="HK2",
	HIAP="BIRC3",
	PI3KP110B="PI3KCB",
	PLCGAMMA1_pY783="PLCG1",
	PLCGAMMA1="PLCG1",
	PLCGAMMA2_pY759="PLCG2",
	PUMA="BBC3",
	RSK="RPS6KA1",
	TWIST="TWIST1",
	ATR_pS428="ATR",
	FLT3FLK2="FLT3",
	CASPASE3CLEAVED="CASP3",
	COMT="COMT",
	CYCLINA="CCNA1",
	CYTOKERATIN17="KRT17",
	ER_pS167="ESR1",
	P90RSK_pS221="RPS6KA1",
	TAZ_pS79="TAZ",
	BLNK_pY84="BLNK",
	BLNK_pY96="BLNK",
	BLNK="BLNK",
	CYCLINE4="CDK4",
	HISTONEH2AX_pS139="H2AFX",
	LAT_pY191="LAT",
	LAT="LAT",
	PLCGAMMA2_pY1217="PLCG1",
	PLCGAMMA2="PLCG1",
	PYK2_pY402="PTK2B",
	SPIB="SPIB",
	DNATOPOISOMERASE2ALPHA="TOP2A",
	DNATOPOISOMERASE2BETA="TOP2B",
	CASEINKINASE="CSNK1G1",
	ERG123="ERG",
	FAK_pY925="PTK2",
	PKCDELTA_pT507="PRKCD",
	PKCDELTA_pS664="PRKCD",
	TIF1ALPHA="TRIM24",
	GAB2="GAB2",
	GAB2_pY452="GAB2",
	PAR="F2R",
	MITOCHONDRIA=NA,
	HER2="ERBB2",
	HER2_pY1248="ERBB2",
	BAX="BAX"
)

# Construct annotation table for protein ids on the format SYMBOL_PTM.
# Takes as input a vector of protein ids.
# Queries uniprot for recommended names. This can fail which is why there
# is a lookup table for proteins.
makeProteinTable = function(prot_ids, rename=std_rename) {
	require(RCurl)

	# Annotate proteins
	proteins = data.frame(mclp_id=prot_ids)

	# Name of protein
	proteins$mclp_name = sapply(
		strsplit(as.character(proteins$mclp_id), "_"),
		function(entry) {
			return(entry[1])  # return first element of split string
	}) 

	# Post translational modification
	proteins$mclp_ptm = sapply(
		strsplit(as.character(proteins$mclp_id), "_"),
		function(entry) {
			return(entry[2])  # return second element of split string
	}) 

	# Get urls for Uniprot searches
	urls = uniprotQueryUrl(proteins$mclp_name)

	# Warning: some mismatches occur with the Uniport search.
	# manual list of gene names which overwrites catched mistakes.
	search_results_tsv = getURL(as.list(urls), async=FALSE)  # async=FALSE takes longer but is safer

	# Get the Uniprot IDs from the search tables
	proteins$uniprot_id = lapply(search_results_tsv, function(tsv) {
		try({
			search_table = read.csv(textConnection(tsv), sep="\t")
			return(as.character(search_table$Entry[1]))
		}, {
			warning("missing data")
			return(NA)
		})
	})

	proteins$uniprot_name = lapply(search_results_tsv, function(tsv) {
		try({
			search_table = read.csv(textConnection(tsv), sep="\t")
			gene_name = strsplit(as.character(search_table$Gene.names[1]), " ")[[1]][1]  # recommended name of top result
			return(gene_name)
		}, {
			warning("missing data")
			return(NA)
		})
	})


	# rename uniprot_name
	for (i in 1:nrow(proteins)) {
		if (proteins$mclp_id[i] %in% names(rename)) {
			proteins$uniprot_name[i] = rename[[as.character(proteins$mclp_id[i])]]
		}
	}

	return(proteins)
}