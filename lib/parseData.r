# Functions that parse specific types of data

# Load TCPA data from table
parseTCPA = function(file_path) {
	# Load data as data frame
	d = read.table(file_path, header=TRUE, sep="\t")

	# Format data matrix, cell lines x protein measurements
	x = d[,4:ncol(d)]
	x = as.matrix(x)
	rownames(x) = d[,1]  #

	return(list(
		data=x,
		meta_col=data.frame(pr_symbol=colnames(x)),
		meta_row=d[,1:3]
	))
}