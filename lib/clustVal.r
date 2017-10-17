# Implementation of the Gap statistic from (Tibsharani, 2001) estimating the supported number of clusters in the dataset.
# ----------------------------------------------

# Calculates the within cluster distance given a distance matrix and a set of class labels.
withinClustSquares = function(dist_mat, groups) {
	if (!class(dist_mat) == "dist") {
		stop("Input not a distance object.")
	}

	if (attr(dist_mat, "Size") != length(groups)) {
		stop("Dimensions of distance matrix and the group labels disagree.")
	}

	within_sum_of_squares_around_mean = 0.0

	for (i in 1:length(unique(groups))) {
		# Data columns associated
		samples = which(groups == unique(groups)[i])

		within_sum_of_squares_around_mean = within_sum_of_squares_around_mean + 
			sum(as.matrix(dist_mat)[samples, samples]) / (2 * length(samples))
	}

	return(within_sum_of_squares_around_mean)
}

# Generate n reference datasets for the provided data_mat datasets.
# The reference data are drawn from the null distribution that have the the ranges in the primary.
# Returns list of reference data matrices.
# mat: data matrix, features x samples
genUniformRef = function(mat, n) {
	# feature ranges for generating uniformly distributed reference data
	ranges = apply(mat, 2, range, na.rm=TRUE) 

	ref_mats = list()  # data structure storing the reference data matrices
	for (i in 1:n) {
		# Init reference data matrix
		ref_mat = matrix(NA, nrow=nrow(mat), ncol=ncol(mat))
		colnames(mat) = colnames(mat)
		rownames(mat) = rownames(mat)

		# Loop over each feature
		for (p in 1:ncol(mat)) {
			ref_mat[,p] = runif(
				nrow(mat),  # observations
				min=ranges[1, p],
				max=ranges[2, p]
				)
		}
		ref_mats[[i]] = ref_mat
	}
	return(ref_mats)
}

# Calculate Gap statistic for a data set using hierarcical clustering.
# mat: data matrix, samples x features
# max_clusters: maximum number of clusters tested.
# n_null: the number of reference (null) sample to draw.
gapStatHc = function(mat, max_clusters, n_null) {
	dmat = dist(mat)  # calculate distance matrix

	# Generate null data sets.
	ref_mats = genUniformRef(mat, n_null)

	# Calcualte reference distances
	ref_dmats = lapply(ref_mats, dist)

	# Hierarcical clustering of primary data
	hc = hclust(dmat, method="average")  # TCPA

	# Cluster data by cuttting the dendrogram
	groups = list()
	for (k in 1:max_clusters) {
		groups[[k]] = cutree(hc, k=k)
	}

	# Calculate groups of reference data
	ref_groups = list()
	for (b in 1:n_null) {
		# Cluster the reference data
		hc_ref = hclust(ref_dmats[[b]], method="average")

		ref_groups[[b]] = list()

		# cut tree at different 
		for (k in 1:max_clusters) {
			ref_groups[[b]][[k]] = cutree(hc_ref, k=k)
		}
	}

	gap = calcGapStat(
		dmat, groups, 
		ref_dmats, ref_groups)

	return(gap)
}

# Calculate the gap statistics using k-means clustering
gapStatKmeans = function(mat, max_clusters, n_null) {
	dmat = dist(mat)


	ref_mats = genUniformRef(mat, n_null)


	ref_dmats = lapply(ref_mats, dist)

	mat[is.na(mat)] = 0.0


	# Cluster data by different ks
	groups = list()
	for (k in 1:max_clusters) {
		kmeans_results = kmeans(mat, centers=k)
		groups[[k]] = kmeans_results$cluster
	}

	ref_groups = list()
	for (b in 1:n_null) {
		ref_groups[[b]] = list()

		for (k in 1:max_clusters) {
			kmeans_result = kmeans(ref_mats[[b]], centers=k)
			ref_groups[[b]][[k]] = kmeans_result$cluster
		}
	}

	gap = calcGapStat(
		dmat, groups, 
		ref_dmats, ref_groups)

	return(gap)
}

# Evaluate the gap statistic
# dmat: primary distance matrix
# groups: group classes for different number of clusters of the objects in distance matrix. List of vectors.
# ref_dmat: list of reference distance matrices generated under some null model.
# ref_groups: list of groups
calcGapStat = function(dmat, groups, ref_dmats, ref_groups) {

	# Calculate within cluster sum of squares for primary data
	within = sapply(groups, function(group) {
		return(withinClustSquares(dmat, group))
	})

	# Cluster each reference dataset and calculate the within-cluster sum of squares for different
	within_ref = matrix(NA, nrow=length(ref_dmats), ncol=length(ref_groups[[1]]))
	for (b in 1:length(ref_dmats)) {
		# hc_ref = hclust(ref_dmats[[b]], method="average")

		# Iterate over cluster numbers
		for (k in 1:max_clusters) {
			# Calculate grouping
			# groups = cutree(hc_ref, k=k)

			# Calculate within-cluster sum of squares
			within_ref[b, k] = withinClustSquares(ref_dmats[[b]], ref_groups[[b]][[k]])
		}
	}

	# Calculate the Gap statistic
	ave_log_within_ref = (1/n_null) * apply(log(within_ref, 2), 2, sum)

	# standard deviation
	sd_log_within_ref = apply(log(within_ref, 2), 2, sd)

	# gap statistic
	gap = ave_log_within_ref - log(within, 2)

	return(list(
			gap=gap,
			sd=sd_log_within_ref,
			within=within
		)
	)
}


# Plot gap statistic output
gapPlot = function(gap, col, ...) {
	seg_width = 0.03

	n = length(gap$gap)

	plot(gap$within, type="l",
		col=col,
		xlab="Clusters",
		ylab="Within-cluster squares",
		...)

	# viewport
	plot(c(1:length(gap$gap), 1:length(gap$gap)), c(gap$gap - gap$sd, gap$gap + gap$sd),
		xlab="Clusters",
		ylab="Gap statistic",
		type="n",
		...)

	# Mean Gap statistic line
	lines(gap$gap, col=col)

	# Sd error bars
	segments(c(1:n), gap$gap - gap$sd, c(1:n), gap$gap + gap$sd, col=col)
	# top and bottom error bar
	segments(c(1:n) - seg_width, gap$gap + gap$sd, c(1:n) + seg_width, gap$gap + gap$sd, col=col)
	segments(c(1:n) - seg_width, gap$gap - gap$sd, c(1:n) + seg_width, gap$gap - gap$sd, col=col)
}
