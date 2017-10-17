require(mixtools)

# Plot histogram of data matrix x using vector of feature ids.
plotHist = function(x, ids, breaks=NULL, main=NULL, col="black", xlim=NULL, ...) {

	# Default name
	if (is.null(main)) {
		mains = colnames(x)
	} else {
		mains = rep(main, ncol(x))
	}

	if (is.null(xlim)) {
		xlim = range(x, na.rm=TRUE)
	}

	if (is.null(breaks)) {
		breaks = seq(range(x, na.rm=TRUE)[1], range(x, na.rm=TRUE)[2], length.out=100)
	}

	for (i in ids) {
		hist(x[,i],
			# breaks=seq(range(x, na.rm=TRUE)[1], range(x, na.rm=TRUE)[2], length.out=100),
			breaks=breaks,
			border=NA,
			main=mains[i],
			ylab=NULL,
			xlab=NULL,
			cex.main=0.9,
			col=col,
			xlim=xlim,
			...)
	}
}

# Calculate the probability that an observation at x belongs to the lowest mixture of a
# two-compnent Gaussian mixture model.
# lambda: vector of the two mixture coefficients
# mu: vector of the two means
# sigma: vector of the two standards deviations
probLowModel = function(x, lambda, mu, sigma) {
	if (length(lambda) != 2) {
		stop("Lambda not a vector of 2.")
	}

	if (length(mu) != 2) {
		stop("mu is not a vector of 2.")
	}

	if (length(sigma) != 2) {
		stop("sigma is not a vector of 2.")
	}

	if (lambda[1] > 1.0 | lambda[1] < 0.0 | lambda[2] > 1.0 | lambda[2] < 0.0) {
		stop("Invalid lambda value.")
	}

	m1 = which.min(mu)  # index of the minimum (lowest) component
	m2 = which.max(mu)

	# calculate the posterior probability model 1 given x
	# Bayes'
	p = lambda[m1] * dnorm(x, mu[m1], sigma[m1]) / 
		(lambda[m1] * dnorm(x, mu[m1], sigma[m1]) + lambda[m2] * dnorm(x, mu[m2], sigma[m2]))

	return(p)
}

# em is an output object from a Gaussian EM fit.
probLowModelMat = function(x, em) {
	if (length(em$results) != ncol(x)) {
		stop("Invalid matrix and em dimension")
	}
	
	prob_low = matrix(NA, nrow=nrow(x), ncol=ncol(x))
	colnames(prob_low) = colnames(x)
	rownames(prob_low) = rownames(x)
	for (i in 1:nrow(x)) {
		for (j in 1:ncol(x)) {
			try({
				prob_low[i, j] = probLowModel(
					x=x[i, j],
					lambda=em$results[[j]]$lambda,
					mu=em$results[[j]]$mu, 
					sigma=em$results[[j]]$sigma
				)
			})
		}
	}

	return(prob_low)
}

# Calculate the loglikelihood from 
# log-likelihoods of single normal distribution fit
normLoglik = function(x) {
	row_mean = apply(x, 2, mean, na.rm=TRUE)
	row_sd = apply(x, 2, sd, na.rm=TRUE)
	norm_loglik = rep(NA, ncol(x))
	for (i in 1:ncol(x)) {
		# norm_loglik[i] = sum(log(pnorm(x[,i], prot_means[i], prot_sd[i])), na.rm=TRUE)
		norm_loglik[i] = sum(log(dnorm(x[,i], row_mean[i], row_sd[i])), na.rm=TRUE)
	}
	return(norm_loglik)
}



# Fit two-component Gaussian mixture models to each protein observation (per column)
fitGMM = function(x) {
	em_results = list()
	for (i in 1:ncol(x)) {
		col = x[,i]

		values = col[!is.na(col)]  # remove missing values

		try({
			em_result = normalmixEM(values)
			em_results[[i]] = em_result  # store result
		})
	}

	# Extract log-likelihood values of two-mixture Gaussian fit
	em_loglik = sapply(em_results, function(res) {
		if (!is.null(res$loglik)) {
			return(res$loglik)
		} else {
			return(NA)
		}
	})
	names(em_loglik) = colnames(x)
	em_loglik_sort = sort(em_loglik, decreasing=TRUE, index.return=TRUE)
	em_loglik_order = match(names(em_loglik_sort$x), colnames(x))  # indicies refering to the column ids


	# Extract lambda (prior) for low mixture component
	em_lambda = sapply(em_results, function(em) {
		if (!is.null(em$lambda)) {
			low = 1
			min_ind = which.min(em$mu)
			return(em$lambda[min_ind])
		} else {
			return(NA)
		}
	})

	return(
		list(
			results=em_results,
			lambda=em_lambda,
			loglik=em_loglik
		)
	)
}

# Calculate contingency tables for each classifier based on threshold
calcTissueEntropy = function(prob_low, tissue, alpha) {
	if (nrow(prob_low) != length(tissue)) {
		stop("Row tissue dimension mismatch")
	}

	class_stat = data.frame(row.names=colnames(prob_low))

	# Loop over each feature in matrix of posterior probabilities
	for (i in 1:ncol(prob_low)) {
		# Contingency table for low bin. Number of various tissues.
		low_conting = table(tissue[prob_low[,i] > 1 - alpha])

		# Contingency table for high bin
		high_conting = table(tissue[prob_low[,i] < alpha])

		# total number of cell lines considered, the cell lines separated.
		m = sum(low_conting) + sum(high_conting) 


		# Calculate the class probabilities. Classical statistics
		low_class_prob = low_conting / sum(low_conting)
		high_class_prob = high_conting / sum(high_conting)

		low_entropy = -sum(low_class_prob * log(low_class_prob, 2), na.rm=TRUE)  # removes log(0) infinites
		high_entropy = -sum(high_class_prob * log(high_class_prob, 2), na.rm=TRUE)  # removes log(0) infinites

		total_entropy = weighted.mean(c(low_entropy, high_entropy), c(sum(low_conting) / m, sum(high_conting) / m), na.rm=TRUE)
		min_entropy = min(low_entropy, high_entropy)

		low_purity = max(low_class_prob)
		high_purity = max(high_class_prob)

		total_purity = weighted.mean(c(low_purity, high_purity), c(sum(low_conting) / m, sum(high_conting) / m), na.rm=TRUE)


		class_stat$low_n[i] = sum(low_conting)
		class_stat$high_n[i] = sum(high_conting)
		class_stat$low_entropy[i] = low_entropy
		class_stat$high_entropy[i] = high_entropy
		class_stat$total_entropy[i] = total_entropy
		class_stat$min_entropy[i] = min_entropy

		class_stat$low_purity[i] =  low_purity
		class_stat$high_purity[i] = high_purity
		class_stat$total_purity[i] = total_purity
	}

	return(class_stat)
}

# Bayesian information criterion test.
# Compares
deltaBIC = function(x, em) {

	norm_loglik = normLoglik(x)  # log likelihood of single Gaussian model

	# Number of observations
	nobs = apply(x, 2, function(row) sum(!is.na(row)))


	em_bic = - 2 * (em$loglik) + 5 * log(nobs)
	norm_bic = - 2 * (norm_loglik) + 2 * log(nobs)

	# delta_bic = em_bic - norm_bic
	delta_bic = norm_bic - em_bic

	return(delta_bic)
}


signPseudoLog = function(x) {
	missing = is.na(x)
	x[missing] = 0  # for convenient vector calculations

	log_trans = x  # init vector

	log_trans[!is.na(x) & x > 0] = log10(x[x > 0] + 1)
	log_trans[!is.na(x) & x < 0] = -log10(-x[x < 0] + 1)

	log_trans[missing] = NA  # reintroduce missing
	return(log_trans)
}

# Plots a color legend vertically on the right side of plotting area.
legend.col = function(col, lev) {
	opar <- par
	 
	n <- length(col)
	 
	bx <- par("usr")
	 
	box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
	bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
	box.cy <- c(bx[3], bx[3])
	box.sy <- (bx[4] - bx[3]) / n
	 
	xx <- rep(box.cx, each = 2)
	 
	par(xpd = TRUE)
	for(i in 1:n) {
		yy <- c(box.cy[1] + (box.sy * (i - 1)),
		box.cy[1] + (box.sy * (i)),
		box.cy[1] + (box.sy * (i)),
		box.cy[1] + (box.sy * (i - 1)))
		polygon(xx, yy, col = col[i], border = col[i])
	}

	par(new=TRUE)
	plot(0, 0, type = "n",
		ylim = c(min(lev), max(lev)),
		yaxt = "n", ylab = "",
		xaxt = "n", xlab = "",
		frame.plot = FALSE
	)
	axis(side=4, las=2, tick=FALSE, line=.25)
	par <- opar
}

plotDeltaBIC = function(delta_bic1, delta_bic2, cross_coupling, xlab="", ylab="") {
	if (length(delta_bic1) != length(delta_bic2)) {
		stop("Invalid delta input lenghts")
	}

	if (length(delta_bic1) != length(cross_coupling)) {
		stop("Invalid cross coupling lenght")
	}

	# Generate point colors
	colorGrad = colorRamp(rev(brewer.pal(9, "Spectral")))

	color_coupling = colorGrad((cross_coupling + 1) / 2)  # [-1, 1] mapping
	color_coupling[is.na(color_coupling)] = 150  # Default grayscale, 255 scale

	# Convert rgb to colors
	pts_color = apply(color_coupling, 1, function(row) {
		rgb(row[1], row[2], row[3], maxColorValue=255)
	})

	par(mar=c(4, 4, 4, 5))

	x = signPseudoLog(delta_bic1)
	y = signPseudoLog(delta_bic2)

	# Set up viewport
	plot(x, y,
		type="n",
		xlab=xlab, ylab=ylab
	)

	abline(h=1, col="grey", lty=3)
	abline(v=1, col="grey", lty=3)
	abline(coef=c(0, 1), col="grey")

	# Fill
	points(x, y,
		# bg=pts_color,  # fill color
		col=pts_color,
		lwd=0.3,
		# col=rgb(0.1, 0.1, 0.1),  # outline color
		pch=16
	)

	# Labels
	text(x, y + 0.06,
		label=names(delta_bic2),
		cex=0.3,
		# col=rgb(0.6, 0.6, 0.6)
		col=rgb(0, 0, 0, 0.5)
	)

	# Vertical color legend
	legend.col(col=rev(brewer.pal(9, "Spectral")), lev=c(-1.0, 1.0))
}

# Get unique gene associations of top-n terms in Enrichr table
getUniqueGenes = function(enrichment, n) {
	genes = unique(
		str_split(
			paste(enrichment$Genes[1:n], collapse=";")  # Combined top-n gene associations
			, ";"
		)[[1]]
	)
	return(genes)
}

