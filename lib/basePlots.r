# Helper functions for barious plots

# Plot error bars
errBar = function(x0, y0, x1, y1, width, ...) {
	segments(x0, y0, x1, y1, ...)
	segments(x0 - width, y0, x1 + width, y0, ...)
	segments(x0 - width, y1, x1 + width, y1, ...)
}

# Plots error bars
# x is a vector of the mean values.
# confidence is a list of associated confidence intervals.
# tail_width is the width of the endpoint segments.
plotErrBars = function(x, confidence, 
	tail_width=0.2,  
	horiz=FALSE,
	...)  # passes inputs to segment, such as col
{
	# Input test
	if (length(x) != length(confidence)) {
		stop("Input length mismatch.")
	}

	if (length(confidence[[1]]) != 2) {
		stop("Not a list of confidence intervals")
	}

	# Error bars
	if (horiz) {
		for (i in 1:length(x)) {
			errBar(
				x[i],
				confidence[[i]][1],
				x[i],
				confidence[[i]][2],
				tail_width, ...)
		}
	} else {
		for (i in 1:length(x)) {
			errBar(
				x[i],
				confidence[[i]][1],
				x[i],
				confidence[[i]][2],
				tail_width, ...)
		}
	}
}
