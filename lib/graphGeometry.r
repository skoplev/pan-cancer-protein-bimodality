# Rescale numbers to [-1, 1] interval based on min and max value.
rescale = function(vec) {
	return(
		2 * (vec - min(vec, na.rm=TRUE)) / 
		(max(vec, na.rm=TRUE) - min(vec, na.rm=TRUE)) - 1
	)
}

# Rotate a xy column matrix around (0, 0)
rotate = function(xy, theta) {
    xy = as.matrix(xy)
    cos.angle = cos(theta)
    sin.angle = sin(theta)
    xy.rot <- xy %*% t(matrix(c(cos.angle, sin.angle, -sin.angle, cos.angle), 2, 2))
    return(xy.rot)
}

# Calculates angle in radians between two vectors
angle = function(x, y) {
	theta = acos(sum(x * y) / (sqrt(sum(x * x)) * sqrt(sum(y * y))))
	return(theta)
}

translate = function(xy, point) {
	xy = as.matrix(xy)
	xy_translate = xy
	xy_translate[,1] = xy[,1] + point[1]
	xy_translate[,2] = xy[,2] + point[2]
	return(xy_translate)


}

# Returns list of rectangle coordinates placed next to each other and 
# mode=c("stacked", "adjacent")
createRectStack = function(vec, width, height_scale, mode="stacked") {
	pen_pos = 0.0  # current rendering position

	poly = list()
	for (i in 1:length(vec)) {
		height = vec[i] * height_scale  # calculate height in stack
		rect_xy = createRect(width, height)  # new rectangle data

		if (mode == "stacked") {
			rect_xy = translate(rect_xy, c(0.0, pen_pos))
			pen_pos = pen_pos + height  # increase pen height

		} else if (mode == "adjacent") {
			rect_xy = translate(rect_xy, c(pen_pos, 0.0))
			pen_pos = pen_pos + width
		} else {
			stop("Invalid mode")
		}

		poly[[i]] = rect_xy

	}
	return(poly)
}

# Plots a rectangle stack at a given angle
polygonRectStack = function(rects, cols=NULL, at=c(0, 0), theta=0) {
	if (length(at) != 2) {
		stop("Argument at is an invalid 2D coordinate.")
	}

	if (is.null(cols)) {
		cols = rep("black", length(rects))
	}

	for (i in 1:length(rects)) {
		rect_xy = rects[[i]]
		rect_xy = rotate(rect_xy, theta)
		rect_xy = translate(rect_xy, at)
		polygon(rect_xy[,1], rect_xy[,2],
			col=cols[i],
			border=NA
		)
	}
}

# Creates polygons points for a rectangle
createRect = function(width, height) {
	# Rectangle buffer
	rect_xy = matrix(c(
		0.0, 0.0,
		0.0, height,
		width, height,
		width, 0.0
		), ncol=2, byrow=TRUE
	)
	return(rect_xy)
}
