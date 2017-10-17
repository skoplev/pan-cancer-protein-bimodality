# A variety of color functionality

# Colors of maximum separability adapted from the colour alphabet.
# https://eleanormaclure.files.wordpress.com/2011/03/colour-coding.pdf
colors = c(rgb(240,163,255,maxColorValue=255),rgb(0,117,220,maxColorValue=255),rgb(153,63,0,maxColorValue=255),rgb(76,0,92,maxColorValue=255),rgb(25,25,25,maxColorValue=255),rgb(0,92,49,maxColorValue=255),rgb(43,206,72,maxColorValue=255),rgb(255,204,153,maxColorValue=255),rgb(128,128,128,maxColorValue=255),rgb(148,255,181,maxColorValue=255),rgb(143,124,0,maxColorValue=255),rgb(157,204,0,maxColorValue=255),rgb(194,0,136,maxColorValue=255),rgb(0,51,128,maxColorValue=255),rgb(255,164,5,maxColorValue=255),rgb(255,168,187,maxColorValue=255),rgb(66,102,0,maxColorValue=255),rgb(255,0,16,maxColorValue=255),rgb(94,241,242,maxColorValue=255),rgb(0,153,143,maxColorValue=255),rgb(224,255,102,maxColorValue=255),rgb(116,10,255,maxColorValue=255),rgb(153,0,0,maxColorValue=255),rgb(255,255,128,maxColorValue=255),rgb(255,255,0,maxColorValue=255),rgb(255,80,5,maxColorValue=255))

# Linear color interpolation factory.
colorGradient = function(x, colors=c("red","yellow","green"), colsteps=100) {
	pal = colorRampPalette(colors) (colsteps)  # the color palette
	return(
		pal[findInterval(x, seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out=colsteps))]
 	)
}

# Plots a gradient scale
# pal: color palette
plotColorBar = function(pal, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(pal)-1)/(max-min)

    # dev.new(width=1.75, height=5)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(2, ticks, las=1)
    for (i in 1:(length(pal)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=pal[i], border=NA)
    }
}
