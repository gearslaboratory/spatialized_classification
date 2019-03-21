#' import circular

plot.spatialized_classification <- function(spatialized_classification_output,
		radius_interval=1, # Units of pixels
		azimuth_interval=45, # Units of degrees
		
		# Layer plot parameters
		layer.main="Layers Used",
		layer.col="lightblue",
		layer.xlab="layer",
		
		# Radius plot parameters
		radius.main="Radii Used",
		radius.col="lightblue",
		
		# Azimuth plot parameters (see ?rose.diag)
		azimuth.main="Azimuth Used",
		azimuth.col="lightblue",
		azimuth.prop=1,
		azimuth.ticks=TRUE,
		azimuth.shrink=1		
)
{
	# library(caret)
	library(circular)
	
	if(is.null(spatialized_classification_output))
	{
		exit("Please set a spatialized_classification_output.")
	}
	
	variable_lookup  <- spatialized_classification_output$variable_lookup
	
	if("rfsrc" %in% class(spatialized_classification_output$model))
	{
		variables_used <- spatialized_classification_output$model$xvar.names
		variable_lookup_used <- spatialized_classification_output$variable_lookup[
				variable_lookup$variable %in% variables_used,,drop=FALSE]
		
		### Band Plot
		layer_breakpoints <- seq(1,(max(variable_lookup$layer)+1))
		used_layer_hist <- hist(variable_lookup_used$layer,breaks=layer_breakpoints,include.lowest=F,right=F,plot=FALSE)
		barplot(used_layer_hist$counts,
				names.arg=layer_breakpoints[1:(length(layer_breakpoints)-1)],
				col=layer.col,
				xlab=layer.xlab,
				ylab="Counts",
				main=layer.main
		)
		box()
		
		### Radius Plot
		radius_breakpoints <- seq(0,(max(variable_lookup$radius)+radius_interval),
				radius_interval)
		all_radius_hist <- hist(variable_lookup$radius,breaks=radius_breakpoints,include.lowest=F,right=F,plot=FALSE)
		used_radius_hist <- hist(variable_lookup_used$radius,breaks=radius_breakpoints,include.lowest=F,right=F,plot=FALSE)
		# Adjust for counts:
		adjusted_counts <- used_radius_hist$counts/all_radius_hist$counts
		barplot(adjusted_counts,
				names.arg=radius_breakpoints[1:(length(radius_breakpoints)-1)],
				xlab="Radius",
				ylab="Chosen / Possible",
				ylim=c(0,1),
				col=radius.col,
				main=radius.main)
		box()
		
		### Azimuth Plot	
		# TODO: Orient the bins on 0.
		# TODO: Adapt http://stackoverflow.com/questions/17266780/wind-rose-with-ggplot-r
		azimuth_circular <- circular(variable_lookup_used$azimuth[variable_lookup_used$radius > 0],
				type="angles",units="degrees",
				zero=pi/2,rotation="clock")	
		rose.diag(azimuth_circular,bins=ceiling(360/azimuth_interval),
				col=azimuth.col,prop=azimuth.prop,ticks=azimuth.ticks,
				shrink=azimuth.shrink,
				main=azimuth.main)
	
		
	}
}

