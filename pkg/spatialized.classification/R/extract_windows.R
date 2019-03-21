#' Local window extraction to create training and testing data for spatialized classification.
#' @param x Raster*. A Raster* used to extract data from.
#' @param y Spatial*. A SpatialPolygonsDataFrame or SpatialPointsDataFrame object used to extract the data.
#' @param window_dim Numeric vector. Two-element vector representing the width and height in pixel units of the window.  Default is 1,1 (a single pixel).
#' @param poly_extract Logical. FALSE currently unsupported, leave as TRUE.
#' @param remove_outside_bbox_points Logical. Remove any points falling outside of the raster bbox.
#' @param verbose Logical. Verbose execution.  Default = FALSE.
#' @author Jonathan A. Greenberg (\email{spatial.tools@@estarcion.net})
#' @details For each pixel falling within a polygon, or under a point, a window around that pixel as given by the window_dim
#' parameter is extracted.  The dataframe variable names are labelled VN, where N is the Nth element of the array surrounding
#' the central pixel, moving line by line, and then band by band.  The output object is a list containing three elements,
#' the "extracted_windows" for use with the classifier, the "variable_lookup" which links the variable names against
#' the local radius, azimuth and band, and "window_dim" which is the same as the input "window_dim". 
#' 
#' For large extractions, if a foreach engine is registered, the extraction will run in parallel, one polygon per worker.
#' 
#' @examples
#' # Extract 5x5 windows from a set of polygons.
#' training_raster <- setMinMax(
#'	brick(system.file("external/tahoe_highrez.tif", package="spatial.tools")))
#' 
#' training_polys <- readOGR(dsn=system.file("external"),
#' 		layer="spatialized.classification.training_v01")
#' 
#' window_dim = c(5,5)
#' extracted_training_polys <- extract_windows(x=training_raster,y=training_polys,window_dim=window_dim)
#' 
#' @import raster rgdal foreach
#' @export
	
extract_windows <- function(x, y,
		window_dim=c(1,1),
		poly_extract="all",
#		return_variable_lookup=TRUE,
		remove_outside_bbox_points=TRUE,
#		window_shape="rectangular",
		verbose=FALSE)
{
	# Helper function:
	
	return_variable_lookup=TRUE
	
	
	in.bbox <- function(x,bbox)
	{
		# Coerce various formats to bbox matrix format:
		if(class(bbox)!="matrix" && !identical(colnames(bbox),c("min","max")))
		{
			# Convert here
		}
		
		if(class(x)!="matrix" && !identical(x,c("min","max")))
		{
			# Convert here
		}
		
		lx_in <- x[1,1] >= bbox[1,1] && x[1,1] <= bbox[1,2]
		rx_in <- x[1,2] >= bbox[1,1] && x[1,1] <= bbox[1,2]
		ly_in <- x[2,1] >= bbox[2,1] && x[2,1] <= bbox[2,2]
		uy_in <- x[2,2] >= bbox[2,1] && x[2,2] <= bbox[2,2]
		
		return(lx_in & rx_in & ly_in & uy_in)
		
	}
	
	
	# TODO: confirm window_dim are odd, not even
	
	projection_x <- projection(x)
	projection_y <- projection(y)
	
	# Fix windows to make sure they are 2-d
	if(length(window_dim)==1) window_dim <- c(window_dim,window_dim)
	
	# Buffer is 1/2 window size rounded down:
	window_buffer <- floor(window_dim/2)
	
	# POINT MODE:
	if(class(y)=="SpatialPointsDataFrame" || class(y)=="SpatialPoints")
	{
		# REPROJECT HERE
		# browser()
		y <- spTransform(x=y,CRSobj=CRS(projection_x))
		
		extracted_windows <- foreach(i = seq(y),.combine=rbind,.packages=c("rgdal","raster")) %dopar%
				{
					# print(i)
					#			browser()
					
					y_point <- y[i,]
					
					central_cell <-  cellFromXY(x,y_point)
					if(is.na(central_cell))
					{
						# Point is outside raster
						if(remove_outside_bbox_points)
						{
							# Point is removed from the output.
							return(NULL)
						}
						else
						{
							# Return the extracted values as NA.
							all_extraction_dataframe <- cbind(
									as.data.frame(y_point),
									data.frame(PolyExtractType=NA),
									as.data.frame(t(rep(NA,prod(window_dim)*nlayers(x))),
											row.names = NULL)
							)
						}
					} else
					{
						
						
						central_rowCol <- rowColFromCell(x,central_cell)
						topleft_rowCol <- central_rowCol-window_buffer
						
						
						all_extraction_raw <- list(getValuesBlock(x=x,row=topleft_rowCol[1],
										nrows=window_dim[1],col=topleft_rowCol[2],ncols=window_dim[2]))	
						
						
						all_extraction_raw <- getValuesBlock(x=x,row=topleft_rowCol[1],
								nrows=window_dim[1],col=topleft_rowCol[2],ncols=window_dim[2])
						
						all_extraction_vector <- as.numeric(all_extraction_raw)
						
						all_extraction_dataframe <- as.data.frame(t(all_extraction_vector))
					}
					# browser()
					
					# Now we need to unspool each window, note that this reads each band after the previous one, 
					# presumably from upper left to lower right.  The central pixel should be in the center position.
					
					#				all_extraction_vector_list <- lapply(all_extraction_raw,function(x) as.numeric(x))
					
					# Now convert it to a vector:
					#				all_extraction_dataframe <- as.data.frame(do.call(rbind,all_extraction_vector_list))
					
					# Finally, append the data if the original Spatial*DataFrame had some, and the extraction
					# type.
					# TODO: If no *DateFrame, this will throw an error.
					#				all_extraction_dataframe <- cbind(y_data,data.frame(PolyExtractType="all"),all_extraction_dataframe)
					
					all_extraction_dataframe <- cbind(as.data.frame(y_point),data.frame(PolyExtractType="point"),all_extraction_dataframe,row.names = NULL)
					
					# Now we fuse all the data frames together.  Note that setting the non-returned
					# items to NULL helps this step:
					return(all_extraction_dataframe)
				}
	}
	
	# POLYGON MODE
	if(class(y)=="SpatialPolygonsDataFrame" || class(y)=="SpatialPolygons")
	{
		# Now we'll loop through each polygon.
		# If a cluster is registered, this will run in parallel.
		extracted_windows <- foreach(i = seq(y),.combine=rbind,.packages=c("rgdal","raster")) %dopar%
				{
					# print(i)
					# if(i==144) browser()
					# TODO: add in support for SpatialPoints*
					
					# Assume SpatialPolygonsDataFrame for now
					# Extract a single polygon to work with:
					y_polygons <- y@polygons[[i]]
					
					# And its data:
					y_data <- y@data[i,,drop=FALSE]
					
					# Coerce back to a SpatialPolygons and then a SpatialPolygonsDataFrame:
					y_polygons_single_sp <- SpatialPolygons(list(y_polygons),
							proj4string = CRS(projection_y))
					y_polygons_single_sp <- SpatialPolygonsDataFrame(y_polygons_single_sp,
							data=y_data)
					
					# Reproject into raster's coordinate system:
					y_polygons_single_sp_transform <- spTransform(
							y_polygons_single_sp,CRS(projection_y))
					
					# Get bbox of the polygon:
					y_polygons_bbox <- bbox(y_polygons_single_sp_transform)
					
					# Buffer the bbox if need be:
					if(length(window_buffer)==1) window_buffer <- c(window_buffer,window_buffer)
					
					if(any(window_buffer != 0))
					{
						# Expand the bbox by the buffer window (in geographic coordinates)
						# First convert window_buffer (pixel units) to geographic:
						
						buffer_geog <- res(x) * window_buffer
						
						# Now expand the bbox on all sides:
						y_polygons_bbox[1,1] <- y_polygons_bbox[1,1] - buffer_geog[1]
						y_polygons_bbox[1,2] <- y_polygons_bbox[1,2] + buffer_geog[1]
						y_polygons_bbox[2,1] <- y_polygons_bbox[2,1] - buffer_geog[2]
						y_polygons_bbox[2,2] <- y_polygons_bbox[2,2] + buffer_geog[2]
					}
					
					# Need to test for in.bboxness:
					if(in.bbox(y_polygons_bbox,bbox(x)))
					{
						
						
						# Now we extract the miniraster (the image surrounding the poly):
						miniraster <- crop(x,y_polygons_bbox)
						
						# Rasterize the polygon to the miniraster's dimensions:
						miniraster_all_poly_mask <- rasterize(y_polygons_single_sp_transform,miniraster)
						
						# Now convert all the places that were in the polygon to points:
						miniraster_all_poly_points <- rasterToPoints(miniraster_all_poly_mask,spatial=FALSE)
						
						# If we need the edge or the inside, we need to rasterize the edge
						# of the polygon:
						if("edge" %in% poly_extract || "inside" %in% poly_extract)
						{
							# Converting to a line didn't work, so...
							
							# We are going to use a focal window to determine the inside.
							# All this is doing is taking a 3x3 window, if the window finds even a single NA
							# (the pixel is at the edge of the poly) it returns NA, otherwise a 1, because
							# na.rm=FALSE.
							miniraster_inside_poly_mask <- focal(miniraster_all_poly_mask,w=matrix(1,nrow=3,ncol=3),
									na.rm=FALSE)
							
							# Convert this mask to a set of inside points (matrix format):
							miniraster_inside_poly_points <- rasterToPoints(miniraster_inside_poly_mask,spatial=FALSE)
						}
						
						# If we want inside, we now use the previous two outputs to figure this out:
						if("edge" %in% poly_extract)
						{
							# TODO: this would probably be faster to use some raster math on the 
							# inside and all masks.
							
							# First, make a single "ID" field to make searching for differences easier.
							# First the entire poly:
							miniraster_all_poly_points_ids <- paste("x_",miniraster_all_poly_points[,"x"],
									"_y_",miniraster_all_poly_points[,"y"],sep="") 
							# Now the inside:
							miniraster_inside_poly_points_ids <- paste("x_",miniraster_inside_poly_points[,"x"],
									"_y_",miniraster_inside_poly_points[,"y"],sep="") 
							
							# Now we just need to determine what differs between the two (those will be the
							# inside points).  First we determine the IDs:
							miniraster_edge_poly_points_ids <- setdiff(miniraster_all_poly_points_ids,
									miniraster_inside_poly_points_ids)
							
							# Now use these IDs to extract from the entire polygon points:
							miniraster_edge_poly_points <- miniraster_all_poly_points[
									miniraster_all_poly_points_ids %in% 
											miniraster_edge_poly_points_ids,]	
						}
						
						# Now we are ready to do the extractions:
						
						# If we just want the entire poly processed:
						if("all" %in% poly_extract)
						{	
							# Determine the row and col of each point from the coordinates:
							all_rowCol <- as.data.frame(
									rowColFromCell(miniraster,
											cellFromXY(miniraster,miniraster_all_poly_points[,c("x","y")])))
							
							# Now we'll use getValuesBlock to pull out the window around each point, using
							# lapply to loop through them:
							all_extraction_raw <- lapply(split(all_rowCol,rownames(all_rowCol)),function(x,miniraster,window_buffer)
									{
										# TODO: if start_row or start_col is <= 0, or if start_row+nrows 
										# or start_col+ncol is larger than the miniraster, we need to 
										# return some NAs instead (the window runs off the image).
										
										start_row <- x$row-window_buffer[1]
										nrows <- window_buffer[1]*2+1
										
										start_col <- x$col-window_buffer[2]
										ncols <- window_buffer[2]*2+1
										
										window <- getValuesBlock(x=miniraster,row=start_row,nrows=nrows,
												col=start_col,ncols=ncols)
										
										return(window)
										
									},
									miniraster=miniraster,window_buffer=window_buffer)
							
							# Now we need to unspool each window, note that this reads each band after the previous one, 
							# presumably from upper left to lower right.  The central pixel should be in the center position.
							
							all_extraction_vector_list <- lapply(all_extraction_raw,function(x) as.numeric(x))
							
							# Now convert it to a vector:
							all_extraction_dataframe <- as.data.frame(do.call(rbind,all_extraction_vector_list))
							
							# Finally, append the data if the original Spatial*DataFrame had some, and the extraction
							# type.
							# TODO: If no *DateFrame, this will throw an error.
							all_extraction_dataframe <- cbind(y_data,data.frame(PolyExtractType="all"),all_extraction_dataframe,row.names = NULL)
							
						} else
						{
							all_extraction_dataframe <- NULL
						}
						
						# Let's do the same for the edge pixels (same notes as above):
						if("edge" %in% poly_extract)
						{	
							# Perform the extraction:
							edge_rowCol <- as.data.frame(rowColFromCell(miniraster,
											cellFromXY(miniraster,miniraster_edge_poly_points[,c("x","y")])))
							
							edge_extraction_raw <- lapply(split(edge_rowCol,rownames(edge_rowCol)),function(x,miniraster,window_buffer)
									{
										start_row <- x$row-window_buffer[1]
										nrows <- window_buffer[1]*2+1
										
										start_col <- x$col-window_buffer[2]
										ncols <- window_buffer[2]*2+1
										
										window <- getValuesBlock(x=miniraster,row=start_row,nrows=nrows,
												col=start_col,ncols=ncols)
										
										return(window)
										
									},
									miniraster=miniraster,window_buffer=window_buffer)
							
							
							# Now we need to unspool each window, note that this reads each band after the previous one, 
							# presumably from upper left to lower right.  The central pixel should be in the center position.
							
							edge_extraction_vector_list <- lapply(edge_extraction_raw,function(x) as.numeric(x))
							
							# Now convert it to a vector:
							edge_extraction_dataframe <- as.data.frame(do.call(rbind,edge_extraction_vector_list))
							
							# Finally, append the data if the original Spatial*DataFrame had some:
							edge_extraction_dataframe <- cbind(y_data,data.frame(PolyExtractType="edge"),edge_extraction_dataframe,row.names = NULL)
						} else
						{
							edge_extraction_dataframe <- NULL
						}
						
						# And inside (same as above):
						if("inside" %in% poly_extract)
						{	
							inside_rowCol <- as.data.frame(rowColFromCell(miniraster,
											cellFromXY(miniraster,miniraster_inside_poly_points[,c("x","y")])))
							
							inside_extraction_raw <- lapply(split(inside_rowCol,rownames(inside_rowCol)),function(x,miniraster,window_buffer)
									{
										start_row <- x$row-window_buffer[1]
										nrows <- window_buffer[1]*2+1
										
										start_col <- x$col-window_buffer[2]
										ncols <- window_buffer[2]*2+1
										
										window <- getValuesBlock(x=miniraster,row=start_row,nrows=nrows,
												col=start_col,ncols=ncols)
										
										return(window)
										
									},
									miniraster=miniraster,window_buffer=window_buffer)
							
							# Now we need to unspool each window, note that this reads each band after the previous one, 
							# presumably from upper left to lower right.  The central pixel should be in the center position.
							
							inside_extraction_vector_list <- lapply(inside_extraction_raw,function(x) as.numeric(x))
							
							# Now convert it to a vector:
							inside_extraction_dataframe <- as.data.frame(do.call(rbind,inside_extraction_vector_list))
							
							# Finally, append the data if the original Spatial*DataFrame had some:
							inside_extraction_dataframe <- cbind(y_data,data.frame(PolyExtractType="inside"),inside_extraction_dataframe,row.names = NULL)
						} else
						{
							inside_extraction_dataframe <- NULL
						}
						
						# Now we fuse all the data frames together.  Note that setting the non-returned
						# items to NULL helps this step:
					} else
					{
						if(remove_outside_bbox_points)
						{
							# Point is removed from the output.
							return(NULL)
						}
						else
						{
							# Return the extracted values as NA.
							all_extraction_dataframe <- cbind(
									y_data,
									data.frame(PolyExtractType=NA),
									as.data.frame(t(rep(NA,prod(window_dim)*nlayers(x))),row.names = NULL)
							)
						}
						
					}
					window_data <- rbind(all_extraction_dataframe,edge_extraction_dataframe,inside_extraction_dataframe)
					return(window_data)
				}
	}
	
	
	if(return_variable_lookup)
	{	
		
		nlayers_x <- nlayers(x)
		window_x <- seq(window_dim[1])-ceiling(window_dim[1]/2)
		window_y <- seq(window_dim[2])-ceiling(window_dim[2]/2)
		
		x_matrix <- matrix(window_x,nrow=window_dim[1],ncol=window_dim[2],byrow=TRUE)
		y_matrix <- matrix(window_y,nrow=window_dim[1],ncol=window_dim[2],byrow=FALSE)
		
		radius_matrix <- sqrt(x_matrix^2 + y_matrix^2)
		rotate = function(mat) t(mat[nrow(mat):1,,drop=FALSE])
		
		azimuth_matrix <- rotate(rotate(rotate(atan2(y_matrix,x_matrix)*180/pi)))
		azimuth_matrix[azimuth_matrix < 0] <- (360+azimuth_matrix)[azimuth_matrix < 0]
		
		# For a single band, unspool and replicate times the number of bands:
		radius_vector = rep(as.numeric(radius_matrix),nlayers_x)
		azimuth_vector = rep(as.numeric(azimuth_matrix),nlayers_x)
		
		num_spatial <- length(radius_vector)/nlayers_x
		
		num_vars <- length(radius_vector)
		var_names <- paste("V",seq(num_vars),sep="")
		
		bands_vector <- rep(seq(nlayers_x),each=num_spatial)
		
		
		
		variable_lookup = data.frame(variable=var_names,radius=radius_vector,azimuth=azimuth_vector,layer=bands_vector,
				stringsAsFactors=FALSE)
#	
	
#		if(window_shape="elliptical")
#		{
#			browser()
#		}
#		
		return(list(extracted_windows=extracted_windows,variable_lookup=variable_lookup,window_dim=window_dim))
		
	} else
	{
		return(extracted_windows)	
	}
	
}




