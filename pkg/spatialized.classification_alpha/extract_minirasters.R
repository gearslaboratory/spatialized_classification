

extract_minirasters <- function(x,y,window_dim,
		return_minirasterstrip=TRUE,
		return_minirasterstrip_vector=TRUE,
		minirasterstrip_divider=0,
		output_directory,
		output_prefix,
		format="ENVI")
{
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
		
		# browser()
		
		if(is.character(output_directory))
		{
			output_filename_prefix <- file.path(output_directory,output_prefix)
			output_filename_ids <- add_leading_zeroes(seq(y),max_number=length(y))
			output_filenames <- paste(output_filename_prefix,"_V",output_filename_ids,sep="")
		} else
		{
			output_filenames=rep("",length(y))
		}
		
		extracted_minirasters <- foreach(i = seq(y),.packages=c("rgdal","raster")) %dopar%
				{
					#	print(i)
					#			browser()
					
					y_point <- y[i,]
					
					central_cell <-  cellFromXY(x,y_point)
					
					if(is.na(central_cell))
					{
						
						# Point is removed from the output.
						return(NULL)
						
					} else
					{
						
						
						central_rowCol <- rowColFromCell(x,central_cell)
						topleft_rowCol <- central_rowCol-window_buffer
						bottomright_rowCol <- central_rowCol+window_buffer
						
						crop_extent <- extent(x,topleft_rowCol[1],bottomright_rowCol[1],topleft_rowCol[2],bottomright_rowCol[2])
						
						miniraster <- crop(x=x,y=crop_extent,filename=output_filenames[i],format=format)
						
					}
					return(miniraster)
				}
		
		
		
		extracted_minirasters_null <- sapply(extracted_minirasters,is.null)
		extracted_minirasters_wrongdims <- sapply(extracted_minirasters,
				function(X,window_dim) 
				{ 
					!(all(dim(X)[1:2]==window_dim))
				}
				,window_dim=window_dim)
		
		
		extracted_minirasters_null <- extracted_minirasters_null | extracted_minirasters_wrongdims
		
		extracted_minirasters[extracted_minirasters_null] <- NULL
		
		#	browser()
		if(return_minirasterstrip)
		{
			# Strip all the projections
			extracted_minirasters_modded <- foreach(i=seq(extracted_minirasters),.packages=c("rgdal","raster")) %do%
					{
						#	print(i)
						temp_miniraster <- extracted_minirasters[[i]]
						projection(temp_miniraster) <- NULL
						extent(temp_miniraster) <- c(0,window_dim[1],
								(i-1)*window_dim[2]+((i-1)*minirasterstrip_divider),
								(i-1)*window_dim[2]+window_dim[2]+((i-1)*minirasterstrip_divider))
						#	print(class(temp_miniraster))
						return(temp_miniraster)
					}
			
			minirasterstrip <- do.call(merge,extracted_minirasters_modded)
		} else
		{
			minirasterstrip <- NULL
		}
#		browser()
		if(return_minirasterstrip && return_minirasterstrip_vector)
		{
			new_x <- rep(window_dim[1]/2,length(extracted_minirasters))
			new_y <- (seq(extracted_minirasters)-1)*window_dim[2]+(window_dim[2]/2)+((seq(extracted_minirasters)-1)*minirasterstrip_divider)
			
			y_subset_data <- y[!extracted_minirasters_null,]@data
			coordinates_y_subset <- cbind(new_x,new_y)
			
			minirasterstrip_vector <- SpatialPointsDataFrame(coordinates_y_subset,data=y_subset_data)
		} else
		{
			minirasterstrip_vector <- NULL
		}
		
	}
	
	output_miniraster <- list(minirasters=extracted_minirasters,minirasterstrip=minirasterstrip,minirasterstrip_vector=minirasterstrip_vector)
	return(output_miniraster)
}