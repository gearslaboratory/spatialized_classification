#' Construct spatialized classification model.
#' 
#' @param training_data List. Output list from ?extract_windows.
#' @param response_variable Character. Column name of the categorical variable to be used.
#' @param training_data_filter_type Character. Extracted windows are rectangular, by default, but they can be adjusted with this setting.  The default is "circular".
#' @param training_data_filter_parameters Vector. Parameters for use with various training_data_filter_types. See Details.
#' @param additional_variables Currently unsupported.  Will eventually allow for indices to be created.
#' @param imbalance_correction Character. Can be "none", "balanced", "upSample", or "downSample".  See Details. Default is "upSample".
#' @param nsamples Numeric. Number of samples to use with balanced correction.  Default is 500.
#' @param model Character. Classification algorithm to use.  Default is randomForestSRC. Other algorithms are not currently supported.
#' @param ntree Numeric. randomForestSRC only. Number of trees to use.  Default is 1000. See ?rfsrc.
#' @param importance Logical. randomForestSRC only. Calculate variable importance?  Default is FALSE.  Setting to TRUE can significantly increase the computation time and memory use. See ?rfsrc.
#' @param nsplit Non-negative integer value. randomForestSRC only. See ?rfsrc.
#' @param mtry Numeric. randomForestSRC only. Number of variables randomly selected as candidates for each node split. See Details and ?rfsrc.
#' @param auto_varselect Logical. If TRUE, will optimize the number of variables used in the final model using the maximal subtree algorithm.  Default is TRUE.
#' @param auto_ntree Logical. If TRUE, will optimze the number of trees used in the final model based on the error_rate_change_threshold setting.  Default is TRUE.
#' @param conservative Logical. Used with auto_varselect. See ?rfsrc. Default=FALSE.
#' @param cores Numeric. randomForestSRC only. How many cores to use in the construction of the model.  Default is half of the available cores (floor(detectCores()/2)).
#' @param do.trace Numeric. Currently unsupported.
#' @param error_rate_change_threshold Numeric. If the change in error rate is less than this amount, the model is deemed "finished" if auto_ntree=T.  Default is 0.001 (1% change).
#' @param verbose Logical. 
#' @author Jonathan A. Greenberg (\email{spatial.tools@@estarcion.net})
#' @details Constructs a classification model using the spatial-spectral data extracted from ?extract_windows. Currently, the model
#' uses randomForestSRC as the "engine".  The algorithm proceeds through several steps:
#' 1) If the user wishes to apply a filter to the extracted windows (e.g. only use pixels within a circle around the central pixel),
#' the algorithm first pre-filters the data.  Currently, only a circular filter is supported (more will be added later), and the radius
#' of this circle is set based on training_data_filter_parameters (training_data_filter_parameters=0 will use only the central pixel, for instance),
#' or, if missing, will be set the smaller of the two window_dim divided by 2.  
#' 2) Random Forests don't behave well with highly imbalanced datasets (more pixels belonging to one class than another), so we have included three
#' imbalancing approaches: "upSample" (default): the minority classes are sampled with replacement up to the size of the majority class, "downSample":
#' the majority classes are sampled without replacement to the size of the minority class, and "balanced": a predetermined sample size is set by
#' the "nsamples" parameter, and all classes with samples smaller than nsamples are upsampled, and all classes larger than nsamples are downsampled. 
#' We currently recommend using "upSample" for smaller training datasets, and "balanced" for much larger datasets.
#' 3) The model is run using all variables not filtered, with the number of trees set to "ntree".  If auto_varselect==TRUE, mtry will be set 
#' by default to the number of variables ^ 4/5.  If not, it will set by default to number of variables ^ 1/2.
#' 4) If auto_varselect=TRUE, the maximal subtree method (?max.subtree) will be used to select the variables needed in the model.  
#' Once this subset is determined, the model is rerun with the reduced variables set.
#' 5) If auto_ntree=TRUE, the error rate of the model will be examined to determine what number of trees is needed to realize an error rate change
#' of less than the error_rate_change_threshold. The model is re-run with only this number of trees.  If the model can't reduce the number of trees,
#' a warning will be thrown.  The model should probably be re-run with ntree set to a higher number.
#' 6) The final model is then returned to the user.
#' 
#' Key optimization settings: 
#' imbalance_correction = "balanced" and set a reasonable number of nsamples
#' importance = FALSE
#' auto_varselect = TRUE
#' auto_ntree = FALSE
#' cores = detectCores()
#' 
#' @examples
#' # Load extracted 5x5 windows from a set of polygons.
#' data("extracted_training_polys")
#'  
#' # Create model:
#' spatialized_classification_model <- spatialized_classification(training_data=extracted_training_polys,
#' 				response="class",training_data_filter_type="circular",model="randomForestSRC",importance=FALSE,
#' 				cores=1,imbalance_correction="upSample",
#' 				auto_varselect=TRUE,auto_ntree=TRUE)
#' 
#' @import randomForestSRC caret parallel
#' @export

spatialized_classification <- function(
		training_data,response_variable,
		training_data_filter_type="circular",
		training_data_filter_parameters=NULL,
		additional_variables=NULL,
		imbalance_correction="upSample", # Artificially balance the training data.
		nsamples=500, # Used with imbalanced_correction = "balanced"
		model="randomForestSRC",
		# randomForestSRC parameters
		ntree=1000,
		importance=FALSE,
		nsplit=1,
		mtry=NULL, # sqrt(p) unless auto_varselect=TRUE, in which case p^(4/5).
		auto_varselect=TRUE, # Try to reduce the number of variables used.
		auto_ntree=TRUE, # Try to reduce the number of trees to produce.
		conservative=FALSE,
		error_rate_change_threshold=0.1, # 0.1% change threshold.
		cores=floor(detectCores()/2),
		do.trace=0, # Won't do anything if cores > 0
		# Execution parameters:
		verbose=FALSE)
{
	
#	browser()
	if(missing(response_variable))
	{
		stop("Please assign a response variable.")
	}
	
	if(auto_ntree) err.block=1 else err.block=NULL
	
	response <- training_data$extracted_windows[,response_variable,drop=FALSE]
	if(class(response[[1]])=="character")
		response[[1]] <- as.factor(response[[1]])
	
	
	# browser()
	
	if(training_data_filter_type != "none")
	{
		if(training_data_filter_type == "circular")
		{
			if(is.null(training_data_filter_parameters))
			{
				# Fit to minimum axis:
				circle_radius = min(training_data$window_dim)/2
			} else
			{
				circle_radius = training_data_filter_parameters
			}
			predictor_variable_ids <- training_data$variable_lookup$variable[training_data$variable_lookup$radius <= circle_radius]
			
		}
#		print(circle_radius)
		
	} else
	{
		predictor_variable_ids <- training_data$variable_lookup$variable
	}
	
	predictors <- training_data$extracted_windows[,predictor_variable_ids,drop=FALSE]
	
	
	# Not functional yet.
	if(!is.null(additional_variables))
	{
		if("center_ratio" %in% additional_variables)
		{
			center_variables <- training_data$variable_lookup[training_data$variable_lookup$radius==0,]
			non_center_variables <- training_data$variable_lookup[training_data$variable_lookup$radius!=0,]
			cr_data_frame <- foreach(i=unique(center_variables$layer),.combine=cbind) %do%
					{
						#print(i)
						center_variable = center_variables[center_variables$layer==i,]$variable
						#print(center_variable)
						center_variable_index <- names(training_data$extracted_windows)==center_variable
						#print(center_variable_index)
						# center_variable_index <- as.numeric(rownames(center_variable))
						variables_in_layer <- non_center_variables[non_center_variables$layer==i,]$variable
						#print(variables_in_layer)
						# variables_in_layer_index <- as.numeric(rownames(variables_in_layer))
						# cr_variable_names <- paste(variables_in_layer$variable,"CR",sep="")
						cr_layer <- foreach(j=variables_in_layer,.combine=cbind)	 %do%
								{
									#print(j)
									temp_variable_index <- names(training_data$extracted_windows)==j
									temp_cr <- as.data.frame(
											(training_data$extracted_windows[,temp_variable_index]-
														training_data$extracted_windows[,center_variable_index])/
													(training_data$extracted_windows[,temp_variable_index]+
														training_data$extracted_windows[,center_variable_index])
									)
									names(temp_cr) <- paste(j,"CR",sep="")
									return(temp_cr)
								}
						return(cr_layer)
					}
			
		}
		
	}
	
	
#	browser()
	
	if(imbalance_correction=="balanced")
	{
		
		library("caret")	
		# Class counts:
		nperclass <- table(response[[1]])
		if(nsamples > max(nperclass))
		{
			nsamples = max(nperclass)
		}
		
		if(verbose) message(paste("Balancing the data with an N of ",nsamples,"...",sep=""))
		
		
		classes_to_upsample <- names(nperclass)[nperclass <= nsamples]
		classes_to_downsample <- names(nperclass)[nperclass > nsamples]
		
		# Dummy response:
		dummy_category <- basename(tempfile(pattern="dummyvar"))
		dummy_category_data.frame <- as.data.frame(rep(dummy_category,nsamples))
		names(dummy_category_data.frame) <- names(response)
		
		# Dummy predictors
		dummy_predictors_data.frame <- as.data.frame(matrix(nrow=nsamples,ncol=ncol(predictors)))
		names(dummy_predictors_data.frame) <- names(predictors)
		
		### Upsampled dataset:
		upsample_ids <- response[[1]] %in% classes_to_upsample
		upsample_response <- response[upsample_ids,,drop=FALSE]
		levels(upsample_response[[1]]) <- c(levels(upsample_response[[1]]),dummy_category)
		upsample_response <- rbind(upsample_response,dummy_category_data.frame)
		upsample_predictors <- rbind(predictors[upsample_ids,,drop=FALSE],dummy_predictors_data.frame)
		
		training_data_upsampled <- upSample(x=upsample_predictors,y=upsample_response[[1]],yname=response_variable)
		
		# Remove the dummy rows:
		upsampled_response <- training_data_upsampled[,dim(training_data_upsampled)[2],drop=FALSE]
		upsampled_response_dummy_ids <- upsampled_response[[1]]!=dummy_category
		upsampled_response <- upsampled_response[upsampled_response_dummy_ids,,drop=FALSE]
		upsampled_predictors <- training_data_upsampled[upsampled_response_dummy_ids,-dim(training_data_upsampled)[2],drop=FALSE]
		
		# Downsampled dataset:
		downsample_ids <- response[[1]] %in% classes_to_downsample
		downsample_response <- response[downsample_ids,,drop=FALSE]
		levels(downsample_response[[1]]) <- c(levels(downsample_response[[1]]),dummy_category)
		downsample_response <- rbind(downsample_response,dummy_category_data.frame)
		downsample_response[[1]] <- factor(downsample_response[[1]])
		downsample_predictors <- rbind(predictors[downsample_ids,,drop=FALSE],dummy_predictors_data.frame)
		
		training_data_downsampled <- downSample(x=downsample_predictors,y=downsample_response[[1]],yname=response_variable)
		
		# Remove the dummy rows:
		downsampled_response <- training_data_downsampled[,dim(training_data_downsampled)[2],drop=FALSE]
		downsampled_response_dummy_ids <- downsampled_response[[1]]!=dummy_category
		downsampled_response <- downsampled_response[downsampled_response_dummy_ids,,drop=FALSE]
		downsampled_predictors <- training_data_downsampled[downsampled_response_dummy_ids,-dim(training_data_downsampled)[2],drop=FALSE]
		
		# Merge them all back together now:
		#response_name <- names(response)
		#predictor_names <- names(predictors)
		response <- rbind(upsampled_response,downsampled_response)
		response[[1]] <- factor(response[[1]])
		predictors <- rbind(upsampled_predictors,downsampled_predictors)
		
	}
	
#	if(imbalance_correction=="SMOTE")
#	{
#		# NOT WORKING
#		library("DMwR")
#		browser()
#		merged_training_data <- cbind(predictors,response)
#		names(merged_training_data)[ncol(merged_training_data)] <- "response"
#		
#		smoted <- SMOTE(response ~ .,data=merged_training_data)
#	}
	
	
	if(imbalance_correction=="upSample")
	{
		if(verbose) message("Upsampling the data...")
		
		library("caret")	
		training_data_upsampled <- upSample(x=predictors,y=response[[1]],yname=response_variable)
		response <- training_data_upsampled[,dim(training_data_upsampled)[2],drop=FALSE]
		predictors <- training_data_upsampled[,-dim(training_data_upsampled)[2],drop=FALSE]
		training_data_upsampled <- NULL
	}
	
	if(imbalance_correction=="downSample")
	{
		if(verbose) message("Downsampling the data...")
		
		library("caret")	
		training_data_downsampled <- downSample(x=predictors,y=response[[1]],yname=response_variable)
		response <- training_data_downsampled[,dim(training_data_downsampled)[2],drop=FALSE]
		predictors <- training_data_downsampled[,-dim(training_data_downsampled)[2],drop=FALSE]
		training_data_downsampled <- NULL
	}
	
	# Haven't tested classic randomForest
#	if(model=="randomForest")
#	{
#		if(verbose) message("Starting randomForest...")
#		
#		nworkers = getDoParWorkers()
#		ntree_per_worker = floor(ntree/nworkers)
#		ntree_worker_vector <- c(rep(ntree_per_worker,nworkers-1),ntree_per_worker+ntree-(ntree_per_worker*nworkers))
#		
	##		training_data_cleaned <- training_data$extracted_windows[,-c(1,3)]
#		# training_data_cleaned <- training_data$extracted_windows
#		library(randomForest)
#		# Parallel version
#		# From http://cran.r-project.org/web/packages/foreach/vignettes/foreach.pdf
#		sfQuickInit(cpus=cores)
#		
#		spatialized_classification_model <- foreach(ntree=ntree_worker_vector, .combine=combine, .packages='randomForest') %dopar%
#				{
#					randomForest(y=response[[1]],x=predictors,ntree=ntree,
#							importance=importance)
#				}
#		
#		sfQuickStop()
#		
#		spatialized_classification_importance <- as.data.frame(importance(spatialized_classification_model))
#		
#		spatialized_classification_importance_joined <- merge(spatialized_classification_importance,training_data$variable_lookup,
#				by.x="row.names",by.y="variable")
#		
#		spatialized_classification_out <- list(model=spatialized_classification_model,importance=spatialized_classification_importance_joined)
#	}
	
	if(model=="randomForestSRC")
	{
		if(verbose) message("Starting randomForestSRC...")
		if(is.null(mtry) && auto_varselect)
		{
			if(auto_varselect)
			{
				mtry <- ceiling(ncol(predictors)^(4/5))
			} else
			{
				mtry <- ceiling(ncol(predictors)^(1/2))
			}
		}
		
		library("randomForestSRC")
		options(rf.cores = cores,mc.cores=cores)
		# Need to see if there is a better way to do this:
		merged_training_data <- cbind(response,predictors)
		names(merged_training_data)[1] <- "response"
		
		if(!importance) 
		{	
			rfsrc_importance="none"
		} else
		{
			rfsrc_importance = c("permute", "random", "permute.ensemble", "random.ensemble", "none")
		}
		
		if(verbose) message("Running initial randomForestSRC...")
		spatialized_classification_model <- rfsrc(response ~ .,data=merged_training_data,
				importance=rfsrc_importance,ntree=ntree,nsplit=nsplit,do.trace=do.trace,
				mtry=mtry,err.block=err.block)
		
		if(auto_varselect)
		{
			#	browser()
			if(verbose) message("Auto-optimizing: running variable selection...")
			spatialized_classification_maxsubtree <- max.subtree(spatialized_classification_model,conservative=conservative)
			
			if(length(spatialized_classification_maxsubtree$topvars)!=0)
			{
				if(verbose) 
				{
					message("Auto-optimizing: rerunning with selected variables...")
					message(paste("Original number of variables:",ncol(predictors)))
					message(paste("Reduced set of variables:",length(spatialized_classification_maxsubtree$topvars)))
					
				}
				merged_training_data <- cbind(response,predictors[spatialized_classification_maxsubtree$topvars])
				names(merged_training_data)[1] <- "response"
				
				# Recalculate mtry back to sqrt(p)...
				# mtry <- ceiling((ncol(merged_training_data)-1)^(1/2))
				
				spatialized_classification_model <- rfsrc(response ~ .,data=merged_training_data,
						importance=rfsrc_importance,ntree=ntree,nsplit=nsplit,do.trace=do.trace,
						mtry=mtry,imbalance_correction="none",err.block=err.block)
				
				
			} else
			{
				if(verbose) message("No variables were selected, try increasing mtry if variable selection is wanted...")
			}			
		}
		
		if(auto_ntree)
		{
			# Now try to reduce the number of trees...
			model_error_rate <- spatialized_classification_model$err.rate
			ntree_optimized <- max(which(abs(diff(model_error_rate[,1])*100) > error_rate_change_threshold))+1
			if(ntree_optimized >= ntree)
			{
				message("ntree is already reduced.  There is a chance your initial ntree was set too low or your error_rate_change_threshold was too high...")
			} else
			{
				# Finally finish the model...
				if(verbose) message(paste("Auto-optimizing: reducing the number of trees to: ",ntree_optimized," trees...",sep=""))
				spatialized_classification_model <- rfsrc(response ~ .,data=merged_training_data,
						importance=rfsrc_importance,ntree=ntree_optimized,nsplit=nsplit,do.trace=do.trace,
						mtry=mtry,imbalance_correction="none",err.block=err.block)	
			}
		}
		
		spatialized_classification_out=list(model=spatialized_classification_model,
				variable_lookup=training_data$variable_lookup,window_dim=training_data$window_dim,
				response_variable=response_variable)
	}
	
	
	return(spatialized_classification_out)
	
}