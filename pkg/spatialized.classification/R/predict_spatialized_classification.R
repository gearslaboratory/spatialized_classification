#' @import spatial.tools randomForestSRC
#' @export

predict_spatialized_classification <- function(x,sc,
		classonly=TRUE,
		parallelmode="none",
		cores=detectCores(),
		ratify_output=TRUE,
		debugmode=FALSE,
		verbose=FALSE)
{
#	browser()
	window_dim <- sc$window_dim
	
	
	if(parallelmode=="openmp")
	{
		registerDoSEQ()
		options(rf.cores=cores, mc.cores=cores)
	} else
	{
		options(rf.cores=1,mc.cores=1)
	}
	
	if(is.Raster(x))
	{
		# User rasterEngine to predict.
		
		sc_prediction_function <- function(x,sc,parallelmode,cores,verbose)
		{
			#	library("randomForestSRC")
#		options(rf.cores=0,mc.cores=1)
			if(parallelmode=="openmp")
			{
				options(rf.cores=cores, mc.cores=cores)
			} else
			{
				options(rf.cores=0,mc.cores=1)
			}
			
			# Need to distinguish randomForest from randomForestSRC
#		options(rf.cores=24,mc.cores=24)
			if(any(class(sc)=="character") && !("sc_loaded" %in% ls())) 
			{
				if(verbose) message("Loading model...")
				load(sc)
			}
			sc_loaded <- sc
			sc <- NULL
			
			
			
			if(all(window_dim==1))
			{
				x_df <- x
				dim(x_df) <- c(prod(dim(x)[1:2]),dim(x)[3])
				x_df <- as.data.frame(x_df)
			} else
			{
				
				x_df <- x
				dim(x_df) <- c(dim(x)[1],prod(dim(x)[2:3]))
				x_df <- as.data.frame(x_df)
				# add names here
			}
			
			
			# Get an index of complete cases:
			
			
			if(any(class(sc_loaded)=="rfsrc"))
			{
				x_df_complete_id <- complete.cases(x_df)
				
				x_predict <- array(dim=c(dim(x)[1],1,1))
				
				if(sum(x_df_complete_id) > 0)
				{
					
					x_df_complete <- x_df[x_df_complete_id,]
					# randomForestSRC
					if(verbose) message("rfsrc prediction started...")
#			system.time(
					browser()
					x_predict_raw <- predict(sc_loaded,newdata=x_df_complete,importance="none",
							membership=FALSE,var.used=FALSE,split.depth=FALSE)
#			)
					if(verbose) message("rfsrc prediction ended...")
					x_predict_class <- as.numeric(x_predict_raw$class)
					x_predict_raw <- NULL
					x_predict[x_df_complete_id,1,1] <- x_predict_class
				}
				#	x_predict <- array(x_predict_class,dim=c(dim(x)[1],1,1))
			} else
			{
				# randomForest
				x_predict <- as.numeric(predict(sc_loaded$model,newdata=x_df,type="response"))
				x_predict <- array(x_predict,dim=c(dim(x)[1],1,1))
			}
			
			# Fix for window_dim==1
			if(all(window_dim==1))
			{
				dim(x_predict) <- c(dim(x)[1:2],1)
			}
			
#			if(dim(x)[1]>2) {
#				browser()
#			}
			
			return(x_predict)
		}
		
		# Disabling classic randomForest:
#		if(any(class(sc)=="character"))
#		{
#			packages=c("randomForestSRC","randomForest")
#		} else
#		{
#			
#			if(any(class(sc)=="rfsrc"))
#			{
#				packages=c("randomForestSRC")
#				sc=sc$model
#			} else
#			{
#				packages=c("randomForest")
#				sc=sc$model
#			}
#		}
		
# rfsrc only:

		packages=c("randomForestSRC")
		sc=sc$model
		
#
		
		sc_prediction <- rasterEngine(x=x,fun=sc_prediction_function,
				args=list(sc=sc,parallelmode=parallelmode,cores=cores,verbose=verbose),
				window_dim=window_dim,processing_unit="chunk",debugmode=debugmode,
				.packages=packages,verbose=verbose,
#			blocksize=1,
				outbands=1,
				datatype="INT1U" #,blocksize=1
		)
		
		# browser()
		
		# Now we need to ratify the output:
		# Model levels:
		if(ratify_output)
		{
			model_levels <- levels(sc$class)
			model_lookup <- data.frame(ID=seq(model_levels),levels=model_levels)
			
			sc_prediction <- as.factor(round(sc_prediction))
			sc_prediction_rat <- levels(sc_prediction)[[1]]
			sc_prediction_rat <- merge(sc_prediction_rat,model_lookup,all.x=TRUE)
			levels(sc_prediction) <- sc_prediction_rat
		}
	} else
	{
		# (Typically) for accuracy assessment.
		if(any(class(sc)=="character") && !("sc_loaded" %in% ls())) 
		{
			if(verbose) message("Loading model...")
			load(sc)
		}
		sc_loaded <- sc
		sc <- NULL
		if(verbose) message("rfsrc prediction started...")
		
		x_predict_raw <- predict(sc_loaded$model,newdata=x$extracted_windows,importance="none",
				membership=FALSE,var.used=FALSE,split.depth=FALSE)
		if(verbose) message("rfsrc prediction ended...")
		if(classonly)
		{
			sc_prediction <- x_predict_raw$class
		} else
		{
			sc_prediction <- x_predict_raw
		}
		
	}
	
	
#	sc_prediction_vector <- rasterToPolygons(sc_prediction,dissolve=TRUE)
	return(sc_prediction)
	
}