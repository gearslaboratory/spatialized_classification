
extract_window_objects <- function(x,object,non_object,window_dim,verbose=FALSE)
{
	if(verbose) message("Extracting objects...")
	object_extract <- extract_windows(x,object,poly_extract=c("edge","inside"),window_dim=window_dim)

	if(verbose) message("Extracting non-objects...")
	non_object_extract <- extract_windows(x,non_object,poly_extract="all",window_dim=window_dim)

	if(verbose) message("Piecing data together...")
	
	
	varnames <- object_extract$variable_lookup$variable
	
	object_extract_predictors <- object_extract$extracted_windows[,
			names(object_extract$extracted_windows) %in% varnames,drop=FALSE]
	
	object_extract_response <- object_extract$extracted_windows[,"PolyExtractType",drop=FALSE]
	
	object_extract_cleaned <- cbind(object_extract_response,object_extract_predictors)
	
	non_object_extract_predictors <- non_object_extract$extracted_windows[,
			names(non_object_extract$extracted_windows) %in% varnames,drop=FALSE]
	
	non_object_extract_response <- non_object_extract$extracted_windows[,"PolyExtractType",drop=FALSE]
	
	non_object_extract_cleaned <- cbind(non_object_extract_response,non_object_extract_predictors)
	
	# Should probably add back in the other data at some point.
	
	object_training <- list()
	object_training$extracted_windows <- rbind(
			object_extract_cleaned,
			non_object_extract_cleaned)
	object_training$variable_lookup <- object_extract$variable_lookup
	
	return(object_training)
	
}