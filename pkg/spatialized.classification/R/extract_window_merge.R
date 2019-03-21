#' @import plyr
#' @export

extract_window_merge <- function(extract_window_list)
{
#	library("plyr")
	
	# Precheck that the variables and lookup are the same:
	reference_variable_lookup <- extract_window_list[[1]]$variable_lookup
	reference_window_dim <- extract_window_list[[1]]$window_dim

	# TODO: window_dim check
	
	variable_lookup_check <- all(sapply(extract_window_list,
			function(X,reference_variable_lookup)
			{
				return(identical(X$variable_lookup,reference_variable_lookup))
			},reference_variable_lookup=reference_variable_lookup))

	if(!variable_lookup_check)
	{
		stop("Variable lookups don't match, can't merge.")
	}
	
	extracted_windows_list <- lapply(extract_window_list,function(X) { return(X$extracted_windows)})
	
	merged_extracted_windows <- rbind.fill(extracted_windows_list)

	
	extract_window_merged <- list(extracted_windows=merged_extracted_windows,
			variable_lookup=reference_variable_lookup,window_dim=reference_window_dim
			)
	
	return(extract_window_merged)
}