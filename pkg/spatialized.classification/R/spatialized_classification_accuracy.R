#' @import caret
#' @export

spatialized_classification_accuracy <- function(sc,testing_data)
{
	library("caret")
	# library("psych")
	
	# Predict on the testing data:
	x_predict <- predict_spatialized_classification(x=testing_data,sc=sc,
			parallelmode="none",debugmode=FALSE,classonly=TRUE)
	
	x_reference <- testing_data$extracted_windows[,names(testing_data$extracted_windows)==sc$response_variable]
	
	# Kappa sometimes fails:
	sc_confusion <- confusionMatrix(x_predict,x_reference)
	# sc_kappa <- cohen.kappa(cbind(x_predict,x_reference))
	
	return(sc_confusion)
	
}