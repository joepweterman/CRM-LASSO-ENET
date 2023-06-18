predict.crm_lasso <- function (object, newdata = NULL, ...) {
  # check inputs
  if (missing(object)) {
    stop("argument 'object' is missing, with no default")
  }
  if (class(object) != "crm_lasso") {
    stop("argument 'object' must be of class 'crm_lasso'")
  }

  # predict response values
  if (is.null(newdata)) {
    return(object$fitted.values)
  } else {
    X <- model.matrix(delete.response(object$terms), newdata)
    predicted.values <- as.numeric(X %*% object$coefficients)
    names(predicted.values) <- rownames(newdata)
    return(predicted.values)
  }
}
