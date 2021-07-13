#' @import logger
NULL

update_appender <- function() {
	logging_file <- RMassBank.env$logging_file
	if (!is.null(logging_file)) {
		appender_obj <- logger::log_appender()
		if (as.character(appender_obj)[1] != "logger::appender_file") {
			appender_obj <- logger::appender_file(logging_file)
			logger::log_appender(appender_obj)
		}
	}
}

#' Update logging file and pass arguments to logger::log_info
#'
#' The logging file to be used must be specified by the user in the \code{logging_file} field of \code{settings.ini}
#' @seealso logger::log_info
#' @author pstahlhofen
#' @export
log_info <- function(...) {
	update_appender()
	logger::log_info(...)
}

#' Update logging file and pass arguments to logger::log_trace
#'
#' The logging file to be used must be specified by the user in the \code{logging_file} field of \code{settings.ini}
#' @seealso logger::log_trace
#' @author pstahlhofen
#' @export
log_trace <- function(...) {
	update_appender()
	logger::log_trace(...)
}

#' Update logging file and pass arguments to logger::log_debug
#'
#' The logging file to be used must be specified by the user in the \code{logging_file} field of \code{settings.ini}
#' @seealso logger::log_debug
#' @author pstahlhofen
#' @export
log_debug <- function(...) {
	update_appender()
	logger::log_debug(...)
}

#' Update logging file and pass arguments to logger::log_warn
#'
#' The logging file to be used must be specified by the user in the \code{logging_file} field of \code{settings.ini}
#' @seealso logger::log_warn
#' @author pstahlhofen
#' @export
log_warn <- function(...) {
	update_appender()
	logger::log_warn(...)
}

#' Update logging file and pass arguments to logger::log_success
#'
#' The logging file to be used must be specified by the user in the \code{logging_file} field of \code{settings.ini}
#' @seealso logger::log_success
#' @author pstahlhofen
#' @export
log_success <- function(...) {
	update_appender()
	logger::log_success(...)
}

#' Update logging file and pass arguments to logger::log_error
#'
#' The logging file to be used must be specified by the user in the \code{logging_file} field of \code{settings.ini}
#' @seealso logger::log_error
#' @author pstahlhofen
#' @export
log_error <- function(...) {
	update_appender()
	logger::log_error(...)
}

#' Update logging file and pass arguments to logger::log_fatal
#'
#' The logging file to be used must be specified by the user in the \code{logging_file} field of \code{settings.ini}
#' @seealso logger::log_fatal
#' @author pstahlhofen
#' @export
log_fatal <- function(...) {
	update_appender()
	logger::log_fatal(...)
}
