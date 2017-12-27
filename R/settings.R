#' Pipe operator
#'
#' @name %>%
#' @noRd
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
NULL


#' Set primer size
#'
#' @param settings a settings list created by `p3Settings`
#' @param min minimal length (PRIMER_MIN_SIZE)
#' @param optimal optimal length (PRIMER_OPT_SIZE)
#' @param max maximal length (PRIMER_MAX_SIZE)
#'
#' @return a settings list with updated values
#' @importFrom magrittr "%>%"
#' @export
#'
#' @examples
#' p3Settings() %>%
#'   primerSize(min = 18, optimal = 20, max = 18)
primerSize <- function(settings, min = NULL, optimal = NULL, max = NULL) {
    settings$PRIMER_MIN_SIZE <- list(min)
    settings$PRIMER_OPT_SIZE <- list(optimal)
    settings$PRIMER_MAX_SIZE <- list(max)
    settings
}
