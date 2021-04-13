#' Print as heading
#'
#' @description
#' heading prints an input string in a "heading" format
#'
#' @details
#' A function primarily for aesthetics within the PathWAS pipeline and for progress checks throughout.
#'
#' @param sentence character. A string of variable length for printing.
heading = function(sentence) {

  cat(paste0("\n\n\n=======================\n\n",
                   sentence, "\n===========================\n\n"))
}
