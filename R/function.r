#' Function to return the GO category codes for a given NOG term.
#'
#' @param nog A character vector with the NOG terms to look up.
#' @details
#' If there are no codes for a given term, the function will return character(0).
#' @keywords None
#' @return Either a vector of GO category codes, or a list of vectors if multiple NOG terms were given.
#' @export
#' @author Chris Field <fieldc@@ethz.ch>
#' @examples
#' None

nog2go <- function(nog){
    gos <- sapply(nog,function(x) nogGoMap[[x]])
    if(length(nog)<2){
        gos <- as.vector(gos)
    }
    return(gos)
}
