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

nog2goEnrichment <- function(nogList,subset,method="fisher"){
    if(!method%in%c("fisher","chisq")){
        stop(paste(method,"is not a recognised method, choose from \"fisher\" or \"chisq\"."))
    }

    goList <- nog2go(nogList)
    names(goList) <- names(nogList)
    goTerms <- unique(unlist(goList))

    missing = sum(!subset%in%names(goList))
    if(missing>0){
        warning(paste(missing,"genes had NOG terms with no matching GO terms and have been discounted."))
    }
    subset <- subset[subset%in%names(goList)]

    mat <- matrix(0,nrow=length(goList),ncol=length(goTerms))
    rownames(mat) <- names(goList)
    colnames(mat) <- goTerms
    mat[cbind(rep(names(goList),lengths(goList)),unlist(goList))] <- 1

    bg <- apply(mat,2,function(x) sum(x))
    bgf <- bg/length(goList)
    fg <- apply(mat[subset,],2,function(x) sum(x))
    fgf <- fg/length(subset)

    if(method=="fisher"){
        pvalues <- sapply(1:length(bg),function(x) fisher.test(matrix(c(bg[x],nrow(mat)-bg[x],fg[x],length(subset)-fg[x]),nrow=2))$p.value)
    }
    if(method=="chisq"){
        pvalues <- sapply(1:length(bg),function(x) chisq.test(matrix(c(bg[x],nrow(mat)-bg[x],fg[x],length(subset)-fg[x]),nrow=2))$p.value)
    }
    directions <- c("Over-represented","Under-represented")[1+as.numeric(fgf<bgf)]

    results = data.frame(Pvalue=pvalues,Direction=directions,fgFreq=fgf,bgFreq=bgf)
    results <- results[order(results$Pvalue),]

    return(results)
}
    



