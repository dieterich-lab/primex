
#' Get exon pairs
#'
#' @param exonsByTx a list of exons returned by ensembl exonsBy tx
#'
#' @return a list of character vectors
#' @export
#'
extractExonPairs <- function(exonsByTx) {
  assertColumns(exonsByTx, "exon_id")
  lapply(exonsByTx, splitToPairs)
}

# works for mcols and data.frames
assertColumns <- function(x, ...) {
  columns <- as.list(...)
  ismissing <- vapply(columns, function(name) is.null(`$`(x, name)), logical(1))
  if (any(ismissing))
    stop(paste("Missing columns: ", paste(name, collapse = ", ")))
}

# create a list of character(2) using the provided order
splitToPairs <- function(exonIds) {
  lapply(seq_len(exonIds), function(i) {
    c(exonIds[i - 1], exonIds[i])
  })
}

#' Given a list of exons returns candidates for primer design
#' 
#' The transcripts of single exon are excluded!
#'
#' @param exonsByTx exon GRanges by transcripts
#' @param tolerance the minimal tolerable difference in splice junctions 
#'   coordinates
#'
#' @return a GRangesList
#' @export
#'
exonsBySJ <- function(exonsByTx, tolerance = 0) {
  # filter single exon tx's
  exonsByTx <- exonsByTx[vapply(exonsByTx, length, integer(1)) > 1]
  expairsbyTx <- extractExonPairs(exonsByTx)
  # transform to pairs of names 
  allExons <- unlist(exonsByTx)
  # transform pairs of exon_id --> GRanges
  expairsbyTx <- lapply(expairsbyTx, function(exonsInTx) {
    res <- lapply(exonsInTx, function(p) {
      pair <- allExons[match(p, allExons$exon_id)]
      names(pair) <- pair$exon_id
      pair
    })
    names(res) <- lapply(exonsInTx, paste, collapse="|")
    res
  })
  expairsbyTx <- filterByDistance(expairsbyTx, tolerance)
  expairsbyTx
}

# returns a filtered list
filterByDistance <- function(expairs, tolerance) {
  sjCoords <- lapply(unlist(expairs), pair2sj)
  sjCoords <- do.call(rbind, unname(sjCoords))
  d <- as.matrix(stats::dist(sjCoords[c("start", "end")], method = "minkowski", p = 1))
  # check if second minimal distance (after self==0) is above tolerance
  isDuplicate <- tolerance >= apply(d, 1, function(x) sort(x)[2])
  toInclude <- names(isDuplicate)[!isDuplicate]
  # filter
  lapply(expairs, function(pairsInTx) {
    pairIds <- vapply(pairsInTx, function(x)
      paste(x$exon_id, collapse = "|"), character(1))
    pairsInTx[pairIds %in% toInclude]
  })
}

# take end of the first and the start of the second exon
# assume that they are ordered by the exon_rank
pair2sj <- function(p) {
  x <- as.data.frame(p)
  res <- x[1,]#c("seqnames", "start", "end", "strand")]
  res[1,c("start", "end")] <- range(x$end[1], x$start[2])
  rownames(res) <- paste0(p$exon_id, collapse = "|")
  res
}


#' Return sequences for the exon pairs
#'
#' @param exonPairs a list of GRanges with the splice junction exon pairs.
#'   For a single pair a single GRanges object may be provided.
#' @param src whether a BSgenome instance, e.g. BSgenome.Hsapiens.NCBI.GRCh38 or
#'   a FaFile object
#'
#' @return a list of GRanges with additional `seq` data column with
#'   the retrieved sequence.
#' @export
#'
exonSeqs <- function(exonPairs, src) {
  # extract all names 
  if (!is.list(exonPairs))
    exonPairs <- list(exonPairs)
  allExons <- exonPairs
  if (!inherits(allExons, "GRanges")) {
    stopifnot(all(vapply(
      allExons, inherits, "GRanges", FUN.VALUE = logical(1)
    )))
    allExons <- do.call(c, unname(allExons))
  }
  allExons <- allExons[!duplicated(allExons$exon_id)]
  if (inherits(src, "FaFile")) {
    exseqs <-  as.character(Rsamtools::getSeq(x = src, allExons))
  } else  if (inherits(src, "BSgenome")) {
    exseqs <- BSgenome::getSeq(bsg, allExons, as.character = TRUE)
  } else {
    stop("src must be FaFile or BSgenome object")
  }
  names(exseqs) <- allExons$exon_id
  # a list of transcripts with seqs of exons
  for (i in seq_along(exonPairs)) {
      ids <- exonPairs[[i]]$exon_id
      S4Vectors::mcols(exonPairs[[i]])$seq <-  exseqs[ids]
  }
  exonPairs
}
