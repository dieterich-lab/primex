.toGenome <- function(x, granges) {
  granges <- as.data.frame(granges)
  up   <- granges[1,, drop = FALSE]
  down <- granges[2,, drop = FALSE]
  strand <- unique(granges$strand)
  inUpstream <- x <= up$width
  if (strand == '+') {
    x[ inUpstream] <- x[ inUpstream] + up$start - 1
    x[!inUpstream] <- x[!inUpstream] - up$width + down$start - 1
  } else {
    x[ inUpstream] <- -x[ inUpstream] + up$end + 1
    x[!inUpstream] <- -x[!inUpstream] + down$end + up$width + 1
  }
  x
}

extractCoords <- function(coords) {
  coords <- strsplit(coords, ",") %>% lapply(as.integer)
  stopifnot(all(lapply(coords, length) == 2))
  do.call(rbind, coords) %>% 
    as.data.frame() %>% 
    setNames(c("start", "width"))
}

#' Return primer coordinates.
#' 
#' In Primer3 LEFT and RIGHT are 0-based.
#' RIGHT start is 5'-end of the primer, we coerce it to IRanges style.
#'
#' @param primers 
#' @noRd
primersToPairs <- function(primers) {
  lefts <- extractCoords(primers$PRIMER_LEFT)
  lefts$start <- lefts$start + 1
  rights <- extractCoords(primers$PRIMER_RIGHT)
  rights$start <- rights$start - rights$width + 2
  lapply(seq_along(lefts), function(i) {
    cbind(rbind(lefts[i,], rights[i,]), 
          direction = c("forward", "reverse"))
  })
}

primersToGRanges <- function(primers, exons) {
  pairs <- primersToPairs(primers) %>%
   lapply(., splitPrimer,IRanges::width(exons))
  pairs
}


splitPrimer <- function(x, upExonWidth) {
  res <- lapply(1:2, function(i) {
    y <- x[i,]
    if (y$start <= upExonWidth && y$start + y$width - 1 > upExonWidth) {
      y <- rbind(y, y, make.row.names = FALSE)
      y$width[1] <- upExonWidth - y$start[1] + 1
      y$start[2] <- upExonWidth + 1
      y$width[2] <- y$width[2] - y$width[1]
    }
    y
  })
  do.call(rbind, res)
}


pairToGenome <- function(pair, seqs) {
  pair$end   <- .toGenome(pair$start + pair$width, granges = seqs)
  pair$start <- .toGenome(pair$start, granges = seqs)
  res <- do.call(rbind, Map(range, res$start, res$end))
  res <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(
      unique(GenomicRanges::seqnames(seqs[1])), lengths = 2),
    IRanges::IRanges(start = res[, 1], end = res[, 2]),
    strand = S4Vectors::Rle(unique(GenomicRanges::strand(seqs[1])),2)
  )
  res
}