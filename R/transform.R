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

primersToPairs <- function(primers) {
  lefts <- extractCoords(primers$PRIMER_LEFT)
  rights <- extractCoords(primers$PRIMER_RIGHT)
  lapply(seq_along(lefts), function(i) {
    cbind(rbind(lefts[i,], rights[i,]), 
    direction =c("forward", "reverse"))
  })
}

pairToGenome <- function(pair, seqs) {
  res <- lapply(
    list(start = IRanges::start(pair), end = IRanges::end(pair)),
    .toGenome,
    granges = seqs
  )
  res <- do.call(rbind, Map(range, res$start, res$end))
  res <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(
      unique(GenomicRanges::seqnames(seqs[1])), lengths = 2),
    IRanges::IRanges(start = res[, 1], end = res[, 2]),
    strand = S4Vectors::Rle(unique(GenomicRanges::strand(seqs[1])),2)
  )
  res
}