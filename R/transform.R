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

primerCoords <- function(coords) {
  coords <- do.call(
    rbind,
    lapply(coords,
           function(x) {
             x <- strsplit(x, ",")[[1]]
             as.numeric(x)
           }))
  colnames(coords) <- c("start", "width")
  as.data.frame(coords)
}

primersToRanges <- function(res) {
  lefts <- do.call(IRanges::IRanges, primerCoords(res$PRIMER_LEFT))
  rights <- do.call(IRanges::IRanges, primerCoords(res$PRIMER_RIGHT))
  lapply(seq_along(lefts), function(i) {
    c(lefts[i], rights[i])
  })
}

pairToGRanges <- function(pair, seqs) {
  seqstrand <- as.character(strand(seqs))[1]
  res <- lapply(
    list(start = IRanges::start(pair), end = IRanges::end(pair)),
    .toGenome,
    upExon = seqs[1],
    downExon = seqs[2],
    strand = seqstrand
  )
  res <- do.call(rbind, Map(range, res$start, res$end))
  res <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(
      unique(GenomicRanges::seqnames(seqs[1])), lengths = 2),
    IRanges::IRanges(start = res[, 1], end = res[, 2]),
    strand = S4Vectors::Rle(unique(GenomicRanges::strand(seqs[1])),2)
  )
}