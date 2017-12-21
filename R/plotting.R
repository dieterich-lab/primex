#' Plot exon structure 
#'
#' @param ex exon structure
#' @param colours colours
#'
#' @export
#' @import  ggplot2
plotSegments <- function(ex, colours = NULL) {
  segs <- gr2segments(ex) 
  segs$colours <- "black"
  if (!is.null(colours)) {
    segs$colours <- colours
  }
  q <- ggplot2::ggplot() + 
    ggplot2::geom_segment(data = segs,
                 mapping = ggplot2::aes(
                   x = start,
                   xend = end,
                   y = track_id,
                   yend = track_id,
                   colour = I(colours)),
                 size = 8) + 
    ggplot2::geom_vline(xintercept = unique(c(segs$start, segs$end)),
               colour = "grey80") + 
    ggplot2::theme_classic() + 
    ggplot2::theme(axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())+
    ggplot2::ylab("")  +
    ggplot2::xlab("")
  print(q)
}

gr2segments <- function(gr, normalise = TRUE) {
  if (class(gr) %in% c("GRangesList", "list")) {
    gr <- unlist(GenomicRanges::GRangesList(gr))
  }
  gr$track_id <- names(gr)
  gr <- unname(gr)
  gr <- as.data.frame(gr)
  if (normalise) {
    norm <- normaliseData(gr = gr)
  }
  gr[c("start", "end")] <- norm$gr
  gr 
}

unifyDiff <- function(x,y) {
  points <- cbind(x,y)
  o <- order(points)
  res <- rle(points[o]) 
  res$values <- seq_along(unique(res$values))
  points[o] <- inverse.rle(res)
  list(points[seq_along(x)],
       points[seq_along(y) + length(x)])
}

normaliseData <- function(...){
  dat <- list(...)
  toProcess <- vapply(dat, Negate(is.null), logical(1))
  columns <- c("start", "end")
  positions <- do.call(rbind,
                       lapply(names(dat)[toProcess], function(x){
                         cbind(id = x, dat[[x]][, columns])
                       })
  )
  result <- as.data.frame(
    do.call(cbind,
            unifyDiff(positions$start, positions$end))
  )
  names(result) <- columns
  dat[toProcess] <- split(result, positions$id)
  dat
}


  
