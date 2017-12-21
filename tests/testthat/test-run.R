context("Running Primer3")

test_that("can execute", {
  seqOpts <-setSeqOptions()
  seqOpts$SEQUENCE_ID <- "bcl6"
  seq <- paste0("TGTGAGAAGTGTAACCTGCATTTCCGTCACAAAAGCCAGCTGC",
           "GACTTCACTTGCGCCAGAAGCATGGCGCCATCACCAACACCAAGGTGCAATACCGCGTGTCA")
  seqOpts$SEQUENCE_TEMPLATE <- seq
  #seqOpts$SEQUENCE_OVERLAP_JUNCTION_LIST <- nchar(sjexons[[1]])
  runPrimer3(seqOpts)
})
