library("EnsDb.Hsapiens.v86")
library(BSgenome.Hsapiens.NCBI.GRCh38)
db <- EnsDb.Hsapiens.v86
bsg <- BSgenome.Hsapiens.NCBI.GRCh38

geneName <- "BCL6"

bcl6exons <- ensembldb::exonsBy(
  db, filter = AnnotationFilter::GenenameFilter(geneName))

devtools::use_data(bcl6exons)
