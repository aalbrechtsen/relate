 
getPed <- function(pedfile){
  cat("Function argument should NOT be a filename \nname but rather a list given by the \"read.snps.pedfile\"\n  supplied by the snpMatrix package\n")
  #do tilde expansion
  pedfullname <- file.path(dirname(pedfile),basename(pedfile))
  return(.Call("getPedfile",pedfullfile,PACKAGE="Relate"))
}
