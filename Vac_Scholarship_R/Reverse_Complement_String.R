library(readr)
DNA <- read_file("rosalind_revc.txt")
DNA

reverse_chars <- function(string) {
  string_split = strsplit(as.character(string), split = "")
  reversed_split = string_split[[1]][nchar(string):1]
  paste(reversed_split, collapse="")
}
DNAREV <- reverse_chars(DNA)
DNAREV
chartr("ATCG", "TAGC", DNAREV)
       