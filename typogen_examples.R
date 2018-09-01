
# Typogenetics, a game by Douglas R. Hofstadter
#   see Hofstadter, Douglas R. "Goedel, Escher, Bach: an eternal golden braid." Harmondsworth: Penguin, 1979.

# code by Andras G. Hubai

# ---

rm(list=ls())

source("typogen.R")
#if (!file.exists("code.txt")) codetrans()

nts <- scan("code.txt", character(), nlines=1, quiet=TRUE)
codons <- c(sapply(nts, function (i) paste0(i, nts)))
operations <- matrix(scan("code.txt", character(), skip=1, quiet=TRUE), 3)
if (!identical(codons, operations[1, ])) stop("corrupted code.txt")

# 1st case
strand <- "ACA"
enzyme <- list(c("del", "mvr", "int"), "A")
operate(strand, enzyme, 1, TRUE)

operate(strand, enzyme, 2, TRUE)

# 2nd case
strand <- "CAAAGAGAATCCTCTTTGAT"
enzyme <- list(c("rpy", "cop", "rpu", "cut"), "A")
operate(strand, enzyme, 2, TRUE)

# 3rd case
strand <- "TAGATCCAGTCCATCGA"
enzyme <- list(c("rpu", "inc", "cop", "mvr", "mvl", "swi", "lpu", "int"), "G")
operate(strand, enzyme, 2, TRUE)

# 4th case
strand <- "TAGATCCAGTCCACATCGA"
translate(strand, nts, operations)

# 5th case
strand <- "CGGATACTAAACCGA"
enzymes <- translate(strand, nts, operations)
separenz(enzymes)

# 6th case
set.seed(1)
strand <- paste0(sample(nts, 1 + rbinom(1, 38, 0.5), TRUE), collapse="")
enzymes <- translate(strand, nts, operations)

for (i in separenz(enzymes)) {
  numst <- numstart(strand, i)
  for (j in seq_len(numst))
    print(operate(strand, i, j))
}

# 7th case
strand <- "CGTCTCTCTATAGAGAGACG"
enzymes <- translate(strand, nts, operations)
(enzyme <- separenz(enzymes)[[1]])
(stranew <- operate(strand, enzyme, 1))
strand == stranew
