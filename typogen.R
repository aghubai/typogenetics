
# Typogenetics, a game by Douglas R. Hofstadter
#   see Hofstadter, Douglas R. "Goedel, Escher, Bach: an eternal golden braid." Harmondsworth: Penguin, 1979.

# code by Andras G. Hubai

# ---

# turn strand into enzymes
translate <- function (strand, nts, operations) {
  
  # turn strand into codons (ie. indices thereof)
  strnum <- c(factor(unlist(strsplit(strand, "")), nts))
  length(strnum) <- length(strnum) %/% 2 * 2
  strmat <- matrix(strnum, 2)
  codind <- (strmat[1, ] -1) * 4 + strmat[2, ]
  
  # mind the punctuation mark
  whstop <- which(codind == 1)
  lens <- diff(c(0, whstop, length(codind)+1)) -1
  lsta <- c(1, whstop +1)
  lend <- c(whstop -1, length(codind))
  whpos <- which(lens > 0)
  enzpcs <- lapply(whpos, function (i) lsta[i]:lend[i])
  enzn <- length(whpos)
  
  # turn codons into enzymes
  ops <- operations[2, codind]
  enzops <- lapply(enzpcs, function (i) ops[i])
  
  # find binding preference
  kns <- c(factor(operations[3, codind], c("l", "s", "r"))) - 2
  enzkns <- sapply(enzpcs, function (i) sum(kns[i[-c(1, length(i))]]) %% 4 +1)
  if (enzn > 0)
    enzpref <- c("A", "G", "T", "C")[enzkns] else character(0)
  
  return(list(enzn, enzops, enzpref))}

# shape strand using enzyme (determ: bind leftmost target)
operate <- function (strand, enzyme, determ=1, show=FALSE) {
  
  # bind enzyme to strand
  pref <- enzyme[[2]]
  stra <- c(" ", unlist(strsplit(strand, "")), " ")
  whpref <- which(stra == pref)
  if (determ > 0)
    pos <- whpref[determ] else
      pos <- whpref[sample(length(whpref), 1)]
  
  # manipulate (tracking substrate, position and copy mode)
  coms <- paste0("f", enzyme[[1]])
  prod <- list(rbind(" ", stra), pos, 0)
  if (show) printst(prod, "start")
  i <- 1
  while (i <= length(coms) && prod[[2]] > 0) {
    prod <- get(coms[i])(prod)
    if (show) printst(prod, enzyme[[1]][i])
    i <- i +1}
  
  # separate strands
  separ <- c(rev(prod[[1]][1, ]), " ", prod[[1]][2, ])
  stranew <- unlist(strsplit(paste(separ, collapse=""), " "))
  stranew <- stranew[which(stranew != "")]
  
  return(stranew)}

fcut <- function (state) {
  # cuts to the right of position, creating separate molecules
  subs <- state[[1]]; pos <- state[[2]]
  state[[1]] <- matrix(c(subs[seq_len(pos*2)], rep(" ", 2), subs[-seq_len(pos*2)]), 2)
  return(state)}

fdel <- function (state) {
  # deletes only on active strand, moves to the right
  subs <- state[[1]]; pos <- state[[2]] +1
  subs[2, pos -1] <- " "
  # ends: if runs off molecule
  if (subs[2, pos] == " ") pos <- 0 else
  # copies: if there is no complementary base (copy mode on)
    if (subs[1, pos] == " " && state[[3]] == 1)
      subs[1, pos] <- compl(subs[2, pos])
  state[[1]] <- subs; state[[2]] <- pos
  return(state)}

fswi <- function (state) {
  # "rotates" strands and repositions
  subs <- state[[1]]; pos <- state[[2]]
  if (subs[1, pos] != " ") {
    state[[1]] <- matrix(rev(subs), 2)
    state[[2]] <- ncol(subs) - pos +1
  } else
    # detaches
    state[[2]] <- 0
  return(state)}

fmvr <- function (state) {
  subs <- state[[1]]; pos <- state[[2]] +1
  # ends: if runs off molecule
  if (subs[2, pos] == " ") pos <- 0 else
    # copies: if there is no complementary base (copy mode on)
    if (subs[1, pos] == " " && state[[3]] == 1)
      subs[1, pos] <- compl(subs[2, pos])
  state[[1]] <- subs; state[[2]] <- pos
  return(state)}

fmvl <- function (state) {
  subs <- state[[1]]; pos <- state[[2]] -1
  # ends: if runs off molecule
  if (subs[2, pos] == " ") pos <- 0 else
    # copies: if there is no complementary base (copy mode on)
    if (subs[1, pos] == " " && state[[3]] == 1)
      subs[1, pos] <- compl(subs[2, pos])
  state[[1]] <- subs; state[[2]] <- pos
  return(state)}

fcop <- function (state) {
  # turns on copy mode, and copies at position
  subs <- state[[1]]; pos <- state[[2]]
  if (subs[1, pos] == " ")
    state[[1]][1, pos] <- compl(subs[2, pos])
  state[[3]] <- 1
  return(state)}

foff <- function (state) {
  # turns off copy mode
  state[[3]] <- 0
  return(state)}

# insert pertains to both strands if copy mode is on
fina <- function (state) {
  subs <- state[[1]]; pos <- state[[2]]
  if (state[[3]] == 1) add <- c("T", "A") else add <- c(" ", "A")
  state[[1]] <- matrix(c(subs[seq_len(pos*2)], add, subs[-(seq_len(pos*2))]), 2)
  state[[2]] <- pos +1
  return(state)}

finc <- function (state) {
  subs <- state[[1]]; pos <- state[[2]]
  if (state[[3]] == 1) add <- c("G", "C") else add <- c(" ", "C")
  state[[1]] <- matrix(c(subs[seq_len(pos*2)], add, subs[-(seq_len(pos*2))]), 2)
  state[[2]] <- pos +1
  return(state)}

fing <- function (state) {
  subs <- state[[1]]; pos <- state[[2]]
  if (state[[3]] == 1) add <- c("C", "G") else add <- c(" ", "G")
  state[[1]] <- matrix(c(subs[seq_len(pos*2)], add, subs[-(seq_len(pos*2))]), 2)
  state[[2]] <- pos +1
  return(state)}

fint <- function (state) {
  subs <- state[[1]]; pos <- state[[2]]
  if (state[[3]] == 1) add <- c("A", "T") else add <- c(" ", "T")
  state[[1]] <- matrix(c(subs[seq_len(pos*2)], add, subs[-(seq_len(pos*2))]), 2)
  state[[2]] <- pos +1
  return(state)}

# search along the strand
frpy <- function (state) {
  subs <- state[[1]]; pos <- state[[2]] +1
  while (pos != 0) {
    # ends: if runs off molecule
    if (subs[2, pos] == " ") pos <- 0 else {
      # copies: if there is no complementary base (copy mode on)
      if (subs[1, pos] == " " && state[[3]] == 1)
        subs[1, pos] <- compl(subs[2, pos])
      # moves: if target is not found
      if (subs[2, pos] %in% c("C", "T")) break() else
        pos <- pos +1}}
  state[[1]] <- subs; state[[2]] <- pos
  return(state)}

frpu <- function (state) {
  subs <- state[[1]]; pos <- state[[2]] +1
  while (pos != 0) {
    # ends: if runs off molecule
    if (subs[2, pos] == " ") pos <- 0 else {
      # copies: if there is no complementary base (copy mode on)
      if (subs[1, pos] == " " && state[[3]] == 1)
        subs[1, pos] <- compl(subs[2, pos])
      # moves: if target is not found
      if (subs[2, pos] %in% c("A", "G")) break() else
        pos <- pos +1}}
  state[[1]] <- subs; state[[2]] <- pos
  return(state)}

flpy <- function (state) {
  subs <- state[[1]]; pos <- state[[2]] -1
  while (pos != 0) {
    # ends: if runs off molecule
    if (subs[2, pos] == " ") pos <- 0 else {
      # copies: if there is no complementary base (copy mode on)
      if (subs[1, pos] == " " && state[[3]] == 1)
        subs[1, pos] <- compl(subs[2, pos])
      # moves: if target is not found
      if (subs[2, pos] %in% c("C", "T")) break() else
        pos <- pos +1}}
  state[[1]] <- subs; state[[2]] <- pos
  return(state)}

flpu <- function (state) {
  subs <- state[[1]]; pos <- state[[2]] -1
  while (pos != 0) {
    # ends: if runs off molecule
    if (subs[2, pos] == " ") pos <- 0 else {
      # copies: if there is no complementary base (copy mode on)
      if (subs[1, pos] == " " && state[[3]] == 1)
        subs[1, pos] <- compl(subs[2, pos])
      # moves: if target is not found
      if (subs[2, pos] %in% c("A", "G")) break() else
        pos <- pos +1}}
  state[[1]] <- subs; state[[2]] <- pos
  return(state)}

# find complementer bases
compl <- function (nt) {
  b <- factor(nt, c("A", "C", "G", "T"))
  levels(b) <- c("T", "G", "C", "A")
  return(as.character(b))}

# find complementer strand
complst <- function (strand) {
  b <- compl(unlist(strsplit(strand, "")))
  return(paste0(b, collapse=""))}

# separate enzymes (change format: enzymes -> enzlist)
separenz <- function (enzymes) {
  n <- enzymes[[1]]
  enzlist <- vector("list", n)
  for (i in seq_len(n))
    enzlist[[i]] <- list(enzymes[[2]][[i]], enzymes[[3]][i])
  return(enzlist)}

# number of possible starting positions for enzyme along strand
numstart <- function (strand, enzyme) {
  return(sum(unlist(strsplit(strand, "")) == enzyme[[2]]))}

# brief output of strands
printst <- function (state, comm=NULL) {
  if (!is.null(comm))
    cat(paste0(c(comm, ":\n"), collapse=""))
  cat(paste0(state[[1]][1, ], collapse=""), "\n")
  cat(paste0(state[[1]][2, ], collapse=""), "\n")
  symb <- ifelse(state[[3]] == 1, "|", "-")
  if (state[[2]] != 0)
    cat(paste0(c(rep(" ", state[[2]] -1), symb), collapse=""), "\n") else
      cat("*\n")}

# transform code_sec.txt and code_thr.txt into code.txt
codetrans <- function (f="code.txt", f2="code_sec.txt", f3="code_thr.txt") {
  ord <- read.table(f2, stringsAsFactors=FALSE)
  knk <- read.table(f3, stringsAsFactors=FALSE)
  
  nuc <- paste(ord[1, -1], collapse=" ")
  cod <- c(sapply(2:nrow(ord), function (i)
    sapply(2:ncol(ord), function (j)
      paste0(as.character(c(ord[i, 1], ord[1, j], " ", ord[i, j],
                            " ", knk[i, j])), collapse=""))))
  
  write(c(nuc, cod), file=f)}
