#### Other utility functions #####

#' Return a timestamp for logging
#'
#' @returns the timestamp
#' @export
#'
timestamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

#' Translate ungapped coordinates back to gapped coordinates
#'
#' @param site.no.gap  the integer site in an ungapped sequence to convert
#' @param gapped.seq the sequence with gaps from an alignment
#'
#' @returns the cooordinate in the gapped sequence
#' @export
#'
convert.to.gapped.coordinate <- function(site.no.gap, gapped.seq) {
  if (is.na(site.no.gap)) {
    return(NA)
  }
  gapped.seq.char <- as.character(gapped.seq)
  # find gaps and stop codons
  gaps <- stringr::str_locate_all(gapped.seq.char, "-|\\*")[[1]][, 1]
  n <- site.no.gap
  for (i in gaps) {
    if (i <= n) n <- n + 1
  }
  n
}

#' Convert FASTA format to clustal style format.
#'
#' Prints and returns as a text object
#'
#' @param alignment seqinr alignment format
#' @param names optional character vector of names (if null, alignment names are used)
#' @param names.length override the default width of the names column (if na, default is used)
#' @param chunksize number of letters per row
#'
#' @returns the formatted text
#' @export
#'
printMultipleAlignment <- function(alignment, names = NULL, names.length = NA, chunksize = 60) {
  # this function requires the Biostrings package
  # find the number of sequences in the alignment
  numseqs <- alignment$nb

  if (is.null(names)) {
    names <- alignment$nam
  }

  if (is.na(names.length)) {
    names.length <- max(nchar(names))
  }

  # find the length of the alignment
  alignmentlen <- nchar(alignment$seq[[1]])

  # Calculate the start position of each line
  line.starts <- seq(1, alignmentlen, by = chunksize)

  # How many blocks are needed
  n.blocks <- length(line.starts)

  # get the alignment for each  sequences
  aln <- unlist(alignment$seq)
  lettersprinted <- rep(0, numseqs)

  create.block <- function(start) {
    block.lines <- rep("", numseqs + 1)
    block.lines[numseqs + 1] <- "\n"
    for (j in 1:numseqs) {
      alnj <- aln[j]
      chunkseq <- toupper(substring(alnj, start, start + chunksize - 1))

      # Calculate how many residues of the sequence we have printed so far in the alignment
      # Total minus gaps
      lettersprinted[j] <<- lettersprinted[j] + chunksize - Biostrings::countPattern("-", chunkseq)
      block.lines[j] <- paste0(sprintf(paste0("%", names.length, "s"), names[j]), "\t", chunkseq, " ", lettersprinted[j])
    }

    paste0(block.lines, "\n")
  }

  result <- unlist(lapply(line.starts, create.block))
  cat(paste0(result, "\n"))
  result
}

#' Extract sequence from a Biostrings alignment.
#'
#' Given a Biostrings alignment, extract the sequence string at the given coordinates
#'
#' @param aln msa in Biostrings format
#' @param sequence.name the name of the sequence from the alignment to extract
#' @param start,end coordinates in the gapped alignment
#'
#' @returns the sequence from the alignment
#' @export
subset.sequence <- function(aln, sequence.name, start, end) {
  as.character(aln@unmasked[[sequence.name]][start:end])
}

#' Reroot the given tree to the given tip labels
#'
#'  If the tip labels contains more than one node, their MRCA will be used as the root.
#'
#' @param tree the tree to reroot (e.g. ape)
#' @param node.labels  vector of tips to root on
#' @param position the position to place the root on the selected branch
#'
#' @returns the rerooted tree
#' @export
#'
reroot.tree <- function(tree, node.labels, position = 0.01) {
  if (!all(node.labels %in% tree$tip.label)) {
    cat("Not all nodes (", node.labels, ") are in the tree, cannot reroot, not changing tree\n")
    return(tree)
  }
  if (length(node.labels) == 1) {
    return(phytools::reroot(tree, which(tree$tip.label == node.labels), position = position))
  }
  root.node <- ape::getMRCA(tree, node.labels)
  return(phytools::reroot(tree, root.node, position = position))
}

#' Group a tree by a column from a data.frame.
#'
#' The metadata data.frame provided should have a 'value' column with the tip labels
#' in the tree. The contents of the column with the name in 'grouping.column' will
#' be used to create a tree group. The tree group will have the same name as the
#' grouping column.
#'
#' @param tree the tree to group
#' @param metadata a data.frame with columns 'value' and grouping.column
#' @param grouping.column the name of the column in metadata with the grouping variables
#'
#' @returns the tree with the added grouping
#' @export
add.tree.grouping <- function(tree, metadata, tip.label.column = "tip.label", grouping.column) {
  clade.data <- unlist(sapply(tree$tip.label, \(x) metadata[metadata[, tip.label.column] == x, grouping.column]))
  clade.groups <- split(tree$tip.label, clade.data)
  tidytree::groupOTU(tree, clade.groups, group_name = grouping.column)
  tree
}


#' Change the tip labels of a tree to a column from a data.frame.
#'
#' The metadata data.frame provided should have a column with the new tip labels
#' in the tree. The contents of the column with the name in 'grouping.column' will
#' be used to create a tree group. The tree group will have the same name as the
#' grouping column.
#'
#' @param tree the tree to group
#' @param metadata a data.frame with columns 'value' and grouping.column
#' @param original.tip.label.column the name of the column in metadata with the original tip labels
#' @param new.tip.label.column the name of the column in metadata with the new tip labels
#'
#' @returns the tree with the added grouping
#' @export
update.tip.labels <- function(tree, metadata,
                              original.tip.label.column = "value",
                              new.tip.label.column = "tip.label") {
  attr(tree, "old.tip.labels") <- tree$tip.label
  new.tip.labels <- unlist(sapply(tree$tip.label, \(x) metadata[metadata[, original.tip.label.column] == x, new.tip.label.column]))
  tree$tip.labels <- new.tip.labels
  tree
}


#' Convert a Biostrings MSA to a character list
#'
#' Convert a Biostrings MSA to a character list format that can be exported as FASTA
#'
#' @param aln the Biostrings MSA
#'
#' @returns a character list
#' @export
#'
biostrings.aln.to.list <- function(aln) {
  as.list(apply(as.matrix(aln), 1, paste, collapse = ""))
}

#' Convert a matrix MSA to a character list format
#'
#' @param aln the MSA as a matrix
#'
#' @returns a character list
#' @export
#'
matrix.aln.to.list <- function(aln) {
  as.list(apply(aln, 1, paste, collapse = ""))
}
