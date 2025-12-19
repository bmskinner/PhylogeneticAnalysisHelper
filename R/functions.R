#' Test if the given binary name is on the PATH
#'
#' Returns true if the binary is on the PATH, false otherwise
#'
#' @param binary the binary to search for
#'
#' @returns true if the binary is on the PATH, false otherwise
#' @export
#'
#' @examples
#' is.on.PATH("macse")
is.on.PATH <- function(binary) {
  return(Sys.which(binary) != "")
}

#

#' If the given binary is not on the PATH, stop the script
#'
#' @param binary the binary to test
#'
#' @returns nothing; halts the script if the binary is not found
#' @export
#'
#' @examples
#' check.binary.exists("macse")
check.binary.exists <- function(binary) {
  if (!is.on.PATH(binary)) stop(paste("Cannot find", binary, "on PATH"))
}

#' Download ortholgoues for the given gene symbol using NCBI datasets
#'
#' Download ortholgoues for the given gene symbol using NCBI datasets. Extracts
#' the downloaded file, and moves sequence files to folder containing the zip
#' file. Expects 'datasets' on the PATH
#'
#' @param symbol the gene symbol to search
#' @param ortholog which orthologues to
#'   retrieve
#' @param include the sequence types to return. Defaults to cds
#' @param overwrite should we overwrite existing downloads. If false, and the outputs  exist, the function will return the output zip file name. Defaults to false.
#' @param filename the path for the downloaded zip file
#' @param ... other arguments to datasets
#'
#' @returns  a vector of FASTA file paths
#' @export
#' @importFrom rlang .data
#'
run.ncbi.datasets <- function(symbol, ortholog = "all", include = "cds", overwrite = FALSE,
                              filename, ...) {
  check.binary.exists("datasets")

  # Define where files will be moved after download
  data.dir <- dirname(filename)
  md5.file <- paste0(data.dir, "/", symbol, ".md5sum.txt")
  metadata.file <- paste0(data.dir, "/", symbol, ".metadata.jsonl")

  # Map the default file paths within the downloaded data to their extracted locations
  include.files <- data.frame(
    input = c(
      paste0("ncbi_dataset/data/", include, ".fna"),
      paste0("ncbi_dataset/data/data_report.jsonl"),
      paste0("md5sum.txt")
    ),
    output = c(
      paste0(data.dir, "/", symbol, ".", include, ".fna"),
      metadata.file,
      md5.file
    )
  )

  # Download only if needed
  if (!file.exists(filename) | overwrite) {
    system2("datasets", paste(
      "download gene symbol",
      symbol,
      "--ortholog", ortholog,
      "--include", include,
      "--filename", filename,
      ...
    ), # any other options
    stdout = paste0(filename, ".download.log"), # logs
    stderr = paste0(filename, ".download.log")
    ) # error logs

    # Extract the data files we need
    utils::unzip(filename, exdir = data.dir)

    for (i in 1:nrow(include.files)) {
      # Add folder string - kept out for name matching in md5sum
      in.file <- paste0(data.dir, "/", include.files[i, "input"])
      fs::file_copy(in.file,
        include.files[i, "output"],
        overwrite = TRUE
      )
    }
  }

  # Check the md5 matches
  check.md5 <- function(file, expected.md5) {
    observed.md5 <- tools::md5sum(file)
    expected.md5 == observed.md5
  }

  # Read the MD5 data and merge with the updated file paths
  md5.data <- utils::read.table(md5.file, sep = " ", col.names = c("MD5", "Empty", "input")) |>
    merge(include.files, y = _, by = "input", all.x = TRUE) |>
    dplyr::select(-"Empty") |>
    stats::na.omit() |>
    dplyr::rowwise() |>
    dplyr::mutate(is.correct = check.md5(.data$output, .data$MD5))

  if (!all(md5.data$is.correct)) {
    print(md5.data)
    stop("Error in data transfer: at least one MD5 is incorrect")
  }

  # Read the JSONL metadata included in teh download
  # Extract only the key metadata we wanted
  json.text <- base::readLines(metadata.file, warn = FALSE, encoding = "UTF-8")
  metadata.elements <- lapply(json.text, jsonlite::fromJSON)

  extract.key.metadata <- function(metadata.item) {
    unnull <- function(x) ifelse(is.null(x), NA, x)

    c(
      chromosome = unnull(metadata.item$chromosomes),
      common.name = unnull(metadata.item$commonName),
      GeneID = unnull(metadata.item$geneId),
      tax.name = unnull(metadata.item$taxname)
    )
  }

  metadata.values <- do.call(rbind, lapply(metadata.elements, extract.key.metadata)) |>
    as.data.frame() |>
    dplyr::rowwise() |>
    dplyr::mutate(display.name = stringr::str_to_sentence(ifelse(is.na(.data$common.name), .data$tax.name, .data$common.name)))


  unlink(paste0(data.dir, "/ncbi_dataset"), recursive = TRUE) # remove the extracted dir

  # Return the paths to the downloaded files and the metadata
  list(
    fasta.files = paste0(data.dir, "/", symbol, ".", include, ".fna"),
    gene.metadata = metadata.values
  )
}

#' Trim non-homologous sequences from an unaligned nucleotide sequence
#'
#' Using MACSE `trimNonHomologousFragments`. Uses default parameters of
#' min_homology_to_keep_seq 0.2, min_trim_in 40 and min_trim_ext 20. Expects
#' 'macse' on the PATH
#'
#' @param fa.file the FASTA file to trim
#' @param out.dir the directory to save the trimmed files. Defaults to 'aln'.
#' @param ... other arguments to macse
#'
#' @returns a list with the trimmed file names for NT and the output stats
#' @export
#'
#' @examples
#' # trim.files <- run.macse.trim.non.homologous.fragments(cds.file)
run.macse.trim.non.homologous.fragments <- function(fa.file, out.dir = "aln", ...) {
  check.binary.exists("macse")
  filesstrings::create_dir(out.dir)
  out.name <- tools::file_path_sans_ext(basename(fa.file))

  nt.out <- file.path(out.dir, paste0(out.name, ".trim.fa"))
  stat.out <- file.path(out.dir, paste0(out.name, ".trim.stats.csv"))
  system2("macse", paste(
    " -prog trimNonHomologousFragments",
    "-seq", fa.file, # input
    "-min_homology_to_keep_seq 0.2", # internal stop codon
    "-min_trim_in 40", # minimum length to internally trim
    "-min_trim_ext 20", # minimum length to externally trim
    "-out_trim_info", stat.out, # generate stat file per site
    "-out_NT", nt.out, # output nt alignment
    ...
  ), # any other options
  stdout = file.path(out.dir, paste0(out.name, ".macse.log")), # logs
  stderr = file.path(out.dir, paste0(out.name, ".macse.log"))
  ) # error logs
  # Return output file names
  list("nt" = nt.out, "stat" = stat.out)
}

#' Run MACSE alingment on an input FASTA file
#'
#' Outputs NT and AA alignment files. Expects 'macse' on the PATH
#'
#' @param fa.file the FASTA file to align
#' @param out.dir the directory to save the alignment files
#' @param ... other arguments to macse
#'
#' @returns a list with the alignment file names for AA and NT
#' @export
#'
run.macse.align.sequences <- function(fa.file, out.dir = "aln", ...) {
  check.binary.exists("macse")
  filesstrings::create_dir(out.dir)
  out.name <- tools::file_path_sans_ext(basename(fa.file))
  aa.out <- file.path(out.dir, paste0(out.name, ".aa.aln"))
  nt.out <- file.path(out.dir, paste0(out.name, ".nt.aln"))
  # Run a codon aware alignment with MACSE
  # Here, macse is an alias on the PATH to java -jar /path/to/macse/jar
  system2("macse", paste(
    " -prog alignSequences",
    "-seq", fa.file, # input
    "-out_NT", nt.out, # output nt alignment
    "-out_AA", aa.out, # output aa alignment
    ...
  ), # any other options
  stdout = file.path(out.dir, paste0(out.name, ".macse.log")), # logs
  stderr = file.path(out.dir, paste0(out.name, ".macse.log"))
  ) # error logs

  # Return output file names
  list("aa" = aa.out, "nt" = nt.out)
}

#' Remove stop codons from NT alignment.
#'
#' Using macse exportAlignment. Expects 'macse' on the PATH.
#'
#' @param nt.aln.file aligned FASTA file containing stop codons and or frameshifts.
#' @param out.dir the directory to save the alignment files
#' @param ... other arguments to macse
#'
#' @returns a list with the alignment file names for AA, NT and the stats file
#' @export
run.macse.export.alignment <- function(nt.aln.file, out.dir = "aln", ...) {
  check.binary.exists("macse")
  filesstrings::create_dir(out.dir)
  out.name <- tools::file_path_sans_ext(basename(nt.aln.file))
  nt.out <- file.path(out.dir, paste0(out.name, ".nt.clean.aln"))
  aa.out <- file.path(out.dir, paste0(out.name, ".aa.clean.aln"))
  stat.out <- file.path(out.dir, paste0(out.name, ".clean.stat.csv"))
  system2("macse", paste(
    " -prog exportAlignment",
    "-align", nt.aln.file, # input
    "-codonForFinalStop ---", # remove final stop codon
    "-codonForInternalStop NNN", # internal stop codon
    "-codonForInternalFS ---", # internal framshift to gap
    "-charForRemainingFS -", # other frameshift to gap

    "-out_stat_per_seq", stat.out, # generate stat file per site
    "-out_NT", nt.out, # output nt alignment
    "-out_AA", aa.out, # output aa alignment
    ...
  ), # any other options
  stdout = file.path(out.dir, paste0(out.name, ".macse.log")), # logs
  stderr = file.path(out.dir, paste0(out.name, ".macse.log"))
  ) # error logs

  # Return output file names
  list("aa" = aa.out, "nt" = nt.out, "stats" = stat.out)
}

#' Run MAFFT on an input FASTA
#'
#' Outputs data/seq.fa to aln/seq.aln
# Uses 6 threads by default. Expects 'mafft' on the PATH.
#'
#' @param fa.file the file of sequences to align
#' @param out.dir   the directory to save the alignment file
#' @param ... other arguments to mafft
#'
#' @returns the alignment file name
#' @export
run.mafft <- function(fa.file, out.dir = "aln", ...) {
  check.binary.exists("mafft")
  filesstrings::create_dir(out.dir)
  out.name <- tools::file_path_sans_ext(basename(fa.file))
  aln.file <- file.path(out.dir, paste0(out.name, ".aln"))
  cat("Executing: mafft", paste(...), fa.file, " > ", aln.file, "\n")
  system2("mafft", paste(
    ...,
    "--thread 6",
    fa.file, " > ", aln.file
  ))

  aln.file
}

#' Create a alignment from two existing alignments
#'
#' Create a alignment of two existing alignments (or an alignment and a single
#' sequence FASTA) This preserves the structure of the existing alignments.
#' Expects 'clustalo' on the PATH.
#'
#' @param aln.file.1 the first alignment file
#' @param aln.file.2 the second alignment file (or single sequence FASTA)
#' @param out.file the combined output alignment file
#'
#' @export
#'
run.clustalo.dual.profile <- function(aln.file.1, aln.file.2, out.file) {
  check.binary.exists("clustalo")
  system2("clustalo", paste(
    "--profile1", aln.file.1,
    "--profile2", aln.file.2,
    "-o", out.file,
    "--force"
  )) # overwrite existing alignments
}

#' Run muscle
#'
#' Expects 'muscle' on the PATH.
#'
#' @param fa.file the FASTA file to align
#' @param out.dir the directory to output the alignment. Defaults to 'aln'
#' @param ... other arguments to muscle
#'
#' @returns the path to the alignment file
#' @export
#'
run.muscle <- function(fa.file, out.dir = "aln", ...) {
  check.binary.exists("muscle")
  filesstrings::create_dir(out.dir)
  out.name <- tools::file_path_sans_ext(basename(fa.file))
  aln.file <- file.path(out.dir, paste0(out.name, ".aln"))
  system2("muscle", paste(
    "-align", fa.file,
    "-output", aln.file
  ))
  aln.file
}

#' Run IQTREE
#'
#' Expects 'iqtree3' on the PATH. By default, runs with 6 threads
#'
#' @param aln.file the alignment
#' @param ... other arguments to IQTREE
#'
#' @returns the path to the treefile
#' @export
#'
run.iqtree <- function(aln.file, ...) {
  check.binary.exists("iqtree3")
  out.dir <- dirname(aln.file)
  out.name <- tools::file_path_sans_ext(basename(aln.file))
  log.file <- file.path(out.dir, paste0(out.name, ".iqtree.log"))

  # Note - cluster default versions IQTREE 1.6 and 2.0.6 reliably segfault using -st CODON
  # Manually added 2.4.0 and 3.0.1 to PATH, appears stable
  system2("iqtree3", paste(
    "-s ", aln.file,
    "-nt 6", # number of threads
    ...
  ), # any other arguments to iqtree
  stdout = log.file,
  stderr = log.file
  )

  # Return file with the tree
  paste0(aln.file, ".treefile")
}

#' Run HyPhy MEME in a Unix environment
#'
#' Tests the given test and background branches in the tree file.
#' Assumes a conda environment named 'hyphy' is available
#'
#' @param nex.file the nexus file with the alignment
#' @param tree.file the tree for the alignment
#'
#' @returns the output file path
#' @export
run.hyphy.meme <- function(nex.file, tree.file) {
  json.file <- paste0(tree.file, ".json")
  bash.file <- paste0(tree.file, ".sh")
  # Create control script for MEME
  brio::write_file(
    paste0(
      "#!/bin/bash\n",
      "source activate hyphy\n\n",
      "# MEME test for positive selection at individual sites\n",
      "hyphy meme --alignment  ", nex.file, " --tree ", tree.file, " --branches Test --output ", json.file, "\n"
    ),
    bash.file
  )
  system2("bash", bash.file)
  json.file
}

#' Run divvier from within a conda environment
#'
#' Run divvier from within a conda environment (named divvier).Creates a divvied
#' alignment.
#' @param aln.file the alignment to divvy
#' @param ... other options to divvier
#'
#' @returns the divvied alignment file path
#' @export
#'
run.divvier <- function(aln.file, ...) {
  bash.file <- paste0(aln.file, ".sh")
  brio::write_file(
    paste0(
      "#!/bin/bash\n",
      "source activate divvier\n\n",
      "# Identify phylogenetically informatve sites with indels\n",
      "divvier -divvygap ", paste(...), " ", aln.file, "\n",
      "# Ensure file extension is aln for use in IQTREE\n",
      "mv ", paste0(aln.file, ".divvy.fas"), " ", paste0(aln.file, ".divvy.aln"), "\n"
    ),
    bash.file
  )
  system2("bash", bash.file)
  # Return the divvied file
  paste0(aln.file, ".divvy.aln")
}

#### Functions to plot or export data #####

#' Create an Excel file containing a data frame
#'
#' @param data the data to export
#' @param file.name the name of the file
#' @param cols.to.fixed.size.font integer vector of columns to set in monospace font
#'
#' @export
create.xlsx <- function(data, file.name, cols.to.fixed.size.font = NULL) {
  oldOpt <- options()
  options(xlsx.date.format = "yyyy-mm-dd") # change date format
  wb <- xlsx::createWorkbook(type = "xlsx")
  sh <- xlsx::createSheet(wb)
  xlsx::addDataFrame(data, sh, row.names = F)
  xlsx::createFreezePane(sh, 2, 2, 2, 2) # freeze top row and first column
  cs <- xlsx::CellStyle(wb) +
    xlsx::Font(wb, heightInPoints = 10, isBold = FALSE, name = "Courier New")
  if (!is.null(cols.to.fixed.size.font)) {
    for (i in xlsx::getCells(xlsx::getRows(sh), colIndex = cols.to.fixed.size.font)) {
      xlsx::setCellStyle(i, cs)
    }
  }

  xlsx::autoSizeColumn(sh, 1:ncol(data))
  xlsx::saveWorkbook(wb, file = file.name)
  options(oldOpt)
}

#' Save a plot to file
#'
#' Saves png and svg versions of the plot
#'
#' @param filename the file to save to
#' @param plot the plot to save. Defaults to the last generated plot.
#' @param width the width in mm
#' @param height the height in mm
#'
#' @export
#'
save.plot <- function(filename, plot = ggplot2::last_plot(), width, height) {
  ggplot2::ggsave(filename, plot, dpi = 600, units = "mm", width = width, height = height)
  ggplot2::ggsave(stringr::str_replace(filename, ".png$", ".svg"), plot,
    dpi = 300,
    units = "mm", width = width, height = height
  )
}

#' Save a double column plot (170mm) at 600dpi
#'
#' @param filename the file to save to
#' @param plot the plot to save. Defaults to the last generated plot.
#' @param height the height in mm. Defaults to 170mm (square plot).
#'
#' @export
#'
save.double.width <- function(filename, plot = ggplot2::last_plot(), height = 170) {
  save.plot(filename, plot, width = 170, height)
}

#' Plot a tree as read by ggtree or ape
#'
#' @param tree.data the tree
#' @param tiplab.font.size the size of tip text in ggtree
#' @param ... other parameters to geom_tiplab
#'
#' @export
#' @importFrom rlang .data
#'
plot.tree <- function(tree.data, tiplab.font.size = 2, ...) {
  # Remove underscores for pretty printing
  tree.data$tip.label <- stringr::str_replace_all(tree.data$tip.label, "_", " ")

  if (is.null(tree.data$node.label)) {
    stop("Cannot display bootstrap values, no node labels in tree")
  }

  # Get the complete node labels
  # Separate out bootstrap info
  # Numbers in parentheses are SH-aLRT support (%) / ultrafast bootstrap support (%)
  node.label.values <- data.frame("label" = tree.data$node.label) |>
    tidyr::separate_wider_delim(.data$label,
      delim = "/", names = c("name", "SHaLRT", "UFBoot"),
      too_few = "align_end", too_many = "merge"
    ) |>
    dplyr::mutate(
      UFBoot = suppressWarnings(as.numeric(.data$UFBoot)), # warning not needed, will be NA if missing
      SHaLRT = suppressWarnings(as.numeric(.data$SHaLRT)),
      isSupportedUFBoot = .data$UFBoot >= 95 & !is.na(.data$UFBoot),
      isSupportedSHalRT = .data$SHaLRT >= 80 & !is.na(.data$SHaLRT),
      colour = dplyr::case_when(.data$isSupportedUFBoot & .data$isSupportedSHalRT ~ "black",
        .data$isSupportedUFBoot | .data$isSupportedSHalRT ~ "grey",
        .default = "white"
      )
    )
  ggtree::ggtree(tree.data) +
    ggtree::geom_tree() +
    ggtree::geom_tiplab(size = tiplab.font.size, ggplot2::aes_string(...)) +
    # geom_nodelab(size=2, nudge_x = -0.003, nudge_y = 0.5, hjust=1,  node = "internal")+
    ggtree::geom_nodepoint(size = 1.25, col = "black") +
    ggtree::geom_nodepoint(size = 0.65, col = node.label.values$colour) +
    ggtree::geom_treescale(fontsize = 1.8, y = -1, width = 0.05) +
    ggplot2::coord_cartesian(
      clip = "off",
      ylim = c(-2, length(tree.data$tip.label) + 1)
    ) +
    ggtree::theme_tree() +
    ggplot2::theme(legend.position = "none")
}


#' Given an alignment, calculate KaKs and return in tidy format
#'
#' @param nt.aln.file the alignment file
#'
#' @returns the KaKs ratio in long format
#' @export
#'
calc.kaks <- function(nt.aln.file) {
  seqin.aln <- seqinr::read.alignment(nt.aln.file, format = "fasta")
  kaks.data <- seqinr::kaks(seqin.aln)

  kaks.ratio <- kaks.data$ka / kaks.data$ks

  # Convert to long format and remove pairwise diagonal
  metagMisc::dist2list(kaks.ratio, tri = F)
}

#' Calculate and plot KaKs
#'
#' Given an alignment and order of species, calculate KaKs and make a pairwise plot
#'
#' @param nt.aln.file  the nucleotide alignment
#' @param species.order a vector with the plotting order for species in the file
#' @param kaks.limits min and max for the plot KaKs scale
#'
#' @export
#' @importFrom rlang .data
#'
plot.kaks <- function(nt.aln.file, species.order, kaks.limits = c(0, 1)) {
  # Convert to long format and remove pairwise diagonal
  kaks.pairwise <- calc.kaks(nt.aln.file) |>
    dplyr::mutate(
      col = gsub("_", " ", .data$col),
      row = gsub("_", " ", .data$row),
      col = factor(.data$col, levels = species.order, ordered = T),
      row = factor(.data$row, levels = species.order, ordered = T),
      colnum = as.integer(.data$col),
      rownum = as.integer(.data$row)
    ) |>
    dplyr::filter(.data$rownum < .data$colnum)

  max.y <- max(kaks.pairwise$value, na.rm = TRUE)

  y.limit <- ifelse(is.infinite(max.y), 1, max.y)

  palette.choice <- ggplot2::scale_fill_viridis_c(limits = kaks.limits, direction = -1)

  ggplot2::ggplot(kaks.pairwise, ggplot2::aes(x = .data$col, y = .data$row)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$value)) +
    palette.choice +
    ggplot2::labs(fill = "dNdS") +
    ggplot2::scale_x_discrete(limits = rev) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 5, angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 5),
      axis.title = ggplot2::element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.9, 0.8),
      legend.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 8)
    )
}

#### Other functions #####

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
