# Functions relating to running external binaries
#' Test if the given binary name is on the PATH
#'
#' Returns true if the binary is on the PATH, false otherwise
#'
#' @param binary the binary to search for
#'
#' @returns true if the binary is on the PATH, false otherwise
#'
#' @examples
#' is.on.PATH("macse")
is.on.PATH <- function(binary) {
  return(Sys.which(binary) != "")
}

#' If the given binary is not on the PATH, stop the script
#'
#' @param binary the binary to test
#'
#' @returns nothing; halts the script if the binary is not found
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

#' Generic wrapper for MACSE
#'
#' Run arbitrary MACSE commands. The subprogram and all options must be
#' specified.
#'
#' @param ... all arguments to macse
#'
#' @export
run.macse <- function(...) {
  check.binary.exists("macse")
  system2("macse", paste(
    ...
  ))
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

#' Run IQTREE2
#'
#' Expects 'iqtree2' on the PATH. By default, runs with 6 threads
#'
#' @param aln.file the alignment
#' @param ... other arguments to IQTREE
#'
#' @returns the path to the treefile
#' @export
#'
run.iqtree2 <- function(aln.file, ...) {
  check.binary.exists("iqtree2")
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

#' Run IQTREE3
#'
#' Expects 'iqtree3' on the PATH. By default, runs with 6 threads
#'
#' @param aln.file the alignment
#' @param ... other arguments to IQTREE
#'
#' @returns the path to the treefile
#' @export
#'
run.iqtree3 <- function(aln.file, ...) {
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
