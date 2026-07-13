#### Functions to plot or export data #####

#' Create an Excel file containing a data frame
#'
#' @param data the data to export
#' @param file.name the name of the file
#'
#' @export
create.xlsx <- function(data, file.name) {
  data <- as.data.frame(data) # ensure not a tibble

  oldOpt <- options()
  options("openxlsx2.dateFormat" = "yyyy-mm-dd")
  options("openxlsx2.datetimeFormat" = "yyyy-mm-dd hh:mm:ss")
  options("openxlsx2.na" = "_openxlsx_NULL")
  options("openxlsx2.maxWidth" = 50)

  openxlsx2::wb_workbook() |>
    openxlsx2::wb_add_worksheet("Sheet 1") |>
    openxlsx2::wb_add_data("Sheet 1", data, row_names = FALSE) |>
    openxlsx2::wb_freeze_pane(sheet = "Sheet 1", first_row = TRUE, first_col = TRUE) |>
    openxlsx2::wb_add_filter(sheet = "Sheet 1", rows = 1, cols = 1:ncol(data)) |>
    openxlsx2::wb_set_col_widths(sheet = "Sheet 1", cols = 1:ncol(data), widths = "auto") |>
    openxlsx2::wb_save(file = file.name)
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
#' @param with.svg also save the plot in SVG format if TRUE
#'
#' @export
#'
save.plot <- function(filename, plot = ggplot2::last_plot(), width, height, with.svg=TRUE) {
  ggplot2::ggsave(filename, plot, dpi = 300, units = "mm", width = width, height = height)
  if(with.svg){
    ggplot2::ggsave(stringr::str_replace(filename, ".png$", ".svg"), plot,
                    dpi = 300,
                    units = "mm", width = width, height = height
    )
  }
}

#' Save a double column plot (170mm) at 300dpi
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

#' Plot a tree as read by ggtree or ape.
#'
#' The tree may or may not have bootstrap
#' values in the node labels. Either zero, one or two bootstrap values per node
#' are parsed, assuming bootstrap values are delimited by `/`.
#'
#' @param tree the tree
#' @param tiplab.font.size the size of tip text in ggtree
#' @param ... other parameters to geom_tiplab
#'
#' @importFrom rlang .data
#'
#' @export plot.tree
plot.tree <- function(tree, tiplab.font.size = 2, ...) {
  # Remove underscores for pretty printing
  tree$tip.label <- stringr::str_replace_all(tree$tip.label, "_", " ")

  base.tree <- suppressWarnings(ggtree::ggtree(tree)) +
    ggtree::geom_tree() +
    # Draw tips. Pass through extra aesthetic arguments in ...
    ggtree::geom_tiplab(size = tiplab.font.size, ggplot2::aes_string(...)) +
    ggtree::geom_treescale(fontsize = 1.8, y = -1, width = 0.05) +
    ggtree::theme_tree() +
    ggplot2::theme(legend.position = "none")

  # No node labels means no bootstrap values. Return the tree as-is.
  if (is.null(tree$node.label)) {
    message("Cannot display bootstrap values, no node labels in tree")
    return(base.tree)
  }

  # Get the complete node labels
  node.label.values <- data.frame("label" = tree$node.label)
  # Bootstraps may be written within the node label before the branch length. There
  # may be more than one bootstrap value for a node (e.g. SH-aLRT and UFBoot in
  # IQTREE). Example: "Node9/98.1/99". The root node will be the empty string.

  has.bootstrap <- any(stringr::str_detect(node.label.values$label, "/"))
  if (!has.bootstrap) {
    message("No bootstrap separator string detected\n")
    return(base.tree)
  }

  n.bootstrap.values <- max(stringr::str_count(node.label.values$label, "/"))+1
  message("Detected ", n.bootstrap.values, " bootstrap values per node\n")

  if (n.bootstrap.values == 1) {
    # Single bootstrap node colour
    bootstrap.colours <- grDevices::grey.colors(n = 101, start = 0, end = 1, rev = TRUE)
    bootstrap.values <- node.label.values |>
      tidyr::separate_wider_delim(.data$label,
        delim = "/", names = c("label", "Bootstrap"),
        too_few = "align_end", too_many = "merge"
      ) |>
      dplyr::mutate(
        Bootstrap = suppressWarnings(as.numeric(.data$Bootstrap)),
        Colour = bootstrap.colours[as.integer(Bootstrap) + 1]
      ) # warning not needed, will be NA if missing

    base.tree <- base.tree +
      # Draw coloured circles over each node with bootstrap support
      ggtree::geom_nodepoint(size = 1.25, col = "black") + # big black circle (will make a black outline)
      ggtree::geom_nodepoint(size = 0.65, col = bootstrap.values$Colour) # small coloured circle on top
  }

  if (n.bootstrap.values == 2) {
    # Two bootstrap values; assume SH-aLRT and UFBoot
    node.label.values <- node.label.values |>
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
    base.tree <- base.tree +
      # Draw coloured circles over each node with bootstrap support
      ggtree::geom_nodepoint(size = 1.25, col = "black") + # big black circle (will make a black outline)
      ggtree::geom_nodepoint(size = 0.65, col = node.label.values$colour) # small coloured circle on top
  }

  return(base.tree)
  # Optional - we can set the y axis to exactly size needed if space is tight
  # ggplot2::coord_cartesian(
  #   clip = "off",
  #   ylim = c(-2, length(tree.data$tip.label) + 1)
  # ) +
}


#' Given a nucleitide alignment file, calculate the KaKs and return in tidy format.
#' The sequence alignment is expected to contain coding sequences only.
#'
#' @param nt.aln.file the alignment file in FASTA format.
#'
#' @returns the KaKs ratio in long format
#' @usage calc.kaks(nt.aln.file)
#' @export calc.kaks
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
#' @usage plot.kaks(nt.aln.file,
#'  species.order = NA,
#'  kaks.limits = c(0, 1))
#' @export plot.kaks
#' @importFrom rlang .data
#'
plot.kaks <- function(nt.aln.file, species.order = NA, kaks.limits = c(0, 1)) {
  # Calculate the KaKs data
  kaks.data <- calc.kaks(nt.aln.file)

  # Set default ordering if not given
  species.order <- ifelse(is.na(species.order), sort(unique(kaks.data$col)), species.order)

  # Convert to long format and remove pairwise diagonal
  kaks.pairwise <- kaks.data |>
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
