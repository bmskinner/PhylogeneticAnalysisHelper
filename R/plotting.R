#### Functions to plot or export data #####

#' Create an Excel file containing a data frame
#'
#' @param data the data to export
#' @param file.name the name of the file
#' @param cols.to.fixed.size.font integer vector of column indexes to set in monospace font
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
