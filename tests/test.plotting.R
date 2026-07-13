# Test cases
# example.tree <- ape::read.tree("tests/example.treefile") |>
#   reroot.tree(c("Xenopus_ZFX.L", "Xenopus_ZFX.S"), 0.15)
#
# plot.tree(example.tree)
#
# # Remove one bootstrap value
# example.tree.single.bootstrap <- ape::read.tree("tests/example.treefile")
# example.tree.single.bootstrap$node.label <- stringr::str_replace_all(
#   example.tree.single.bootstrap$node.label,
#   "([\\w\\d]+/[\\d\\.]+)/[\\d\\.]+", "\\1"
# )
#
# example.tree.single.bootstrap |>
#   reroot.tree(c("Xenopus_ZFX.L", "Xenopus_ZFX.S"), 0.15) |>
#   plot.tree()
