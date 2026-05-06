make.example.alignment <- function(file.name) {
  cat(">No305",
    "NTTCGAAAAACACACCCACTACTAAAANTTATCAGTCACT",
    ">No304",
    "ATTCGAAAAACACACCCACTACTAAAAATTATCAACCACT",
    ">No306",
    "ATTCGAAAAACACACCCACTACTAAAAATTATCAATCACT",
    file = file.name, sep = "\n"
  )
}


test_that("iqtree3 errors without input", {
  testthat::expect_error(run.iqtree3())
})

test_that("iqtree3 does not overwrite existing files by default", {
  on.exit(sapply(list.files(pattern = "example.*"), file.remove), add = TRUE)
  file.create("example.aln")
  file.create("example.aln.treefile")
  testthat::expect_warning(run.iqtree3("example.aln"), "Tree 'example.aln.treefile' exists, not overwriting. Use overwrite = TRUE to force")
})

test_that("iqtree3 does overwrite existing files with overwrite option", {
  on.exit(sapply(list.files(pattern = "example.*"), file.remove), add = TRUE)

  make.example.alignment("example.aln")

  # Make an empty file. Should be size 0
  file.create("example.aln.treefile")
  testthat::expect_equal(file.size("example.aln.treefile"), 0)

  # Now the file should be overwritten
  run.iqtree3("example.aln", overwrite = TRUE)
  testthat::expect_gt(file.size("example.aln.treefile"), 0)
})

test_that("iqtree3 combines overwrite option with varargs", {
  on.exit(sapply(list.files(pattern = "example.*"), file.remove), add = TRUE)

  make.example.alignment("example.aln")

  # Make an empty file. Should be size 0
  file.create("example.aln.treefile")
  testthat::expect_equal(file.size("example.aln.treefile"), 0)

  # Now the file should be overwritten
  run.iqtree3("example.aln", "-m JC", overwrite = TRUE)
  testthat::expect_gt(file.size("example.aln.treefile"), 0)
})

test_that("iqtree3 overwrite takes default value when varargs included", {
  on.exit(sapply(list.files(pattern = "example.*"), file.remove), add = TRUE)

  make.example.alignment("example.aln")

  # Now the file should be overwritten
  run.iqtree3("example.aln", "-m JC")
  testthat::expect_true(file.exists("example.aln.treefile"))
})
