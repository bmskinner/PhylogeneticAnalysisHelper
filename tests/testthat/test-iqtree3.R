
test_that("iqtree3 errors without input", {
  testthat::expect_error(iqtree3())
})

test_that("iqtree3 does not overwrite existing files by default", {

  test.dir <- temp.test.dir()
  dir.create(test.dir)
  on.exit(unlink(test.dir, force = TRUE), add = TRUE)

  aln.file <- file.path(test.dir, "example.aln")
  file.create(aln.file)

  tree.file <- paste0(aln.file, ".treefile")
  file.create(tree.file)

  testthat::expect_warning(iqtree3(aln.file), regexp = "not overwriting. Use overwrite = TRUE to force")

})

test_that("iqtree3 does overwrite existing files with overwrite option", {

  test.dir <- temp.test.dir()
  dir.create(test.dir)
  on.exit(unlink(test.dir, force = TRUE), add = TRUE)

  aln.file <- file.path(test.dir, "example.aln")
  make.example.alignment(aln.file)

  tree.file <- paste0(aln.file, ".treefile")
  file.create(tree.file)

  # Make an empty file. Should be size 0
  testthat::expect_equal(file.size(tree.file), 0)

  # Now the file should be overwritten
  iqtree3(aln.file, overwrite = TRUE)
  testthat::expect_gt(file.size(tree.file), 0)
})

test_that("iqtree3 combines overwrite option with previous varargs", {
  test.dir <- temp.test.dir()
  dir.create(test.dir)
  on.exit(unlink(test.dir, force = TRUE), add = TRUE)

  aln.file <- file.path(test.dir, "example.aln")
  make.example.alignment(aln.file)

  tree.file <- paste0(aln.file, ".treefile")
  file.create(tree.file)

  # Make an empty file. Should be size 0
  testthat::expect_equal(file.size(tree.file), 0)

  # Now the file should be overwritten
  iqtree3(aln.file, "-m JC", overwrite = TRUE)
  testthat::expect_gt(file.size(tree.file), 0)
})

test_that("iqtree3 overwrite takes default value when varargs are included", {
  test.dir <- temp.test.dir()
  dir.create(test.dir)
  on.exit(unlink(test.dir, force = TRUE), add = TRUE)

  aln.file <- file.path(test.dir, "example.aln")
  make.example.alignment(aln.file)

  tree.file <- paste0(aln.file, ".treefile")
  file.create(tree.file)

  # The existing file should not be overwritten
  testthat::expect_warning(iqtree3(aln.file, "-m JC"), regexp = "not overwriting. Use overwrite = TRUE to force")
})
