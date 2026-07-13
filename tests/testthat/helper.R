# Common test functions
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

# Create a temporary directory path within the session temp dir
temp.test.dir <- function(){
  file.path(tempdir(), paste(letters[runif(7, min=1, max=26)], collapse=""))
}
