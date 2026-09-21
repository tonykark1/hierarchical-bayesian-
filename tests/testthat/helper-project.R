# Load project R functions once for the test suite.
# testthat evaluates files from tests/testthat, so resolve the repository root
# explicitly instead of relying on the current working directory inside tests.

project_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

r_files <- list.files(
  file.path(project_root, "R"),
  pattern = "\\.[Rr]$",
  full.names = TRUE
)

for (f in sort(r_files)) {
  sys.source(f, envir = globalenv())
}
