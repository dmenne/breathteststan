set.seed(4711)
Sys.unsetenv("R_TESTS")

skip_on_32bit <- function ()
{
  if (.Machine$sizeof.pointer != 4L) {
    return(invisible(TRUE))
  }
  skip("32 bit")
}
