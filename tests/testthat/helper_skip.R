skip_on_32bit <- function ()
{
  if (.Machine$sizeof.pointer != 4L) {
    return(invisible(TRUE))
  }
  skip("32 bit")
}
