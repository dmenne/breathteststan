skip_on_32bit <- function ()
{
  if (.Machine$sizeof.pointer != 8L) {
    return(invisible(TRUE))
  }
  skip("32 bit")
}
