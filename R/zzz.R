.onLoad <- function(libname, pkgname){
  global <- new.env() # make a global environment just for loxcoder internal functions
  wrapper_fill_tables() 
}
