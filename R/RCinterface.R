 
.onLoad = function(libname, pkgname) {
  library.dynam("Relate", pkgname)
  methods:::bind_activation(TRUE)
}

.onAttach <- function(libname, package) {
  methods:::bind_activation(FALSE)
}
