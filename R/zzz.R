.onLoad <- function(libname, pkgname) {
    require("methods")
}

.onAttach <- function(libname, pkgname) {
  addVigs2WinMenu("tilingArray")
}
