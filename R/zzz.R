.onAttach <- function(libname, pkgname) {
  ## load the compiled code
  if(.Platform$OS.type == "windows" && require("Biobase") && interactive()
        && .Platform$GUI ==  "Rgui"){
        addVigs2WinMenu("tilingArray")
    }
}
