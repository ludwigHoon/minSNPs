.onAttach <- function(libname, pkgname) {
    packageStartupMessage("The minSNPs version loaded is: ",
        utils::packageVersion("minSNPs"))
}