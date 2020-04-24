check_dependencies <- function(func_name, pkg) {
  need_MSEtool <- !requireNamespace("MSEtool", quietly = TRUE)
  need_pkg <- !requireNamespace(pkg, quietly = TRUE)

  if(need_pkg) {
    if(pkg == "spict") pkg_repo <- "DTUAqua/spict/spict"
    if(pkg == "stockassessment") pkg_repo <- "fishfollower/SAM/stockassessment"
    if(pkg == "LIME") pkg_repo <- "merrillrudd/LIME"
  }

  if(need_MSEtool || need_pkg) {
    msg_out <- paste0("The following package(s) are needed to run ", func_name, " function:")
    if(need_MSEtool) msg_out <- paste0(msg_out, '\nMSEtool - download from CRAN via install.packages("MSEtool")')

    if(need_pkg) msg_out <- paste0(msg_out, "\n", pkg, ' - download from Github via devtools::install_github("', pkg_repo, '")')

    stop(msg_out, call. = FALSE)
  }
  return(invisible())
}
