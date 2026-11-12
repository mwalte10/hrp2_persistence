# =============================================================================
# provision.R
#
# Installs project package dependencies into the hipercow cluster environment.
# Invoked via:
#   hipercow_provision(method = "script", script = "provision.R")
# (see src/04_calibrate_cluster.R, Section A).
#
# Refs are read from pkgdepends.txt so the two files can't drift apart;
# blank lines and lines starting with "#" are ignored.
# =============================================================================

if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages(
    "pak",
    repos = sprintf(
      "https://r-lib.github.io/p/pak/stable/%s/%s/%s",
      .Platform$pkgType, R.Version()$os, R.Version()$arch
    )
  )
}

refs <- readLines("pkgdepends.txt")
refs <- trimws(refs)
refs <- refs[nzchar(refs) & !startsWith(refs, "#")]

pak::pkg_install(refs)
