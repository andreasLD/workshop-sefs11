pkgs <- readLines(file("https://raw.githubusercontent.com/andreasLD/workshop-sefs11/master/Prep/packages_to_install.txt?token=AGEV6IWVS6XCZEIR2WPIYMC5CUWWW", "r"))
str(pkgs)
install.packages(pkgs)

## piecewiseSEM from GitHub, needs R-Version 3.5 and newer
devtools::install_github("jslefche/piecewiseSEM@devel", build_vignette = TRUE)
devtools::install_github("mfasiolo/mgcViz")
