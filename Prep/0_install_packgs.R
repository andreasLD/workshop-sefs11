pkgs <- readLines(file("https://raw.githubusercontent.com/andreasLD/workshop-sefs11/master/Prep/packages_to_install.txt?token=AGEV6IWVS6XCZEIR2WPIYMC5CUWWW", "r"))
str(pkgs)
install.packages(pkgs)
