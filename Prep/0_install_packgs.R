pkgs <- readLines(file("https://raw.githubusercontent.com/andreasLD/workshop-sefs11/master/Prep/packages_to_install.txt", "r"))
str(pkgs)
install.packages(pkgs)
