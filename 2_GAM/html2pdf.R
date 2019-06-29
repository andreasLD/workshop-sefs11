# script to convert .html presentation to .pdf
require(webshot)

fl = file.path(getwd(), '2_GAM/GAM.html')
file_name = paste0("file://", normalizePath(fl))
webshot(file_name, gsub('html', 'pdf', fl))
