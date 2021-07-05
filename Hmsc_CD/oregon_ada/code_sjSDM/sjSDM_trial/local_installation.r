

devtools::install_github("https://github.com/TheoreticalEcology/s-jSDM", subdir = "sjSDM")
library(sjSDM)


# Make a pdf
pack <- "sjSDM"
path <- find.package(pack)
system(
  paste(shQuote(file.path(R.home("bin"), "R")),"CMD", "Rd2pdf", shQuote(path)))
devtools::build_manual(path)
