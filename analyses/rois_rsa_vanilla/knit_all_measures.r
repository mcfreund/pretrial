library(rmarkdown)
library(here)

wd <- here("analyses", "rois_rsa_vanilla")

for (measure.i in c("corr", "eucl", "neuc")) {
  
  dir.create(file.path(wd, measure.i))
  
  render(
    input = file.path(wd, "rois_rsa_vanilla.rmd"),
    output_file = file.path(wd, paste0("rois_rsa_vanilla_", measure.i, ".html")),
    knit_root_dir = file.path(wd, measure.i)
  )
  
}

