library(rmarkdown)
library(here)

wd <- here("analyses", "rois_rsa_vanilla")

for (measure.i in c("corr", "eucl", "neuc")) {
  
  for (norma.i in c("raw", "prw")) {
    
    # if (paste0(measure.i, "_", norma.i) %in% c("corr_raw", "corr_prw", "eucl_raw")) next
    
    dir.create(file.path(wd, paste0(measure.i, "_", norma.i)))
    
    render(
      input = file.path(wd, "rois_rsa_vanilla.rmd"),
      output_file = file.path(wd, paste0("rois_rsa_vanilla_", measure.i, "_", norma.i, ".html")),
      knit_root_dir = file.path(wd, paste0(measure.i, "_", norma.i))
    )
    
  }
  
}

