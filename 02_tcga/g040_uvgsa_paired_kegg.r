##g040_uvgsa_paired_kegg.r
##2016-03-23 fgarcia@cipf.es, dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script performs the GSA analysis of the transferred miRNA

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.2.1 (2015-06-18)"
library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"1.2.0"

try (source (".job.r")); try (.job)

setwd (file.path (.job$dir$proces))

## ANNOTATION
load (file.path (.job$dir$proces, "kegg", "annot_kegg_for_paired_data.RData"))
ls ()

## INDEX
load (file.path (.job$dir$proces, "rindex_paired.RData"))
ls ()

tags <- names (rindex)
TAGS <- toupper (tags)
names (TAGS) <- tags
tags
TAGS

################################################################################

length (tags)

## tags <- tags[1:3]    ## 1
## tags <- tags[4:6]    ## 2
## tags <- tags[7:9]    ## 3
## tags <- tags[10:12]  ## 4
## tags <- tags[13:15]  ## 5
## tags <- tags[16:17]  ## 6

################################################################################

## Checking for annot
longitudes <- as.numeric(lapply(annot, length))
summary(longitudes)
hist(longitudes)
table(longitudes == 1)
table(longitudes > 300)


## GSA
setwd (file.path (.job$dir$proces, "kegg"))
dir.create("res_uvgsa_paired", showWarnings = TRUE, recursive = FALSE, mode = "0777")
setwd (file.path (.job$dir$proces, "kegg",  "res_uvgsa_paired"))

for (tag in tags) {
  anotacion <- annotFilter (annot, index = rindex[[tag]], minBlockSize = 2,
                            maxBlockSize = 300)
  res <- try (uvGsa (rindex[[tag]], anotacion, fulltable = TRUE))
  fichero <- paste (tag,".RData", sep = "")
  save (list = "res", file = fichero)
}


###EXIT
warnings ()
sessionInfo ()
q ("no")
