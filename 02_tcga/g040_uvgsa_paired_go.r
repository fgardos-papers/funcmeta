##g040_uvgsa_paired.r
##2016-03-23 fgarcia@cipf.es, dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script performs the GSA analysis of the transferred miRNA

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.2.1 (2015-06-18)"
library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"1.2.0"
#library (mirbaseID); packageDescription ("mirbaseID", fields = "Version") #"0.0.2"
#help (package = mdgsa)
#help (package = mirbaseID)

try (source (".job.r")); try (.job)

## ANNOTATION
load (file.path (.job$dir$proces, "go", "annot_for_paired_data.RData"))
ls ()

ontologias <- names (annot)
ontologias

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

## GSA
setwd (file.path (.job$dir$proces, "go",  "res_uvgsa_paired"))

for (tag in tags) {
    for (onto in ontologias) {
        cat ("\n=============== ", TAGS[tag], ":", onto, " ===============\n")
        
        anotacion <- annotFilter (annot[[onto]], index = rindex[[tag]], minBlockSize = 10, maxBlockSize = 300)

        if (.job$testmode) anotacion <- anotacion[1:3]
        
        res <- try (uvGsa (rindex[[tag]], anotacion, fulltable = TRUE))
        
        fichero <- paste (tag, "_", onto, ".RData", sep = "")
        save (list = "res", file = fichero)
    }
}

###EXIT
warnings ()
sessionInfo ()
q ("no")
