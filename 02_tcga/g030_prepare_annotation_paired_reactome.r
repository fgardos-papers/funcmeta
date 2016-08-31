##g030_prepare_annotation_paired.r
##2016-03-23 fgarcia@cipf.es, dmontaner@cipf.es
##Meta-analysis of the miRNA data from The Cancer Genome Atlas
##Functional annotation: Reactome

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.0.2 (2013-09-25)"
library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"1.2.0"

try (source (".job.r")); try (.job)

##DATOS
setwd (file.path (.job$dir$proces))
load ("rindex_paired.RData")
ls ()

sapply (rindex, class)
sapply (rindex, head)
sapply (rindex, length)
sapply (rindex, summary)

sapply (rindex, function (x) table (x == 0))

genes <- unique (unlist (lapply (rindex, names)))
length (genes)


## ANNOTATION
reactome <- read.table (file.path (.job$dir$annotation, "mart_export_reactome.txt"), 
                     header = TRUE, sep = "\t", quote = "", na.strings = "", colClasses = "character")
dim (reactome)
reactome[1:3,]

#checkings
table (duplicated (reactome))
table (reactome$Reactome.ID == "")
table (reactome$Associated.Gene.Name == "")

length (unique (reactome$Associated.Gene.Name))
length (unique (reactome$Reactome.ID))

table (reactome$Reactome.ID == "")
table (reactome$Associated.Gene.Name == "")

na.react   <- is.na (reactome$Reactome.ID)
na.gene <- is.na (reactome$Associated.Gene.Name)

table (na.react, na.gene) ## some missing

touse <- !na.react & !na.gene
table (touse)

reactome <- reactome[touse,]
dim (reactome)


## keep just genes in the dataset
touse <- reactome[,"Associated.Gene.Name"] %in% genes  ##most of them are present
table (touse)

reactome <- reactome[touse,]
reactome[1:3,]
length (unique (reactome[,"Associated.Gene.Name"]))
length (unique (reactome[,"Reactome.ID"]))
dim(reactome)



## LIST FORMAT
reactome <- reactome[, c("Associated.Gene.Name", "Reactome.ID")]
reactome[1:3,]

system.time (annot <- annotMat2list (reactome))
length (annot)

size.reactome <- annotList2mat(lapply(annot,length))
colnames(size.reactome) <- c("size", "ID")
size <- as.numeric(size.reactome[,"size"])
hist(size)


## SAVING
save (list = "annot", file = file.path (.job$dir$proces, "reactome", "annot_reactome_for_paired_data.RData"))

###EXIT
warnings ()
sessionInfo ()
q ("no")
