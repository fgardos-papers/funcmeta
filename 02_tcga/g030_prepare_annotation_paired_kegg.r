##g030_prepare_annotation_paired_kegg.r
##2016-03-23 fgarcia@cipf.es, dmontaner@cipf.es
##Meta-analysis of the miRNA data from The Cancer Genome Atlas
##Functional annotation: KEGG

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.2.1 (2015-06-18)"
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
kegg <- read.table (file.path (.job$dir$annotation, "hsa_KEGG.txt"), 
                    header = TRUE, sep = "\t", quote = "", na.strings = "", colClasses = "character")
dim (kegg)
kegg[1:3,]

#checkings
table (duplicated (kegg))
table (kegg$KEGG.ID == "")
table (kegg$Associated.Gene.Name == "")

length (unique (kegg[,"Associated.Gene.Name"]))
length (unique (kegg[,"KEGG.ID"]))


na.react   <- is.na (kegg$KEGG.ID)
na.gene <- is.na (kegg$Associated.Gene.Name)

table (na.react, na.gene) ## some missing

touse <- !na.react & !na.gene
table (touse)

kegg <- kegg[touse,]
dim (kegg)


## keep just genes in the dataset
touse <- kegg[,"Associated.Gene.Name"] %in% genes  ##most of them are present
table (touse)

kegg <- kegg[touse,]
kegg[1:3,]
length (unique (kegg[,"Associated.Gene.Name"]))
length (unique (kegg[,"KEGG.ID"]))
dim(kegg)


## LIST FORMAT
kegg <- kegg[, c("Associated.Gene.Name", "KEGG.ID")]
kegg[1:3,]

system.time (annot <- annotMat2list (kegg))
length (annot)

size.kegg <- annotList2mat(lapply(annot,length))
colnames(size.kegg) <- c("size", "ID")
size <- as.numeric(size.kegg[,"size"])
hist(size)




## SAVING
save (list = "annot", file = file.path (.job$dir$proces, "kegg", "annot_kegg_for_paired_data.RData"))

###EXIT
warnings ()
sessionInfo ()
q ("no")
