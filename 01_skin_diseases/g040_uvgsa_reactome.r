##g040_uvgsa.r
##2016-03-23 fgarcia@cipf.es, dmontaner@cipf.es
##Meta-analysis of skin diseases 
##Functional annotation: Reactome

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.2.1 (2015-06-18)"
library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"1.2.0"
library(Biobase); packageDescription ("Biobase", fields = "Version")

try (source (".job.r")); try (.job)


##ANNOTATION
reactome <- read.table (file.path (.job$dir$annotation, "mart_export_reactome_ensembl_gene84.txt"), 
                    header = TRUE, sep = "\t", quote = "", na.strings = "", colClasses = "character")
dim (reactome)
reactome[1:3,]

#checkings
table (duplicated (reactome))
table (reactome$ReactomeID == "")
table (reactome$EnsemblGeneID == "")

length (unique (reactome$ReactomeID))
length (unique (reactome$EnsemblGeneID))

na.react   <- is.na (reactome$ReactomeID)
na.gene <- is.na (reactome$EnsemblGeneID)
table (na.react, na.gene) ## some missing
touse <- !na.react & !na.gene
table (touse)

reactome <- reactome[touse,]
dim (reactome)
length (unique (reactome$ReactomeID))
length (unique (reactome$EnsemblGeneID))

anot.mat <- reactome


# ## keep just genes in the dataset
# touse <- kegg[,"Associated.Gene.Name"] %in% genes  ##most of them are present
# table (touse)
# 
# kegg <- kegg[touse,]
# kegg[1:3,]
# length (unique (kegg[,"Associated.Gene.Name"]))
# length (unique (kegg[,"KEGG.ID"]))
# dim(kegg)



##DATOS
setwd(file.path(.job$dir$proces, "difexp"))
ficheros <- dir(pattern= "GSE")
ficheros

for (fi in ficheros){
  setwd(file.path(.job$dir$proces, "difexp"))
  cat ("\n ============================= \n")
  print (fi)
  load (fi)
  
  table (rownames (fit$p.value) == rownames (fit$t))
  
  ##RANKING INDEX
  rindex0 <- pval2index (pval=fit$p.value[,2], sign= fit$t[,2])
  #transforma p.value y su valor asociado de t en un ranking
  summary (rindex0)
  rindex0[1:3]
  length (rindex0)
  
  table (rownames (fit$p.value) == names (rindex0))
  
  ##TRANSFORM
  rindex <- indexTransform (index = rindex0, method = "normalize")
  #transforma el ranking para que su distribucion sea adecuada como variable independiente de un modelo de regresion logistica
  
  touse <- anot.mat[,1] %in% names (rindex)
  table (touse)
  anot.mat <- anot.mat[touse,]
  
  anot <- annotMat2list (anot.mat)
  #convierte una matriz de anotacion con dos columnas (la primera con los id de los genes y la segunda
  #con los id de la anotacion) en una lista de anotacion. Cada GO con los id_ensembl asociados
  head (anot)
  
#   ##propagate annotation 
#   anot <- propagateGO (anot, verbose = TRUE)
#   #los genes anotados bajo un GO heredan la anotacion de todos los ancestros
#   
#   ##filter annotation
#   anot <- annotFilter (annot = anot, index = rindex, minBlockSize = 10, maxBlockSize = 500)
#   #filtra los GO estableciendo un minimo y maximo de id_ensembl asociados a cada GO
#   length (anot)
#   
  ## gene set analysis 21h
  res <- uvGsa (index = rindex, annot = anot, fulltable = TRUE)
  #realiza un analisis con un modelo de regresion logistica
  #Warning: The analysis did not converge for some blocks. You may re-run uvGsa using 'fulltable = TRUE' to find them.
  res[1:3,]
  save  (res, file = file.path (.job$dir$proces, "gsa", fi))
}


###EXIT
warnings ()
sessionInfo ()
q ("no")
