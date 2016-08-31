##g050_meta_analysis_reactome.r
##2016-01-04 fgarcia@cipf.es, antoi@alumni.uv.es
##Functional Meta-Analysis of transcriptomic studies
##Annotation: REACTOME 



# STEP 1. Preparing input for meta-analysis: LOR and SD matrix
#################################################################
## starting
date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
R.version.string 
try (source (".job.r")); try (.job)

# clean the working space
rm (list = ls ())


### load libraries
# We need to install and load these packages:
library(Biobase); packageDescription("Biobase", fields = "Version")
library(metafor); packageDescription("metafor", fields = "Version")
library(mdgsa); packageDescription("mdgsa", fields = "Version")
#library(RamiGO); packageDescription("RamiGO", fields = "Version")
library(ggplot2); packageDescription("ggplot2", fields = "Version")
library(reshape); packageDescription("reshape", fields = "Version")


### previously for each study, we need results GSA from mdgsa in a file .RData including this structure (dataframe): 
####            N        lor        pval      padj        sd          t conv
#### GO:0000002 12  0.4373106 0.128203007 1.0000000 0.2874434  1.5213797    1
#### GO:0000018 29  0.2177107 0.238623945 1.0000000 0.1847325  1.1785188    1



### load results from gene set analysis for each study
setwd (file.path (.job$dir$proces, "gsa", "reactome"))
ficheros <- dir(pattern= ".RData")
ficheros

# we search a list including all unique reactome IDs for all studies
reactomeuniq<-NULL
for (fi in ficheros){
  load (fi)
  reactomeuniq <- c(reactomeuniq, rownames (res))
}
length (reactomeuniq)
reactomeuniq <- unique (reactomeuniq)
length (reactomeuniq)
reactomeuniq <- sort (reactomeuniq)


### generating matrix with all LOR for all studies
mat.lor <- matrix (NA, nrow = length (reactomeuniq), ncol = length (ficheros))
rownames (mat.lor) <- reactomeuniq
colnames (mat.lor) <- strsplit (ficheros, ".RData")
head (mat.lor[, 1:5])

for (fi in ficheros){
  load (fi)
  co <- as.character(strsplit (fi, ".RData"))
  lor <- res$lor
  names (lor) <- (rownames (res))
  mat.lor[, co] <- lor[rownames(mat.lor)] 
}

head (mat.lor[, 1:5])
table (is.na(mat.lor))
dim (mat.lor)
summary(mat.lor)


### generating matrix with all SD for all studies
mat.sd <- matrix (NA, nrow = length (reactomeuniq), ncol = length (ficheros))
rownames (mat.sd) <- reactomeuniq
colnames (mat.sd) <- strsplit (ficheros, ".RData")
head (mat.sd[, 1:5])

for (fi in ficheros){
  load (fi)
  co <- as.character(strsplit (fi, ".RData"))
  sd <- res$sd
  names (sd) <- (rownames (res))
  head (sd)
  mat.sd[, co] <- sd[rownames(mat.sd)] 
}

head (mat.sd)
table (is.na(mat.sd))
dim (mat.sd)
summary(mat.sd)


### generating matrix with all adjusted p values for all studies
mat.adjp <- matrix (NA, nrow = length (reactomeuniq), ncol = length (ficheros))
rownames (mat.adjp) <- reactomeuniq
colnames (mat.adjp) <- strsplit (ficheros, ".RData")
head (mat.adjp[, 1:5])

for (fi in ficheros){
  load (fi)
  co <- as.character(strsplit (fi, ".RData"))
  adjp <- res$padj
  names (adjp) <- (rownames (res))
  head (adjp)
  mat.adjp[, co] <- adjp[rownames(mat.adjp)] 
}

head (mat.adjp)
table (is.na(mat.adjp))
dim (mat.adjp)
table(mat.adjp < 0.05)

mat.adj2 <- (mat.adjp < 0.05)
head(mat.adj2)
total <- apply(mat.adj2, 1, sum)
head(total)
summary(total)
table(total)






# STEP 2. Meta-analysis for functional terms
#################################################################

#In general, we suppose between-study variance is non-zero.
#there are different methods to estimate this variance:
#DL (Dersimonian-Laird), REML (Restricted maximum-likelihood, default)....
#result.lor <- rma(yi = mat.lor[1, ], vi = mat.sd[1, ], method = "DL") 
# DerSimonian-Laird. Y -> log(OR) V -> var(log(OR)

#We also evaluate Fixed Effects (FE)


# explore the function to do the meta-analysis
#?rma

# parameters to select
# methods = c("DL", "HE", "HS", "SJ", "ML", "REML", "EB", "PM", "FE")  #methods for meta-analysis
methods = c("DL", "HE", "HS", "SJ",                 "FE")  #methods for meta-analysis

corte = 0.05  #threshold to detect significant results
or    = 0.5   #threshold to detect OR 
adj.p.value = "fdr"  # method to adjust p-values

# preparing a matrix for results
res <- matrix (NA, nrow = length (methods), ncol = 6)
colnames(res) <- c("over", "under", "sig.over", "sig.under", "sig.or.over", "sig.or.under")
rownames(res) <- methods
res
# important: we have to indicate our variability measure is standard deviance and not variance
# from "sei" (by defalt is "vi"):

# meta-analysis
for (i in methods){
  meta_analisis <- lapply(1:length(rownames(mat.lor)),
                          function(x){yi = rma(mat.lor[x, ], sei =mat.sd[x, ],
                                          method = i)})
  # sometimes for   methods there are converge problems: I have to increase the number of iter
  #(http://www.metafor-project.org/doku.php/tips:convergence_problems_rma)
  #   length(meta_analisis)
  #   meta_analisis[[1]]
  names (meta_analisis) <- rownames(mat.lor)
  result_meta <- as.data.frame(do.call("rbind", lapply(meta_analisis,function(x)
    {c(x$ci.lb, x$b, x$ci.ub, x$pval, x$QE, x$QEp, x$se, x$tau2, x$I2, x$H2) })))
  colnames(result_meta) <- c("lower_bound", "summary_LOR", "upper_bound", "pvalue", 
                             "QE", "QEp", "SE", "tau2", "I2", "H2")
  p.adjust <- p.adjust(result_meta[, 4], method= adj.p.value)  
  result_meta <- round(cbind (result_meta, p.adjust),3)
  
  if (file.exists(file.path (.job$dir$code, "results",  "reactome", "files"))){
    setwd(file.path(.job$dir$code, "results",  "reactome", "files"))
  } else {
    dir.create(file.path(.job$dir$code, "results",  "reactome", "files"), recursive = T)
    setwd(file.path(.job$dir$code, "results",  "reactome", "files"))    
  }
  #save  (list = c("meta_analisis", "result_meta"), file = sprintf("meta_result_%s.RData", i))
  
  res[i, "over"]      <-  sum(result_meta[, "summary_LOR"] > 0)
  res[i, "under"]    <-  sum(result_meta[, "summary_LOR"] < 0)
  res[i,"sig.over"]   <-  sum(result_meta[, "summary_LOR"] > 0 & result_meta[,"p.adjust"] < corte)
  res[i,"sig.under"] <-  sum(result_meta[, "summary_LOR"] < 0 & result_meta[,"p.adjust"] < corte)
  res[i,"sig.or.over"]   <-  sum(result_meta[, "summary_LOR"] >  or & result_meta[,"p.adjust"] < corte)
  res[i,"sig.or.under"] <-  sum(result_meta[, "summary_LOR"] < -or & result_meta[,"p.adjust"] < corte)
  

  # Getting common IDs for all models 
  if (i == "DL"){
    DL <- result_meta    
    DL[, "ID"]   <- rownames(DL)
    DL <- DL[c("ID", "lower_bound", "summary_LOR", "upper_bound", "pvalue",
               "p.adjust", "QE", "QEp", "SE", "tau2", "I2", "H2" )]
    write.table(DL, "all.results.DL.txt", sep = "\t", quote = F, row.names = F)
    DLs = subset(DL, DL[,"p.adjust"] < corte)
    write.table(DLs, "sig.results.DL.txt", sep = "\t", quote = F, row.names = F)
    idDL = rownames(DLs)
   #these lines, only for DL to include in final report
    DLs_or = DLs[abs(DLs$summary_LOR) > or,]
    write.table(DLs_or, "sig.or.results.DL.txt", sep = "\t", quote = F, row.names = F)
  
  }
  if (i == "HE"){
    HE <- result_meta    
    HE[, "ID"]   <- rownames(HE)
    HE <- HE[c("ID", "lower_bound", "summary_LOR", "upper_bound", "pvalue",
               "p.adjust","QE", "QEp", "SE", "tau2", "I2", "H2" )]
    write.table(HE, "all.results.HE.txt", sep = "\t", quote = F, row.names = F)
    HEs = subset(HE, HE[,"p.adjust"] < corte)
    write.table(HEs, "sig.results.HE.txt", sep = "\t", quote = F, row.names = F)
    idHE = rownames(HEs)
  }
  if (i == "HS"){
    HS <- result_meta
    HS[, "ID"]   <- rownames(HS)
    HS <- HS[c("ID",  "lower_bound", "summary_LOR", "upper_bound", "pvalue",
               "p.adjust","QE", "QEp", "SE", "tau2", "I2", "H2" )]
    write.table(HS, "all.results.HS.txt", sep = "\t", quote = F, row.names = F)
    HSs = subset(HS, HS[,"p.adjust"] < corte)
    write.table(HSs, "sig.results.HS.txt", sep = "\t", quote = F, row.names = F)
    idHS = rownames(HSs)
  }
  if (i == "SJ"){
    SJ <- result_meta
    SJ[, "ID"]   <- rownames(SJ)
    SJ <- SJ[c("ID", "lower_bound", "summary_LOR", "upper_bound", "pvalue",
               "p.adjust","QE", "QEp", "SE", "tau2", "I2", "H2" )]
    write.table(SJ, "all.results.SJ.txt", sep = "\t", quote = F, row.names = F)
    SJs = subset(SJ, SJ[,"p.adjust"] < corte)
    write.table(SJs, "sig.results.SJ.txt", sep = "\t", quote = F, row.names = F)
    idSJ = rownames(SJs)
  }
#   if (i == "ML"){
#     ML <- result_meta
#     ML[,"name"] <- getKEGGnames(substr(rownames(ML),4,8))
#     ML[, "ID"]   <- rownames(ML)
#     ML <- ML[c("ID", "name", "lower_bound", "summary_LOR", "upper_bound", "pvalue",
#                "p.adjust","QE", "QEp", "SE", "tau2", "I2", "H2" )]
#     write.table(ML, "all.results.ML.txt", sep = "\t", quote = F, row.names = F)
#     MLs = subset(ML, ML[,"p.adjust"] < corte)
#     write.table(MLs, "sig.results.ML.txt", sep = "\t", quote = F, row.names = F)
#     idML = rownames(MLs)
#   }
#   if (i == "REML"){
#     REML <- result_meta
#     REML[,"name"] <- getKEGGnames(substr(rownames(REML),4,8))
#     REML[, "ID"]   <- rownames(REML)
#     REML <- REML[c("ID", "name", "lower_bound", "summary_LOR", "upper_bound", "pvalue",
#                    "p.adjust","QE", "QEp", "SE", "tau2", "I2", "H2" )]
#     write.table(REML, "all.results.REML.txt", sep = "\t", quote = F, row.names = F)
#     REMLs = subset(REML, REML[,"p.adjust"] < corte)
#     write.table(REMLs, "sig.results.REML.txt", sep = "\t", quote = F, row.names = F)
#     idREML = rownames(REMLs)
#   }
#   if (i == "EB"){
#     EB <- result_meta
#     EB[,"name"] <- getKEGGnames(substr(rownames(EB),4,8))
#     EB[, "ID"]   <- rownames(EB)
#     EB <- EB[c("ID", "name", "lower_bound", "summary_LOR", "upper_bound", "pvalue",
#                "p.adjust","QE", "QEp", "SE", "tau2", "I2", "H2" )]
#     write.table(EB, "all.results.EB.txt", sep = "\t", quote = F, row.names = F)
#     EBs = subset(EB, EB[,"p.adjust"] < corte)
#     write.table(EBs, "sig.results.EB.txt", sep = "\t", quote = F, row.names = F)
#     idEB = rownames(EBs)
#   }
  if (i == "FE"){
    FE <- result_meta
    FE[, "ID"]   <- rownames(FE)
    FE <- FE[c("ID", "lower_bound", "summary_LOR", "upper_bound", "pvalue",
               "p.adjust","QE", "QEp", "SE", "tau2", "I2", "H2" )]
    write.table(FE, "all.results.FE.txt", sep = "\t", quote = F, row.names = F)
    FEs = subset(FE, FE[,"p.adjust"] < corte)
    write.table(FEs, "sig.results.FE.txt", sep = "\t", quote = F, row.names = F)
    idFE = rownames(FEs)
  }
}



##GLOBAL RESULTS
setwd (file.path (.job$dir$code, "results","reactome","files"))


# A. All results for different models
print(res)
cat ("method\t", file = "res_all.txt")
write.table (res, file = "res_all.txt",
             append = TRUE, quote = FALSE, sep = "\t", 
             row.names = TRUE, col.names = TRUE)


# B. Intersection
# intersectmodels = Reduce(intersect, list(idDL, idHE, idHS, idSJ, idML, idREML, idEB, idPM, idFE))
 intersectmodels =  Reduce(intersect,  list(idDL, idHE, idHS, idSJ,                  idFE))
length(intersectmodels)
write.table (intersectmodels, file = "res_intersection.txt",
             quote = FALSE, sep = "\t", row.names = F, col.names = F)
# # to print in a file:
# cat(intersectmodels, file = "intersection_models.txt", sep = "\n", append = TRUE) 
# cat("\nList of identifiers generated from the intersection of  generated models.")



# C. Union
#unionmodels = Reduce(c, list(idDL, idHE, idHS, idSJ, idML, idREML, idEB, idPM, idFE))
 unionmodels = Reduce(c, list(idDL, idHE, idHS, idSJ,                idFE))

mat.union <- as.data.frame(table(unionmodels))
dim(mat.union)
mat.union <- mat.union[order(mat.union$Freq, decreasing = T),]
cat ("ID\tFreq\n", file = "res_union.txt")
write.table (mat.union, file = "res_union.txt",append = TRUE,
             quote = FALSE, sep = "\t", row.names = F, col.names = F)




# STEP 3. GRAPHICAL REPRESENTATION
#################################################################

if (file.exists(file.path (.job$dir$code, "results",  "reactome",  "plots"))){
  setwd(file.path(.job$dir$code, "results",  "reactome",  "plots"))
} else {
  dir.create(file.path(.job$dir$code, "results",  "reactome",  "plots"), recursive = T)
  setwd(file.path(.job$dir$code, "results",  "reactome",  "plots"))    
}




## 3.1. plot to evaluate heterogeneity for all models 
######################################################

# methods<-  c(rep("DL", nrow(DL)),rep("HE", nrow(HE)),rep("HS", nrow(HS)),rep("SJ", nrow(SJ)), 
#             rep("ML", nrow(ML)),rep("REML", nrow(REML)),rep("EB", nrow(EB)),
#             rep("PM", nrow(PM)),rep("FE", nrow(FE)))
methods<-  c(rep("DL", nrow(DL)),rep("HE", nrow(HE)),rep("HS", nrow(HS)),rep("SJ", nrow(SJ)), 
                                                    rep("FE", nrow(FE)) )
# QEp   <- c(DL$QEp, HE$QEp, HS$QEp, SJ$QEp, ML$QEp, REML$QEp, EB$QEp, PM$QEp, FE$QEp)
QEp      <- c(DL$QEp, HE$QEp, HS$QEp, SJ$QEp,                     FE$QEp)
SE       <- c(DL$SE, HE$SE, HS$SE, SJ$SE,                     FE$SE)
I2       <- c(DL$I2, HE$I2, HS$I2, SJ$I2,                     FE$I2)
H2       <- c(DL$H2, HE$H2, HS$H2, SJ$H2,                     FE$H2)
tau2     <- c(DL$tau2, HE$tau2, HS$tau2, SJ$tau2,                    FE$tau2)

datos <- data.frame (cbind(as.factor(methods), QEp, SE, I2, H2, tau2))
head(datos)
dim(datos)
summary(datos)
param <- c("QEp", "SE", "I2", "H2", "tau2")
x.por <- 3
y.por <- 2


# We drawn separately QEp because we need abline for 0.05
jpeg (filename = "KEGG_het_QEp.jpeg",  
      width = 480 *x.por ,   height = 480*y.por,
      pointsize = 150, quality = 100) 
p <- ggplot(datos, aes(methods, datos[, "QEp"]))
p <- p + geom_boxplot(col = "blue")
p <- p + xlab("Métodos") 
p <- p + ylab("QEp")
p <- p + theme(axis.title.x= element_text(colour= "black", size = 40 )) #to change size of title in axis x
p <- p + theme(axis.title.y= element_text(colour= "black", size = 40 )) #to change size of title in axis y
p <- p + theme(axis.text.x = element_text(colour= "black", size = 30 )) #to change size of numbers
p <- p + theme(axis.text.y = element_text(colour= "black", size = 30 )) #to change size of numbers
p <- p + geom_hline(yintercept = 0.05, col = "red")
print(p)
dev.off ()

#Now the rest:
for (i in 2:length(param)){
  jpeg (filename = paste("KEGG_het_", param[i], ".jpeg", sep = ""),  
        width = 480 *x.por ,   height = 480*y.por,
        pointsize = 150, quality = 100) 
  p <- ggplot(datos, aes(methods, datos[, param[i]]))
  p <- p + geom_boxplot(col = "blue")
  p <- p + xlab("Métodos") 
  p <- p + ylab(param[i])
  p <- p + theme(axis.title.x= element_text(colour= "black", size = 40 )) #to change size of title in axis x
  p <- p + theme(axis.title.y= element_text(colour= "black", size = 40 )) #to change size of title in axis y
  p <- p + theme(axis.text.x = element_text(colour= "black", size = 30 )) #to change size of numbers
  p <- p + theme(axis.text.y = element_text(colour= "black", size = 30 )) #to change size of numbers
#   p <- p + geom_hline(yintercept = 0.05, col = "red")
  print(p)
  dev.off ()
}





## 3.2. plots to explore OUTPUT data: LOR vs. FDR (volcano plot) 
######################################################

# This plot shows final results for each method
# We want to know the relationship between combined LOR and its SE

# Select results for a specific method:
res.method = DL #  results from Dersimonian-Laird method

summary(res.method[, "summary_LOR"])
summary(res.method[, "p.adjust"])
log10 <- -log(res.method[, "p.adjust"], base = 10)
summary(log10)
corte.logfdr <- 3
table(log10 > 3)
head(res.method)

#I assign colors but the plot show other colors by default. I have to adjust it!

res.method$threshold = rep("white", nrow(res.method))
sig.up   <- (res.method[, "p.adjust"] < 0.05) & (res.method[, "summary_LOR"] > 0)
sig.down <- (res.method[, "p.adjust"] < 0.05) & (res.method[, "summary_LOR"] < 0)

res.method$threshold [sig.up]    <- "blue"
res.method$threshold [sig.down]  <- "red"
table(res.method$threshold)
res.method2 <- res.method[log10 < corte.logfdr,]
table(res.method2$threshold)

dim(res.method2)
res.method$threshold <- as.factor(res.method$threshold)

##Construct the plot object
x.por <- 3
y.por <- 3
jpeg (filename = "volcano_reactome.jpeg",  width = 480 *x.por ,   height = 480*y.por, pointsize = 100, quality = 100) 
g <- ggplot(data=res.method2, aes(x=summary_LOR, y=-log10(p.adjust), colour= threshold)) #generate data
g <- g + geom_point(alpha=0.8, size=4) 
g <- g + xlim(c(-0.51, 0.51)) + ylim(c(0, 3.1))
g <- g + xlab("log2 Odds Ratio") 
g <- g + theme(axis.title.x= element_text(colour= "black", size = 35 )) #to change size of title in axis x
g <- g + ylab("-log10 FDR") 
g <- g + theme(axis.title.y= element_text(colour= "black", size = 35 )) #to change size of title in axis y
g <- g + labs(title = "Volcano plot")
g <- g + theme(title= element_text(colour= "black", size = 35 )) #to change size of title 
g <- g + theme(axis.text.x = element_text(colour= "black", size = 35 )) #to change size of numbers
g <- g + theme(axis.text.y = element_text(colour= "black", size = 35 )) #to change size of numbers
g <- g + geom_vline(xintercept= c(0), colour = "black", size = 0.5)
g <- g + geom_hline(yintercept= -log(0.05, base= 10), colour = "black", size = 0.5)
g <- g + theme(legend.position = "none") 
print(g)
dev.off ()


## 3.3. plot to evaluate heterogeneity for all models 
######################################################

# methods<-  c(rep("DL", nrow(DL)),rep("HE", nrow(HE)),rep("HS", nrow(HS)),rep("SJ", nrow(SJ)), 
#             rep("ML", nrow(ML)),rep("REML", nrow(REML)),rep("EB", nrow(EB)),
#             rep("PM", nrow(PM)),rep("FE", nrow(FE)))
methods<-  c(rep("DL", nrow(DL)),rep("HE", nrow(HE)),rep("HS", nrow(HS)),rep("SJ", nrow(SJ)), 
             rep("FE", nrow(FE)) )
# QEp   <- c(DL$QEp, HE$QEp, HS$QEp, SJ$QEp, ML$QEp, REML$QEp, EB$QEp, PM$QEp, FE$QEp)
QEp      <- c(DL$QEp, HE$QEp, HS$QEp, SJ$QEp,                    FE$QEp)
SE       <- c(DL$SE, HE$SE, HS$SE, SJ$SE,                     FE$SE)
I2       <- c(DL$I2, HE$I2, HS$I2, SJ$I2,                     FE$I2)
H2       <- c(DL$H2, HE$H2, HS$H2, SJ$H2,                     FE$H2)
tau2     <- c(DL$tau2, HE$tau2, HS$tau2, SJ$tau2,                   FE$tau2)

datos <- data.frame (cbind(as.factor(methods), QEp, SE, I2, H2, tau2))
head(datos)
dim(datos)
summary(datos)
param <- c("QEp", "SE", "I2", "H2", "tau2")
x.por <- 3
y.por <- 2


# We drawn separately QEp because we need abline for 0.05
jpeg (filename = "reactome_het_QEp.jpeg",  
      width = 480 *x.por ,   height = 480*y.por,
      pointsize = 150, quality = 100) 
p <- ggplot(datos, aes(methods, datos[, "QEp"]))
p <- p + geom_boxplot(col = "blue")
p <- p + xlab("Métodos") 
p <- p + ylab("QEp")
p <- p + theme(axis.title.x= element_text(colour= "black", size = 40 )) #to change size of title in axis x
p <- p + theme(axis.title.y= element_text(colour= "black", size = 40 )) #to change size of title in axis y
p <- p + theme(axis.text.x = element_text(colour= "black", size = 30 )) #to change size of numbers
p <- p + theme(axis.text.y = element_text(colour= "black", size = 30 )) #to change size of numbers
p <- p + geom_hline(yintercept = 0.05, col = "red")
print(p)
dev.off ()

#Now the rest:
for (i in 2:length(param)){
  jpeg (filename = paste("reactome_het_", param[i], ".jpeg", sep = ""),  
        width = 480 *x.por ,   height = 480*y.por,
        pointsize = 150, quality = 100) 
  p <- ggplot(datos, aes(methods, datos[, param[i]]))
  p <- p + geom_boxplot(col = "blue")
  p <- p + xlab("Métodos") 
  p <- p + ylab(param[i])
  p <- p + theme(axis.title.x= element_text(colour= "black", size = 40 )) #to change size of title in axis x
  p <- p + theme(axis.title.y= element_text(colour= "black", size = 40 )) #to change size of title in axis y
  p <- p + theme(axis.text.x = element_text(colour= "black", size = 30 )) #to change size of numbers
  p <- p + theme(axis.text.y = element_text(colour= "black", size = 30 )) #to change size of numbers
  #   p <- p + geom_hline(yintercept = 0.05, col = "red")
  print(p)
  dev.off ()
}







## 3.4. plot to evaluate heterogeneity and biases for each specific model and function
###################################################### 

# Select results for a specific method:
res.method = DLs # significant results from Dersimonian-Laird method
metodo <- "DL"   # method to estimate the variability

sig.fun <- rownames(res.method)
x.por = 2.5; y.por = 2.5


if (length(sig.fun) == 0){print ("Not significant results")} else {
  for (i in 1:length(sig.fun)){
    #fit the model. It will be used for several plots
    res <- rma(yi= mat.lor[sig.fun[i],], sei =mat.sd[sig.fun[i],], method = metodo)
    
    ## A. FOREST PLOT (detailed infomation of size effect from each study)
    #   jpeg (filename = paste("reactome_",sig.fun[i],".png",sep =""),
    #         width = 480 *x.por ,   height = 480*y.por, pointsize = 100, quality = 100) 
    
    png (filename = paste("reactome_forest_", sig.fun[i],".png",sep =""),  
         width = 480 * x.por, height = 480 * y.por, res = 200)
    par(mar=c(4,4,1,2))
    forest(res, slab = toupper(substr(colnames(mat.lor),1,4)),
           xlab="Odds Ratio", cex=0.7,
           mlab=paste(metodo, "Model for All Studies", sep = " "), col = "red", 
           main = paste("\n", sig.fun[i], sep=""))    
    text( 9,-3, "Odds Ratio [IC 95%]", pos=2, cex = 0.7)
    dev.off()
    
    ## B. FUNNEL  PLOT (detailed infomation of size effect from each study)
    # A funnel plot shows the observed effect sizes or outcomes on the x-axis against 
    # some measure of precision of the observed effect sizes or outcomes on the y-axis. 
    # Based on Sterne and Egger (2001)
    # http://www.metafor-project.org/doku.php/plots:funnel_plot_variations
    png (filename = paste("REACTOME_funnel_",sig.fun[i],".png",sep =""),  
         width = 480 * x.por, height = 480 * y.por, res = 200)
    ### set up 2x2 array for plotting
    par(mfrow=c(2,2))
    ### draw funnel plots
    funnel(res, main="Standard Error", back ="darkslategray1",
           xlab = paste("LOR (", sig.fun[i], ")",sep =""))
    funnel(res, yaxis="vi", main="Sampling Variance", back ="darkslategray1",
           xlab = paste("LOR (", sig.fun[i], ")",sep =""))
    funnel(res, yaxis="seinv", main="Inverse Standard Error", back ="darkslategray1",
           xlab = paste("LOR (", sig.fun[i], ")",sep =""))
    funnel(res, yaxis="vinv", main="Inverse Sampling Variance", back ="darkslategray1",
           xlab = paste("LOR (", sig.fun[i], ")",sep =""))
    dev.off()
    
    #     Some problems with this plot!!
    #     ## C. FUNNEL PLOT with TRIM and FILL
    #     # This funnel plot helps to detect asymmetries
    #     # http://www.metafor-project.org/doku.php/plots:funnel_plot_with_trim_and_fill
    #     png (filename = paste("go_mf_funnel_tf_",sig.fun[i],".png",sep =""),
    #          width = 480 * x.por, height = 480 * y.por, res = 200)
    #     par(mar=c(5,4,1,2))
    #     ### carry out trim-and-fill analysis
    #     taf <- trimfill(res)
    #     ### draw funnel plot with missing studies filled in
    #     funnel(taf, back ="darkslategray1",xlab = paste("LOR (", sig.fun[i], ")",sep =""))
    #     dev.off()
    
    ## C. RADIAL PLOT 
    # Radial plots were introduced by Rex Galbraith (1988a, 1988b, 1994) 
    # and can be useful in the meta-analytic context to examine the data for heterogeneity.
    # http://www.metafor-project.org/doku.php/plots:radial_plot
    
    ### to save as png file
    png (filename = paste("reactome_radial_",sig.fun[i],".png",sep =""),  
         width = 480 * x.por, height = 480 * y.por, res = 200)
    ### adjust margins so the space is better used
    par(mar=c(5,4,0,2))
    ### draw radial plot
    radial(res,back ="darkslategray1",
           main = paste("\n", sig.fun[i], " (", res.method[i,"name"], ")",sep=""))
    dev.off()
    
    
    ## D. INFLUENCE PLOT 
    # That shows various diagnostic measures
    
    ### to save as png file
#     jpeg (filename = paste("reactome_influ_",sig.fun[i],".jpeg",sep =""),  
#           width = 480 *x.por ,   height = 480*y.por,
#           pointsize = 50, quality = 100) 
    png (filename = paste("REACTOME_influ_",sig.fun[i],".png",sep =""),  
         width = 480 * x.por, height = 480 * y.por, res = 200)
    inf <- influence(res)
    plot(inf, plotfb = T)
    dev.off()
  }
}





## 3.5. sensitivity analysis for a specific model 
###################################################### 

# First we calculate sd of "estimates" for each method:

methods = c("DL", "HE", "HS", "SJ",  "FE")  

mat.sen <- as.data.frame(matrix (NA, nrow = nrow(mat.lor), ncol = length (methods)))
colnames(mat.sen) <- methods
dim(mat.sen)
head(mat.sen)


for (i in methods){
  meta_analisis <- lapply(1:length(rownames(mat.lor)),
                          function(x){yi = rma(mat.lor[x, ], sei =mat.sd[x, ],
                                               method = i)})
  for (j in 1:length(meta_analisis)){
    mat.sen[j,i]  <- sd(leave1out(meta_analisis[[j]])$se)
    }
  print(i)
}

summary(mat.sen)
dim(mat.sen)
boxplot(mat.sen)
#only random effects model:
mat.sen <- mat.sen[, c("DL", "HE", "HS", "SJ" )]

gmat.sen  <- mat.sen
g2mat.sen <- melt(gmat.sen)
head(g2mat.sen)
colnames(g2mat.sen) <- c("metodo", "ee_lor")
head(g2mat.sen); dim(g2mat.sen)

setwd(file.path(.job$dir$code, "results",  "reactome",  "plots"))
x.por <- 3
y.por <- 2

jpeg (filename = "reactome_sens_allmethods.jpeg",  
      width = 480 *x.por ,   height = 480*y.por,
      pointsize = 150, quality = 100) 
p <- ggplot(g2mat.sen, aes(metodo, g2mat.sen[, "ee_lor"]))
p <- p + geom_boxplot(col = "blue")
p <- p + xlab("Método") 
p <- p + ylab("Error estándar del LOR")
p <- p + theme(axis.title.x= element_text(colour= "black", size = 40 )) #to change size of title in axis x
p <- p + theme(axis.title.y= element_text(colour= "black", size = 40 )) #to change size of title in axis y
p <- p + theme(axis.text.x = element_text(colour= "black", size = 20 )) #to change size of numbers
p <- p + theme(axis.text.y = element_text(colour= "black", size = 20 )) #to change size of numbers
print(p)
dev.off ()









#selecting one method:
res.method = DLs # significant results from Dersimonian-Laird method
metodo <- "DL"   # method to estimate the variability
sig.fun <- rownames(res.method)



setwd(file.path(.job$dir$code, "results",  "reactome",  "files"))   

if (length(sig.fun) == 0){print ("Not significant results")} else {
  for (i in 1:length(sig.fun)){
    #fit the model. It will be used for several plots
    res <- rma(yi= mat.lor[sig.fun[i],], sei =mat.sd[sig.fun[i],], method = metodo)
    res.l1out  <- leave1out(meta_analisis[[i]])
    resu.l1out <- round(print.list.rma(res.l1out),3)
    studies <- toupper(substr(colnames(mat.lor),1,4))
    resu.l1out <- cbind(studies,resu.l1out)
    #head(resu.l1out)
    write.table(resu.l1out, 
                paste("sensi_reactome_", sig.fun[i],".txt", sep =""),
                quote = F, row.names = F)
  }
}
    
  




### EXIT
warnings ()
sessionInfo ()
q ("no")


