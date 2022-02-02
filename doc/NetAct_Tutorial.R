## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----quickYo, eval = FALSE----------------------------------------------------
#  if(!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("NetAct")
#  install_github("suke18/NetAct", dependencies=T, build_vignettes = T)

## ----load data, eval=FALSE----------------------------------------------------
#  rm(list = ls())
#  library("GEOquery")
#  library("affy")
#  accessionID="GSE17708"
#  geo=getGEO(accessionID)[[1]]
#  pd=pData(geo)
#  filePaths = getGEOSuppFiles(accessionID)
#  system(paste("tar xvf", sprintf("%s/%s_RAW.tar", accessionID, accessionID)))
#  dat=ReadAffy(filenames=basename(as.character(pd$supplementary_file)),phenoData=pd)
#  eset = rma(dat)
#  edata = exprs(eset)
#  sample_names = c (paste("Un",seq(1, 3),sep=""),
#                    paste("0.5h",seq(1, 3),sep=""),
#                    paste("1h",seq(1, 3),sep=""),
#                    paste("2h",seq(1, 2),sep=""),
#                    paste("4h",seq(1, 3),sep=""),
#                    paste("8h",seq(1, 3),sep=""),
#                    paste("16h",seq(1, 3),sep=""),
#                    paste("24h",seq(1, 3),sep=""),
#                    paste("72h",seq(1, 3),sep=""))
#  colnames(edata) = sample_names

## ----match gene names, eval=FALSE---------------------------------------------
#  boxplot(edata)
#  edata = as.data.frame(edata)
#  library("hgu133plus2.db")
#  x = hgu133plus2SYMBOL
#  # Get the probe identifiers that are mapped to a gene symbol
#  mapped_probes = mappedkeys(x)
#  # Convert to a list
#  xx = as.list(x[mapped_probes])
#  mapped_genes = lapply(rownames(edata), function(x) xx[[x]])
#  mapped_genes[sapply(mapped_genes, is.null)] = NA
#  mapped_genes = unlist(mapped_genes)
#  edata$Symbol = mapped_genes
#  edata = edata[!is.na(edata[,"Symbol"]),]   # remove unannotated genes
#  edata = aggregate(. ~ Symbol, edata, mean)  # mean expression for redudant probes
#  rownames(edata) = edata[,"Symbol"]
#  edata$Symbol = NULL

## ----store the data, eval=FALSE-----------------------------------------------
#  celltype = c (rep("Early", each=9), rep("Middle", each=8), rep("Late", each=9))
#  names(celltype) = sample_names
#  phenoData = new("AnnotatedDataFrame", data = data.frame(celltype = factor(celltype)))
#  rownames(phenoData) = colnames(edata)
#  neweset = ExpressionSet(assayData = as.matrix(edata), phenoData = phenoData)
#  save(neweset, file = "EMT.Rdata")

## ----pressure, eval=FALSE-----------------------------------------------------
#  library(NetAct)
#  data("hDB")
#  data("EMT")
#  compare_list = c("Early-Middle", "Middle-Late", "Early-Late")
#  DErslt = MicroProcess(neweset, compare_list)

## ----permutation of the GSEA, eval=FALSE--------------------------------------
#  gsearslt_1 = TF_GSEA(hDB, DErslt$`Early-Middle`, minSize=8, nperm = 10000, qval = T)
#  #gsearslt_2 = TF_GSEA(hDB, DErslt$`Middle-Late`, minSize=8, nperm = 10000, qval = T)
#  #gsearslt_3 = TF_GSEA(hDB, DErslt$`Early-Late`, minSize=8, nperm = 10000, qval = T)

## ----aggregate the gene list, eval=FALSE--------------------------------------
#  gsearslt_1 = read.csv("vignettes/results_v3_1.csv")
#  gsearslt_2 = read.csv("vignettes/results_v3_2.csv")
#  gsearslt_3 = read.csv("vignettes/results_v3_3.csv")
#  gsearslts = list(gsearslt_1, gsearslt_2, gsearslt_3)
#  tfs = TF_Selection(gsearslts, qval = 0.01, compList = compare_list, ntop = NULL)

## ----calculation activities, eval = FALSE-------------------------------------
#  act.me <- TF_Activity(tfs, hDB, neweset, DErslt$Overall)
#  act.me$all_list
#  acts_mat = act.me$all_activities
#  Combine_heatmap(acts_mat, neweset)

## ----filter links, eval=FALSE-------------------------------------------------
#  tf_links = TF_Filter(acts_mat, hDB, miTh = .05, nbins = 8, corMethod = "spearman", DPI = T)
#  cor(acts_mat["CTNNB1",], acts_mat["ZEB1",], method = "spearman")
#  dim(tf_links)
#  length(unique(c(tf_links$Source, tf_links$Target)))
#  unique(c(tf_links$Source, tf_links$Target))
#  "SNAI1" %in% unique(c(tf_links$Source, tf_links$Target))
#  plot_network_v(tf_links)

## ----switch signs, eval=FALSE-------------------------------------------------
#  nodes = unique(c(tf_links$Source, tf_links$Target))
#  #idenify DEs and non DEs
#  non.DEGs = row.names(acts_mat)[row.names(acts_mat) %in% row.names(DErslt$Overall$table)[DErslt$Overall$table$adj.P.Val >= .05]]
#  non.in.net = non.DEGs[non.DEGs %in% nodes]
#  DEGs = row.names(acts_mat)[row.names(acts_mat) %in% row.names(DErslt$Overall$table)[DErslt$Overall$table$adj.P.Val < .05]]
#  de.in.net = DEGs[DEGs %in% nodes]
#  #identify DE and non DE targets
#  de.in.net.targets = vector(mode = "list", length = length(de.in.net))
#  for (i in 1:length(de.in.net)){
#    de.in.net.targets[[i]] = names(act.me$all_list[de.in.net[[i]]][[1]])
#  }
#  names(de.in.net.targets) = de.in.net

## ----find non-DE genes targets, eval=FALSE------------------------------------
#  non.in.net.targets = vector(mode = "list", length = length(non.in.net))
#  for (i in 1:length(non.in.net)){
#    non.in.net.targets[[i]] = names(act.me$all_list[non.in.net[[i]]][[1]])
#  }
#  names(non.in.net.targets) = non.in.net

## ----create the matrix, eval=FALSE--------------------------------------------
#  a = length(non.in.net)
#  b = length(de.in.net)
#  C = expand.grid(1:a, 1:b)
#  combos.vec = vector(length = nrow(C))
#  o.mat = as.data.frame(matrix(nrow = length(combos.vec), ncol = 10))
#  colnames(o.mat) = c("pair", "tf1.targets", "tf2.targets", "tf1.non", "tf2.non", "A","B","C","D", "P.value")
#  for (i in 1:nrow(C)){
#    combos.vec[[i]] = paste(C[i,][[1]], C[i,][[2]], sep = " ")
#  }
#  P.mat = matrix(nrow = length(non.in.net), ncol = length(de.in.net))
#  rownames(P.mat) = non.in.net
#  colnames(P.mat) = de.in.net
#  ex.genes = names(act.me$all_list)

## ----build up relations, eval=FALSE-------------------------------------------
#  for (i in 1:length(combos.vec)){
#    tf1 <- non.in.net[[as.integer(strsplit(combos.vec[[i]], split = " ")[[1]][1])]]
#    tf2 <-de.in.net[[as.integer(strsplit(combos.vec[[i]], split = " ")[[1]][2])]]
#    o.mat[i,1] <- paste0(tf1, " and ", tf2)
#    tf1.targets <- non.in.net.targets[[tf1]]
#    tf2.targets <- de.in.net.targets[[tf2]]
#    tf1.non <- setdiff(ex.genes, tf1.targets)
#    tf2.non <- setdiff(ex.genes, tf2.targets)
#    A <- length(intersect(tf1.targets, tf2.targets))
#    B <- length(intersect(tf1.targets, tf2.non))
#    C <- length(intersect(tf1.non, tf2.targets))
#    D <- length(intersect(tf1.non, tf2.non))
#    fish.mat <- matrix(c(A,C,B,D), nrow = 2, ncol=2)
#    #fish <- fisher.test(fish.mat)
#    fish <- fisher.test(fish.mat, alternative = "greater")
#    #fish <- fisher.test(fish.mat, alternative = "less")
#    o.mat[i,2] <- length(tf1.targets)
#    o.mat[i,3] <- length(tf2.targets)
#    o.mat[i,4] <- length(tf1.non)
#    o.mat[i,5] <- length(tf2.non)
#    o.mat[i,6] <- A
#    o.mat[i,7] <- B
#    o.mat[i,8] <- C
#    o.mat[i,9] <- D
#    o.mat[i,10] <- fish$p.value
#    P.mat[tf1,tf2] <- fish$p.value
#  }

