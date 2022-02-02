##################################################
################### Test cases ###################
##################################################
rm(list = ls())
library("GEOquery")
library("affy")
​
accessionID="GSE17708"
geo=getGEO(accessionID)[[1]]
pd=pData(geo)
filePaths = getGEOSuppFiles(accessionID)
system(paste("tar xvf", sprintf("%s/%s_RAW.tar", accessionID, accessionID)))
dat=ReadAffy(filenames=basename(as.character(pd$supplementary_file)),phenoData=pd)
eset = rma(dat)
edata = exprs(eset)
​
sample_names = c (paste("Un",seq(1, 3),sep=""),
                  paste("0.5h",seq(1, 3),sep=""),
                  paste("1h",seq(1, 3),sep=""),
                  paste("2h",seq(1, 2),sep=""),
                  paste("4h",seq(1, 3),sep=""),
                  paste("8h",seq(1, 3),sep=""),
                  paste("16h",seq(1, 3),sep=""),
                  paste("24h",seq(1, 3),sep=""),
                  paste("72h",seq(1, 3),sep=""))
colnames(edata) = sample_names
​
boxplot(edata)
edata = as.data.frame(edata)
​
library("hgu133plus2.db")
x = hgu133plus2SYMBOL
# Get the probe identifiers that are mapped to a gene symbol
mapped_probes = mappedkeys(x)
# Convert to a list
xx = as.list(x[mapped_probes])
mapped_genes = lapply(rownames(edata), function(x) xx[[x]])
mapped_genes[sapply(mapped_genes, is.null)] = NA
mapped_genes = unlist(mapped_genes)
edata$Symbol = mapped_genes
​
edata = edata[!is.na(edata[,"Symbol"]),]   # remove unannotated genes
edata = aggregate(. ~ Symbol, edata, mean)  # mean expression for redudant probes
rownames(edata) = edata[,"Symbol"]
edata$Symbol = NULL
​
######### ######### 
# three groups
celltype = c (rep("Early", each=9),
              rep("Middle", each=8),
              rep("Late", each=9))
names(celltype) = sample_names
compare_list = c("Early-Middle", "Middle-Late", "Early-Late")
​
phenoData = new("AnnotatedDataFrame", data = data.frame(celltype = celltype))
rownames(phenoData) = colnames(edata)
neweset = ExpressionSet(assayData = as.matrix(edata), phenoData = phenoData)
​
save(neweset, phenoData, file = "EMT.Rdata")
​
library(NetAct)
data("hDB")
DErslt = MultiMicroDegs(neweset, compare_list)
​
#gsearslt_1 = TF_GSEA(hDB_v2, DErslt$`Early-Middle`, minSize=8, nperm = 10000, qval = T)
##gsearslt_3 = TF_GSEA(hDB_v2, DErslt$`Early-Late`, minSize=8, nperm = 10000, qval = T)
#write.csv(gsearslt_1, "results_v3_1.csv")
#write.csv(gsearslt_2, "results_v3_2.csv")
#write.csv(gsearslt_3, "results_v3_3.csv")
​
gsearslt_1 = read.csv("/Users/c-clausb/Downloads/NetAct-master 2/results_v3_1.csv")
gsearslt_2 = read.csv("/Users/c-clausb/Downloads/NetAct-master 2/results_v3_2.csv")
gsearslt_3 = read.csv("/Users/c-clausb/Downloads/NetAct-master 2/results_v3_3.csv")
gsearslts = list(gsearslt_1, gsearslt_2, gsearslt_3)
​
tfs = TF_Selection(gsearslts, qval = 0.01, compList = compare_list, ntop = NULL)
​
#tfs_1 = as.character(gsearslt_1$tf[gsearslt_1$qvals <0.05])
#tfs_2 = as.character(gsearslt_2$tf[gsearslt_2$qvals <0.05])
#tfs_3 = as.character(gsearslt_3$tf[gsearslt_3$qvals <0.05])
#tfs = sort(unique(as.character(c(tfs_1, tfs_2, tfs_3))))
#lengths(list(tfs_1,tfs_2,tfs_3, tfs))
​
​
act.me <- TF_Activity(tfs, hDB_v2, neweset, DErslt$Overall)
act.me$all_list
acts_mat = act.me$all_activities
#exprs_mat = exprs(neweset)[rownames(acts_mat), ]
Combine_heatmap(acts_mat, neweset)
​
tf_links = TF_Filter(acts_mat, hDB_v2, miTh = .05, nbins = 8, corMethod = "spearman", DPI = T)
​
#cor(acts_mat["CTNNB1",], acts_mat["ZEB1",], method = "spearman")
dim(tf_links)
length(unique(c(tf_links$Source, tf_links$Target)))
unique(c(tf_links$Source, tf_links$Target))
"SNAI1" %in% unique(c(tf_links$Source, tf_links$Target))
plot_network_v(tf_links)
save(tf_links, file = "/Users/c-clausb/Desktop/Netact/tf_links.emt")
​
​
############ switch signs ###############
​
nodes <- unique(c(tf_links$Source, tf_links$Target))
​
#idenify DEs and non DEs
non.DEGs <- row.names(acts_mat)[row.names(acts_mat) %in% row.names(DErslt$Overall$table)[DErslt$Overall$table$adj.P.Val >= .05]]
non.in.net <- non.DEGs[non.DEGs %in% nodes]
DEGs <- row.names(acts_mat)[row.names(acts_mat) %in% row.names(DErslt$Overall$table)[DErslt$Overall$table$adj.P.Val < .05]]
de.in.net <- DEGs[DEGs %in% nodes]
​
#identify DE and non DE targets
de.in.net.targets <- vector(mode = "list", length = length(de.in.net))
for (i in 1:length(de.in.net)){
  de.in.net.targets[[i]] <- row.names(act.me$all_list[de.in.net[[i]]][[1]])
}
names(de.in.net.targets) <- de.in.net
#View(de.in.net.targets)
​
non.in.net.targets <- vector(mode = "list", length = length(non.in.net))
for (i in 1:length(non.in.net)){
  non.in.net.targets[[i]] <- row.names(act.me$all_list[non.in.net[[i]]][[1]])
}
names(non.in.net.targets) <- non.in.net
#View(non.in.net.targets)
​
​
ex.genes <- row.names(DErslt$Overall$table$AveExpr)
#### now need to find overlap between targets
​
#make pairwise combos of two lists
a <- length(non.in.net)
b <- length(de.in.net)
C <- expand.grid(1:a, 1:b)  
combos.vec <- vector(length = nrow(C))
o.mat <- as.data.frame(matrix(nrow = length(combos.vec), ncol = 10))
colnames(o.mat) <- c("pair", "tf1.targets", "tf2.targets", "tf1.non", "tf2.non", "A","B","C","D", "P.value")
​
for (i in 1:nrow(C)){
  combos.vec[[i]] <-paste(C[i,][[1]], C[i,][[2]], sep = " ")
}
​
P.mat <- matrix(nrow = length(non.in.net), ncol = length(de.in.net))
rownames(P.mat) <- non.in.net
colnames(P.mat) <- de.in.net
ex.genes <- names(act.me$all_list)
​
for (i in 1:length(combos.vec)){
  tf1 <- non.in.net[[as.integer(strsplit(combos.vec[[i]], split = " ")[[1]][1])]]
  tf2 <-de.in.net[[as.integer(strsplit(combos.vec[[i]], split = " ")[[1]][2])]]
  o.mat[i,1] <- paste0(tf1, " and ", tf2)
  tf1.targets <- non.in.net.targets[[tf1]]
  tf2.targets <- de.in.net.targets[[tf2]]
  tf1.non <- setdiff(ex.genes, tf1.targets)
  tf2.non <- setdiff(ex.genes, tf2.targets)
  A <- length(intersect(tf1.targets, tf2.targets))
  B <- length(intersect(tf1.targets, tf2.non))
  C <- length(intersect(tf1.non, tf2.targets))
  D <- length(intersect(tf1.non, tf2.non))
  fish.mat <- matrix(c(A,C,B,D), nrow = 2, ncol=2)
  #fish <- fisher.test(fish.mat)
  fish <- fisher.test(fish.mat, alternative = "greater")
  #fish <- fisher.test(fish.mat, alternative = "less")
  o.mat[i,2] <- length(tf1.targets)
  o.mat[i,3] <- length(tf2.targets)
  o.mat[i,4] <- length(tf1.non)
  o.mat[i,5] <- length(tf2.non)
  o.mat[i,6] <- A
  o.mat[i,7] <- B
  o.mat[i,8] <- C
  o.mat[i,9] <- D
  o.mat[i,10] <- fish$p.value
  P.mat[tf1,tf2] <- fish$p.value
}
#View(P.mat)
#im(P.mat)
# for each non DE find top 5 most signifcant P values
# first make a matrix where each row is a non DE and each column is a DE it pairs with and the value is the p value
​
P.mat <- t(P.mat)
P.mat <- as.data.frame(P.mat)
​
fish.list <- vector(mode = "list")
for(i in 1:length(colnames(P.mat))){
  tmp <-P.mat[,i]
  names(tmp) <-rownames(P.mat)
  tmp <- tmp[tmp < .05]
  fish.list[[i]] <-names(tmp[order(tmp)])
}
length(fish.list)
names(fish.list) <- colnames(P.mat)
View(fish.list)
​
setdiff(names(fish.list),non.in.net)
lengths(fish.list)
names(fish.list)
​
fish.list$BRCA1
fish.list$CREBBP
fish.list$E2F4
fish.list$ESR1
fish.list$HOXB7
fish.list$IRF9
fish.list$KLF4
fish.list$MYC
fish.list$PAX6
fish.list$RB1
fish.list$RELA
fish.list$SMAD2
fish.list$SMAD4
fish.list$SMAD7
fish.list$SNAI1
fish.list$SP1
fish.list$STAT3
fish.list$TFAP2A
fish.list$TP53
fish.list$TWIST1
fish.list$WT1
​
############ switch signs based on literature evidence ##########
#switch signs that DONT agree 
NDE <- "CREBBP"
DE <-"CTNNB1"
sign(acts_mat[NDE,])
sign(acts_mat[DE,])
acts_mat[NDE,][sign(acts_mat[NDE,]) != sign(acts_mat[DE,])] <- -(acts_mat[NDE,][sign(acts_mat[NDE,]) != sign(acts_mat[DE,])])
​
​
#switch signs that DO agree 
NDE <- "E2F4"
DE <-"SMAD3"
sign(acts_mat[NDE,])
sign(acts_mat[DE,])
acts_mat[NDE,][sign(acts_mat[NDE,]) == sign(acts_mat[DE,])] <- -(acts_mat[NDE,][sign(acts_mat[NDE,]) == sign(acts_mat[DE,])])
​
​
#switch signs that DONT agree 
NDE <- "ESR1"
DE <-"JUN"
sign(acts_mat[NDE,])
sign(acts_mat[DE,])
acts_mat[NDE,][sign(acts_mat[NDE,]) != sign(acts_mat[DE,])] <- -(acts_mat[NDE,][sign(acts_mat[NDE,]) != sign(acts_mat[DE,])])
​
#switch signs that DONT agree 
NDE <- "IRF9"
DE <-"STAT1"
sign(acts_mat[NDE,])
sign(acts_mat[DE,])
acts_mat[NDE,][sign(acts_mat[NDE,]) != sign(acts_mat[DE,])] <- -(acts_mat[NDE,][sign(acts_mat[NDE,]) != sign(acts_mat[DE,])])
​
​
#switch signs that DO agree 
NDE <- "KLF4"
DE <-"CTNNB1"
sign(acts_mat[NDE,])
sign(acts_mat[DE,])
acts_mat[NDE,][sign(acts_mat[NDE,]) == sign(acts_mat[DE,])] <- -(acts_mat[NDE,][sign(acts_mat[NDE,]) == sign(acts_mat[DE,])])
​
​
#switch signs that DO agree 
NDE <- "RB1"
DE <-"E2F1"
sign(acts_mat[NDE,])
sign(acts_mat[DE,])
acts_mat[NDE,][sign(acts_mat[NDE,]) == sign(acts_mat[DE,])] <- -(acts_mat[NDE,][sign(acts_mat[NDE,]) == sign(acts_mat[DE,])])
​
#switch signs that DONT agree 
NDE <- "SMAD4"
DE <-"SMAD3"
sign(acts_mat[NDE,])
sign(acts_mat[DE,])
acts_mat[NDE,][sign(acts_mat[NDE,]) != sign(acts_mat[DE,])] <- -(acts_mat[NDE,][sign(acts_mat[NDE,]) != sign(acts_mat[DE,])])
​
#switch signs that DONT agree 
NDE <- "SP1"
DE <-"JUN"
sign(acts_mat[NDE,])
sign(acts_mat[DE,])
acts_mat[NDE,][sign(acts_mat[NDE,]) != sign(acts_mat[DE,])] <- -(acts_mat[NDE,][sign(acts_mat[NDE,]) != sign(acts_mat[DE,])])
​
#switch signs that DO agree 
NDE <- "STAT3"
DE <-"STAT1"
sign(acts_mat[NDE,])
sign(acts_mat[DE,])
acts_mat[NDE,][sign(acts_mat[NDE,]) == sign(acts_mat[DE,])] <- -(acts_mat[NDE,][sign(acts_mat[NDE,]) == sign(acts_mat[DE,])])
​
#switch signs that DONT agree 
NDE <- "TP53"
DE <-"HIF1A"
sign(acts_mat[NDE,])
sign(acts_mat[DE,])
acts_mat[NDE,][sign(acts_mat[NDE,]) != sign(acts_mat[DE,])] <- -(acts_mat[NDE,][sign(acts_mat[NDE,]) != sign(acts_mat[DE,])])
​
​
############### now to remake network and remove nodes ##############
​
Combine_heatmap(acts_mat, neweset)
​
tf_links = TF_Filter(acts_mat, hDB_v2, miTh = .05, nbins = 8, corMethod = "spearman", DPI = T)
​
​
​
#,"HOXB7", "MYC", "PAX6", "SMAD2", "SMAD7", "SNAI1", "TFAP2A", "TWIST1", "WT1")
​
tf.re <- "BRCA1"
tf_links <-tf_links[tf_links$Source != tf.re,]
tf_links <-tf_links[tf_links$Target != tf.re,]
​
tf.re <- "HOXB7"
tf_links <-tf_links[tf_links$Source != tf.re,]
tf_links <-tf_links[tf_links$Target != tf.re,]
​
tf.re <- "MYC"
tf_links <-tf_links[tf_links$Source != tf.re,]
tf_links <-tf_links[tf_links$Target != tf.re,]
​
​
tf.re <- "PAX6"
tf_links <-tf_links[tf_links$Source != tf.re,]
tf_links <-tf_links[tf_links$Target != tf.re,]
​
tf.re <- "SMAD2"
tf_links <-tf_links[tf_links$Source != tf.re,]
tf_links <-tf_links[tf_links$Target != tf.re,]
​
tf.re <- "SMAD7"
tf_links <-tf_links[tf_links$Source != tf.re,]
tf_links <-tf_links[tf_links$Target != tf.re,]
​
tf.re <- "SNAI1"
tf_links <-tf_links[tf_links$Source != tf.re,]
tf_links <-tf_links[tf_links$Target != tf.re,]
​
​
tf.re <- "TFAP2A"
tf_links <-tf_links[tf_links$Source != tf.re,]
tf_links <-tf_links[tf_links$Target != tf.re,]
