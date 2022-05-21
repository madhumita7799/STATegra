library(NOISeq)
library(sva)
library(biomaRt)
library(edgeR)
library(cqn)
library(scales)
library(limma)
library(maSigPro)
library(readxl)

#################################################
## for miRNA
####################################################

### loading data
mydata0 = read.csv("STATegra.miRNAseq.allSamples.counts.csv",sep="\t")
des.mat = read.csv("miRNAseq_ExpDesign.txt", row.names = 1)

# low count filtering 
## miRNAs with low average counts in all conditions are removed With CPM method from NOISeq package

Fdata = filtered.data(mydata0, 
                      factor = apply(des.mat[,c("cond","time")], 1, paste, collapse = ""),
                      norm = FALSE, method = 1, cv.cutoff = 500, cpm = 1)
## Data oredering
DESIGN.ORD = des.mat[order(des.mat$cond,as.numeric(des.mat$time),des.mat$repB,des.mat$batch),]
Fdata.ORD = Fdata[,rownames(DESIGN.ORD)]

# Defining factors for batch effect correction
DESIGN.ORD$time = as.numeric(DESIGN.ORD$time)
Time = DESIGN.ORD$time
IKAROS = DESIGN.ORD$cond
BATCH = DESIGN.ORD$batch  
ik.con = paste(Time,"h",IKAROS,sep="")

### visualizing the data before normalization
count_mirna = as.matrix(Fdata.ORD)
par(mfrow = c(1,1))
boxplot(count_mirna~col(count_mirna), names = colnames(count_mirna), 
        col = rep(c(2:7), each = 3), las = 2, main = "miRNA_expression", cex.axis = 0.7)

# TMM normalization:Normalizing data to remove bias related with different count distribution across samples
datosnorm = tmm(Fdata.ORD)

# Batch effect correction

design.cb = model.matrix(~ ik.con)

datosnorm = ComBat(log2(datosnorm+1), batch=BATCH, mod=design.cb, par.prior = FALSE) 

# Because of negative values after the batch correction, a constant is added to the data
min(datosnorm)
datosCONrepl = datosnorm + 0.9

# The technical replicates are avaraged

condrep = gsub(" ", "", apply(DESIGN.ORD[,c("cond", "time", "repB")], 1, paste, collapse = "-"))
condrepunique = unique(condrep)

numtec = table(condrep)

sintec = NULL
for (ccc in condrepunique) {  
  if (numtec[ccc] == 1) {
    sintec = cbind(sintec, datosCONrepl[,grep(ccc, condrep)])
  } else {
    sintec = cbind(sintec, rowMeans(datosCONrepl[,grep(ccc, condrep)]))
  }  
}

colnames(sintec) = condrepunique

##### visualization
par(mfrow= c(1,1))
boxplot(sintec~ col(sintec), names = colnames(sintec), 
        col = rep(c(2:7), each = 3), las = 2, main = "miRNA_expression", cex.axis = 0.8)
par(mfrow= c(1,3))
plotIndiv(pca(t(sintec), center = TRUE, scale = TRUE), col = mycol, style = "graphics" , 
          title = "miRNA_Expression", size.title = 1)
plotIndiv(pca(t(sintec[,1:18]), center = TRUE, scale = TRUE), 
          col = mycol, style = "graphics" , title = "miRNA_expression (Control)", size.title = 1)
plotIndiv(pca(t(sintec[,19:36]), center = TRUE, scale = TRUE), 
          col = mycol, style = "graphics" , title = "miRNA_expression (IKAROS)", size.title = 1)

DESIGN.notec = t(as.data.frame(strsplit(condrepunique, "-")))
rownames(DESIGN.notec) = condrepunique
colnames(DESIGN.notec) = c("cond", "time", "repB")

##### saving the final matrix for downstream analysis
write.table(sintec, "miRNA_final.txt", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE) 
write.table(DESIGN.notec, "miRNAseq_ExpDesign.txt", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

####################################################################
###### for mRNA
###############################################################
M.intersect.36 = read.csv("RNA_seq/STATegra.RNAseq.allSamples.counts.csv",row.names=1)
DESIGN.36 = read.csv("RNA_seq/DESIGN.36.csv",row.names = 1)
rownames(DESIGN.36) == colnames(M.intersect.36) ### check gene organization

## Annotation downloaded from BiomaRt 
gene.INFO = read.csv("RNA_seq/ENS_MOUSE_LENGHT_GC.txt",sep="\t",row.names=1)


### Low count gene filtration (Decision criteria:at least 10 reads in at least 20 samples)
M.intersect.n10 = apply(M.intersect.36,1,function(x){sum(x>10)})
M.intersect.36F = M.intersect.36[M.intersect.n10>=20,]
M.go = M.intersect.36F
M.go = M.go[rownames(M.go) %in% rownames(gene.INFO),] 
gene.INFO.GO = gene.INFO[rownames(gene.INFO) %in% rownames(M.go),]

## Data oredering
M.go.ORD=M.go[order(rownames(M.go)),]
gene.INFO.GO.ORD=gene.INFO.GO[order(rownames(gene.INFO.GO)),]
sizeFactors.GO=apply(M.go.ORD,2,sum)

DESIGN.ORD=DESIGN.36[colnames(M.go.ORD),]
rownames(DESIGN.ORD)==colnames(M.go.ORD)

Time = DESIGN.ORD$TIME
IKAROS = DESIGN.ORD$IKAROS
BATCH = DESIGN.ORD$EXP_BATCH
LIB = DESIGN.ORD$LIBRARY_PREP

## 12 groups: 6 times x 2 conditions are considered
TIME_IK = paste("T",Time,"_","Exp",IKAROS,sep="")


### visualizing the data before normalization
count = as.matrix(M.go.ORD)
par(mfrow= c(1,1))
boxplot(count ~ col(count), names = colnames(count), 
        col = rep(c(2:7), each = 3), las = 2, main = "mRNA_expression", cex.axis = 0.5)


#### Conditional quantile normalization (Cqn)
cqn.subset = cqn(M.go.ORD, lengths = gene.INFO.GO.ORD$LENGTH,
                  x = gene.INFO.GO.ORD$GC, sizeFactors = sizeFactors.GO,
                  verbose = TRUE)
M.go.ORD.cqn = cqn.subset$y + cqn.subset$offset
head(M.go.ORD.cqn)

#### Batch effect correction

design.cb = model.matrix(~ TIME_IK )
quant.batch.adj.CQN.LIB = ComBat(M.go.ORD.cqn, batch=LIB, mod=design.cb) 
n.sv.CQN = num.sv(M.go.ORD.cqn,design.cb,method="leek")
n.sv.CORRECTED.CQN = num.sv(quant.batch.adj.CQN.LIB,design.cb,method="leek")


#### Correction to set minimun value to 1 (a linear correction)
DE.set2 = quant.batch.adj.CQN.LIB+abs(min(quant.batch.adj.CQN.LIB)) + 1

##### visualization
par(mfrow = c(1,1))
boxplot(DE.set2~ col(DE.set2), names = colnames(DE.set2), 
        col = rep(c(2:7), each = 3), las = 2, main = "mRNA_expression", cex.axis = 0.5)

par(mfrow= c(1,3))
plotIndiv(pca(t(DE.set2), center = TRUE, scale = TRUE), col = mycol, style = "graphics" , 
          title = "mRNA_Expression", size.title = 1)

plotIndiv(pca(t(DE.set2[,1:18]), center = TRUE, scale = TRUE), 
          col = mycol, style = "graphics" , title = "mRNA_expression (Control)", size.title = 1)

plotIndiv(pca(t(DE.set2[,19:36]), center = TRUE, scale = TRUE), 
          col = mycol, style = "graphics" , title = "mRNA_expression (IKAROS)", size.title = 1)

#### Saving the file
write.csv(DE.set2,"mRNA_final.csv")

####################################################################
## for metabolomics
###########################################################

### defining functions
view.reporter = function (data, reporter, one.series = F, cond , mylegend = "Ik/Ctrl") {
  repo.data = data[reporter,]
  ave.data = t(apply(repo.data, 1, function (x) tapply(x,cond,mean, na.rm = T)))
  ave.data = ave.data[,c(1,4,6,2,3,5,7,10,12,8,9,11)]
  cols = colnames(ave.data)
  print(ave.data)
  time = c(0,2,5,12,18,24)
  colors = c("red", 'blue')
  for (i in 1:length(reporter)) {
    plot(time, y = rep(1, length(time)), ylim = c(min(ave.data[i,], na.rm = T), max(ave.data[i,], na.rm = T)), 
         ty = "l", ylab = "Expression", col = "white",,main = reporter[i])
    x.data = ave.data[i,1:6]
    lines(time[!is.na(x.data)], x.data[!is.na(x.data)], col = colors[1])
    if (one.series == F) {
      y.data = ave.data[i,7:12]
      lines(time[!is.na(y.data)], y.data[!is.na(y.data)], col = colors[2])
      legend("topright", legend = c( "Control", "Ikaros"), text.col = colors)
    } else {
      legend("topright", legend = mylegend)
    }
  }  
}

process = function (data, title = NULL) {
  batch = rep(c(9,10,11), each = 12)
  time = rep(rep(c(0,2,6,12,18,24), each = 2), 3)
  treatment = rep (c("Control", "Ikaros"), 18)
  data = data[,-1]
  rownames(data) = paste(treatment, "_", time, "_h_batch_", batch, sep = "")
  data = t(data)
  data = data[,sort(colnames(data))]
  reorder = c(1:3, 13:18, 4:12, 19:21, 31:36,22:30)
  data = data[,reorder]
  boxplot(log(data)~ col(data), names = colnames(data), 
          col = rep(c(2:7), each = 3), las = 2, main = title, cex.axis = 0.6)
  data
}

visualization =function (data, reporter, cond, main, ikaros = c(19:36)) {
  par(mfrow = c(2,3))
  mycol = rep(rainbow(6), each = 3)
  mycol[mycol == "#FFFF00FF"] = "orange"
  plotIndiv(pca(t(data), center = TRUE, scale = TRUE), col = mycol, style = "graphics" , title = main, size.title = 1)
  plotIndiv(pca(t(data[,ikaros]), center = TRUE, scale = TRUE), col = mycol, style = "graphics" , title = main, size.title = 1)
  boxplot(data ~ col(data), names = colnames(data), 
          col = mycol, las = 2, main = main, cex.axis = 0.6)
  view.reporter (data, reporter, cond = index)
  
}

##### Data loading
LC = as.data.frame(read_excel("LC_Processed.xls", sheet = "Study Samples", range = "A1:AK37"))
GC = as.data.frame(read_excel("GC_Processed.xls", sheet = "Study Samples", range = "A1:W37"))

###### Data ordering
par(mfrow = c(2,1))
LC.data =process(data = LC, title = "LC data")
GC.data=process(data = GC, title = "GC data")
dev.off()

####### Fusing data from the two platforms

fused.data = rbind(LC.data, GC.data)
fused.data.log = log(fused.data)


####### mean centering
index = paste(rep(c("Control", "Treatment"), each = 18), rep(c(0,2,6,12,18,24), each = 3), sep= "_") 
median = apply(fused.data.log, 2, median, na.rm = T) # median per sample
means.of.medians = tapply(median, index, mean)[c(1,4,6,2,3,5,7,10,12,8,9,11)]
means.of.medians.vector = rep(means.of.medians, each = 3)
means.of.medians.deviation = median-means.of.medians.vector
means.of.medians.deviation.matrix = matrix(data = rep(means.of.medians.deviation, nrow(fused.data.log)), 
                                            nrow = nrow(fused.data.log), ncol = ncol(fused.data.log), byrow = T)
fused.data.log.mean = fused.data.log - means.of.medians.deviation.matrix

## Futher visualisation and check with reporter metabolites
reporter = c("lacticacid", "glucose", "pyruvicacid")
visualization (fused.data, reporter, cond = index, main = "Fused data")
visualization (fused.data.log, reporter, cond = index, main = "Fused data.log")
visualization (fused.data.log.mean, reporter, cond = index, main = "Fused data.log.mean")

#### Saving the file
write.table (fused.data.log.mean, "Metabolomics_final.txt", quote = F, row.names = T, col.names = T, sep = "\t")
















