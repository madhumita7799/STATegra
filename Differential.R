library(org.Mm.eg.db)
library(AnnotationDbi)
library(limma)
library(biomaRt)
library(clusterProfiler)
library(DOSE)
library(pathview)
library(enrichplot)
library(ggplot2)
library(pathview)
BiocManager::install("clusterProfiler")

################## loading data
mirna = read.table(file="miRNA/miRNA_final.txt", sep="\t", header=T)
mrna = read.table(file="mRNA/mRNA_final.csv", sep=",", header=T, row.names=1)
met = read.table(file="Metabolomics/Metabolomics_final.txt", sep="\t", header=T, row.names=1)

###### differential expression analysis
exp = mrna
name = rownames(mrna)

hh = as.data.frame(mapIds(org.Mm.eg.db, name, keytype="ENSEMBL", column="SYMBOL", multiVals = "first"))
symbol = hh$`mapIds(org.Mm.eg.db, name, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")`

#si= str_replace_all(name, "mmu-miR-", "") 

control = exp[,1:18]
ikiros = exp[,19:36]

C_0 = exp[,1:3]
I_0 = exp[,19:21]

C_2 = exp[,4:6]
I_2 = exp[,22:24]

C_6 = exp[,7:9]
I_6 = exp[,25:27]

C_12 = exp[,10:12]
I_12 = exp[,28:30]

C_18 = exp[,13:15]
I_18 = exp[,31:33]

C_24 = exp[,16:18]
I_24 = exp[,34:36]

group1 = C_0
group2 = C_24
dat = cbind(group1,group2)

#group = c(rep("1",18),rep("2",18))
group = c(rep("1",3),rep("2",3))
design = model.matrix(~0+group)

contrast.matrix1 = makeContrasts(group1-group2,levels=design)
fit1=lmFit(dat, design)
fit1=eBayes(fit1)
fit1= contrasts.fit(fit1, contrast.matrix1)
fit1= eBayes(fit1)
results1=decideTests(fit1)
#vennDiagram(results1)
probeset.list1=topTable(fit1, number= "all",adjust="BH",p.value = 0.01)
dim(probeset.list1)
#probe1=rownames(probeset.list1)
#probe1=as.data.frame(probe1)
write.table(probeset.list1,"/home/cblab/Documents/Barmingam/B3_Data/Enrichment/Diff_exp/diff_mrna_I0-24H_6391", sep="\t", quote=FALSE)

par(mfrow= c(1,1))
volcanoplot(fit1, coef = 1, highlight=10, hl.col="red",names = symbol, pch=5, cex=0.32, main="24H")

#ao_10 = fit1
#a2_215 = fit1
#a6_481 = fit1
#a12_1982 = fit1
#a18_4426 = fit1
#a24_6098 = fit1
#par(mfrow= c(1,1))

################# gene enrichment analysis

organism = "org.Mm.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

############## biological processes
log_fold = probeset.list1$logFC
names(log_fold)= rownames(probeset.list1)
gene_list = na.omit(log_fold)
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list
length(gene_list)

gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.01, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")
write.table(gse@result,"/home/cblab/Documents/Barmingam/B3_Data/Enrichment/gene/diff_I0-24_BP", sep="\t", quote=FALSE)

png(file="/home/cblab/Documents/Barmingam/plots/Enrichment/Gene/I0-24_BP_dotplot.png",width=1748, height=912)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()

png(file="/home/cblab/Documents/Barmingam/plots/Enrichment/Gene/I0-24_BP_map.png",width=994, height=912)
emapplot(gse, showCategory = 10)
dev.off()

png(file="/home/cblab/Documents/Barmingam/plots/Enrichment/Gene/I0-24_BP_ridge.png",width=1657, height=912)
ridgeplot(gse) + labs(x = "enrichment distribution")
dev.off()

png(file="/home/cblab/Documents/Barmingam/plots/Enrichment/Gene/I0-24_BP_lit.png",width=994, height=912)
terms <- gse$Description[1:3]
pmcplot(terms, 2015:2021, proportion=FALSE)
dev.off()

############# molecular function
gse <- gseGO(geneList=gene_list, 
             ont ="MF", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.01, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")
write.table(gse@result,"/home/cblab/Documents/Barmingam/B3_Data/Enrichment/gene/diff_I0-24_MF", sep="\t", quote=FALSE)

png(file="/home/cblab/Documents/Barmingam/plots/Enrichment/Gene/I0-24_MF_dotplot.png",width=1748, height=912)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()

png(file="/home/cblab/Documents/Barmingam/plots/Enrichment/Gene/I0-24_MF_map.png",width=994, height=912)
emapplot(gse, showCategory = 10)
dev.off()

png(file="/home/cblab/Documents/Barmingam/plots/Enrichment/Gene/I0-24_MF_ridge.png",width=1657, height=912)
ridgeplot(gse) + labs(x = "enrichment distribution")
dev.off()

png(file="/home/cblab/Documents/Barmingam/plots/Enrichment/Gene/I0-24_MF_lit.png",width=994, height=912)
terms <- gse$Description[1:3]
pmcplot(terms, 2015:2021, proportion=FALSE)
dev.off()


############### kegg pathway enrichment analysis
exp = mrna
name = rownames(mrna)

hh = as.data.frame(mapIds(org.Mm.eg.db, name, keytype="ENSEMBL", column="SYMBOL", multiVals = "first"))
symbol = hh$`mapIds(org.Mm.eg.db, name, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")`

control = exp[,1:18]
ikiros = exp[,19:36]

C_0 = exp[,1:3]
I_0 = exp[,19:21]

C_2 = exp[,4:6]
I_2 = exp[,22:24]

C_6 = exp[,7:9]
I_6 = exp[,25:27]

C_12 = exp[,10:12]
I_12 = exp[,28:30]

C_18 = exp[,13:15]
I_18 = exp[,31:33]

C_24 = exp[,16:18]
I_24 = exp[,34:36]

group1 = I_0
group2 = I_24
dat = cbind(group1,group2)

#group = c(rep("1",18),rep("2",18))
group = c(rep("1",3),rep("2",3))
design = model.matrix(~0+group)

contrast.matrix1 = makeContrasts(group1-group2,levels=design)
fit1=lmFit(dat, design)
fit1=eBayes(fit1)
fit1= contrasts.fit(fit1, contrast.matrix1)
fit1= eBayes(fit1)
results1=decideTests(fit1)
#vennDiagram(results1)
probeset.list1=topTable(fit1, number= "all",adjust="BH",p.value = 0.01)
dim(probeset.list1)

ids <- bitr(rownames(probeset.list1), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
df2 = probeset.list1[rownames(probeset.list1) %in% dedup_ids$ENSEMBL,]
df2$Y = dedup_ids$ENTREZID
kegg_gene_list <- df2$logFC
names(kegg_gene_list) <- df2$Y
kegg_gene_list<-na.omit(kegg_gene_list)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "mmu"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")


write.table(kk2@result,"/home/cblab/Documents/Barmingam/B3_Data/Enrichment/gene/I_0-24H_KEGG", sep="\t", quote=FALSE)

png(file="/home/cblab/Documents/Barmingam/plots/Enrichment/Gene/I_0-24H_KEGG_dotplot.png",width=1748, height=912)
dotplot(kk2, showCategory = 10, title = "Enriched Pathways", split=".sign") + facet_grid(.~.sign)
dev.off()

png(file="/home/cblab/Documents/Barmingam/plots/Enrichment/Gene/I_0-24H_KEGG_map.png",width=994, height=912)
emapplot(kk2)
dev.off()

png(file="/home/cblab/Documents/Barmingam/plots/Enrichment/Gene/I_0-24H_KEGG_ridge.png",width=1657, height=912)
ridgeplot(kk2) + labs(x = "enrichment distribution")
dev.off()
dme <- pathview(gene.data=kegg_gene_list, pathway.id="mmu04662", species = kegg_organism)

################## metabolomics enrichment (kegg pathway and modules)

pathway = read.table(file="/home/cblab/Documents/Barmingam/B3_Data/Enrichment/metabolomics/I-0_24_meta.tsv", sep="\t", header = T)
module = read.table(file="/home/cblab/Documents/Barmingam/B3_Data/Enrichment/metabolomics/I-0_24_module.tsv", sep="\t", header = T)

library(ggplot2)   # for plotting


Terms = as.character(pathway$Annotation[1:10])
Gene_ratio = pathway$n[1:10]/pathway$N[1:10]
p.adjust = round(-log(pathway$FDR.correction[1:10],10),0)
p_adjust = p.adjust

dat = as.data.frame(cbind(Terms,Gene_ratio,p.adjust))

png(file="/home/cblab/Documents/Barmingam/plots/Enrichment/Metabolome/I-0_24H_Kegg.png", width=912, height=912)
ggplot(data = dat, aes(x = p.adjust, y = Terms, 
                       color = p_adjust, size = Gene_ratio)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("-log10(p.adjust)") + 
  ggtitle("Pathway Enrichment")
dev.off()


modules = as.character(module$Annotation[1:10])
GeneRatio = round(module$n[1:10]/module$N[1:10],2)
p.adjust = round(-log(module$FDR.correction[1:10],10),2)
p_adjust = p.adjust

dat = as.data.frame(cbind(modules,GeneRatio,p.adjust))

png(file="/home/cblab/Documents/Barmingam/plots/Enrichment/Metabolome/I-0_24H_Module.png", width=912, height=912)
ggplot(data = dat, aes(x = p.adjust, y = modules, 
                       color = p_adjust, size = GeneRatio)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("-log10(p.adjust)") + 
  ggtitle("Module Enrichment")
dev.off()





















