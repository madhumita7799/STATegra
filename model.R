library(openintro)
library(tidyverse)
library(rsample)
library(doParallel)
library(caret)
library(iml)
library(patchwork)
library(mlbench)
library(Boruta)
library(TH.data)
library(xgboost)
library(shapr)

mirna = read.table(file="miRNA/miRNA_final.txt", sep="\t", header=T)
mrna = read.table(file="mRNA/mRNA_final.csv", sep=",", header=T, row.names=1)
met = read.table(file="Metabolomics/Metabolomics_final.txt", sep="\t", header=T, row.names=1)

class_label = c(rep(0,18),rep(1,18))
class = c(rep("C",18),rep("I",18))

mirna1 = mirna
mrna1 = mrna
met1 = met

colnames(met1)=colnames(mrna1)=colnames(mirna1)

integrated_data = rbind(mirna1, mrna1, met1)
data_wt_cl = cbind(t(integrated_data), class_label)

###############################################################################
#### feature selection
gene = cbind(t(mrna),class)
micro = cbind(t(mirna),class)
metabo = cbind(t(met),class)

train = as.data.frame(gene,row.names = F)
train_micro = as.data.frame(micro,row.names = F)
train_metabo = as.data.frame(metabo,row.names = F)

boruta_output1 <- Boruta(class ~ ., data=na.omit(train), doTrace=0) 
signif_gene <- getSelectedAttributes(boruta_output1, withTentative = FALSE)

boruta_output_micro <- Boruta(class ~ ., data=na.omit(train_micro), doTrace=0) 
signif_micro <- getSelectedAttributes(boruta_output_micro, withTentative = FALSE)
library("stringr")
signif_micro= str_replace_all(signif_micro, "`", "")

boruta_output_metabo <- Boruta(class ~ ., data=na.omit(train_metabo), doTrace=0) 
signif_met <- getSelectedAttributes(boruta_output_metabo, withTentative = FALSE)


gene_subset =  gene[,signif_gene]
mirna_subset =  micro[,signif_micro]
met_subset = metabo[,signif_met]

integrated_subset = cbind(gene_subset, mirna_subset, met_subset)

write.table (integrated_subset, "Integrated_subset", quote = F, row.names = T, col.names = T, sep = "\t")
write.table (gene_subset, "gene_subset", quote = F, row.names = T, col.names = T, sep = "\t")
write.table (mirna_subset, "mirna_subset", quote = F, row.names = T, col.names = T, sep = "\t")
write.table (met_subset, "met_subset", quote = F, row.names = T, col.names = T, sep = "\t")

############## classification ###########3
#### for gene 
#########################
gene_subset = cbind(gene_subset, class_label)
class(gene_subset) <- "numeric"
gene_subset = as.data.frame(gene_subset)

set.seed(1993)
intrain <- createDataPartition(y = gene_subset$class_label, p= 0.7, list = FALSE)
training <- gene_subset[intrain,]
testing <- gene_subset[-intrain,]
training[["class_label"]] = factor(training[["class_label"]])
testing[["class_label"]] = factor(testing[["class_label"]])
summary(testing)

trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
svm_Linear <- train(class_label ~., data = training, method = "svmLinear",
                    trControl=trctrl,
                    preProcess = c("center", "scale"),
                    tuneLength = 10)

test_pred <- predict(svm_Linear, newdata = testing)
confusionMatrix(table(test_pred, testing$class_label))

############# for mirna
micro_subset = cbind(mirna_subset, class_label)
class(micro_subset) <- "numeric"
micro_subset = as.data.frame(micro_subset)

set.seed(1993)
intrain <- createDataPartition(y = micro_subset$class_label, p= 0.7, list = FALSE)
training <- micro_subset[intrain,]
testing <- micro_subset[-intrain,]
training[["class_label"]] = factor(training[["class_label"]])
testing[["class_label"]] = factor(testing[["class_label"]])
summary(testing)

trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
svm_Linear <- train(class_label ~., data = training, method = "svmLinear",
                    trControl=trctrl,
                    preProcess = c("center", "scale"),
                    tuneLength = 10)

test_pred <- predict(svm_Linear, newdata = testing)
confusionMatrix(table(test_pred, testing$class_label))

############# for metabolomics
met_subset = cbind(met_subset, class_label)
class(met_subset) <- "numeric"
met_subset = as.data.frame(met_subset)

set.seed(1993)
intrain <- createDataPartition(y = met_subset$class_label, p= 0.7, list = FALSE)
training <- met_subset[intrain,]
testing <- met_subset[-intrain,]
training[["class_label"]] = factor(training[["class_label"]])
testing[["class_label"]] = factor(testing[["class_label"]])
summary(testing)

trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
svm_Linear <- train(class_label ~., data = training, method = "svmLinear",
                    trControl=trctrl,
                    preProcess = c("center", "scale"),
                    tuneLength = 10)

test_pred <- predict(svm_Linear, newdata = testing)
confusionMatrix(table(test_pred, testing$class_label))

############### for fused dataset
integrated_subset= cbind(integrated_subset, class_label)
class(integrated_subset) <- "numeric"
integrated_subset = as.data.frame(integrated_subset)

set.seed(1993)
intrain <- createDataPartition(y = integrated_subset$class_label, p= 0.7, list = FALSE)
training <- integrated_subset[intrain,]
testing <- integrated_subset[-intrain,]
training[["class_label"]] = factor(training[["class_label"]])
testing[["class_label"]] = factor(testing[["class_label"]])
summary(testing)

trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
svm_Linear <- train(class_label ~., data = training, method = "svmLinear",
                    trControl=trctrl,
                    preProcess = c("center", "scale"),
                    tuneLength = 10)
test_pred <- predict(svm_Linear, newdata = testing)
confusionMatrix(table(test_pred, testing$class_label))


######################## 
library(knitr)
library(printr)
library(pvRsuper)












