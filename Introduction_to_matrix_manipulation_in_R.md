
# Introduction to matrix manipulation in R

- A variable can hold a single value

```{r}
x <- 1; # This variable has name 'x'. It is assigned to hold the number 1
print(x);
```
- A vector holds an ordered list of values 

```
x <- 1:10; #In this case, x is a vector of cosecurtive natural numbers 1, 2, 3,..., 10
print(x);
```

- We can perform simple arithmetic operations on a vector 

```
y <- x^2 - 2*x + 4;  # Simple algebraic operations can be applied to all the values in a vector. In this case, y is also a vector of 10 values.
print(y);
plot(x, y, type="b", main="Relationship between x and y");   # The 'plot' function generates a 2-dimensional x-y plot 

```

- A matrix holds an ordered arrangement of n x m numbers (n is number of rows, and m is number of columns) 

```
x <- matrix(1:15, nrow=5, ncol=3)  #x is assigned to be a matrix consisting of 2 rows and 5 columns, containing the values 1, 2, 3,..., 15
print(x);
```  

- **PCA plot**  

```
#install.packages('ggfortify')
library(ggfortify)
dat_pca=t(dat)
#head(dat_pca[,1:100])

##apply PCA - scale. = TRUE is highly advisable, but default is FALSE. 
out_pca <- prcomp(dat_pca,scale= TRUE)
plot(out_pca,type="l")
autoplot(out_pca,data=dat_pca,size=0.1,label=FALSE,label.size=5)
```  
  
  
- read in sample information  
Sample ID, "event", "time to event" and "clinical stage" gain  

see file *GSE102349_series_matrix_survival.csv*, clinical file should involved Sample ID, "event", "time to event" and "clinical_stage" extracted from file *GSE102349_series_matrix.txt.gz*. Here we also group samples based on clinical stage. 

```
info <- as.data.frame(read.table("GSE102349_series_matrix_survival.csv",header = TRUE,sep = ",", dec = ".",na.strings = "NA",stringsAsFactors=FALSE,check.names = FALSE))
row.names(info) <- info[,1]
info<-info[,-1]
dim(info)

#remove missing data and filter out clinical stage I and II samples
info=info[!info[, "time"] == "N/A",]
info=info[info[, "clinical_stage"] == "III" | info[, "clinical_stage"] == "IV",]

#set time to numeric datatype
info$time=as.numeric(as.character(info$time))

dim(info)
head(info)
table(info$`clinical_stage`)
```

- **survival plot**  

```
#install.packages("survminer")
#install.packages("survival")

library(survival)
library(ggplot2)
library(survminer)
library(dplyr)

#set time value to numeric 
info$time=as.numeric(as.character(info$time))

#set status value to numeric
a <- sub("Disease progression",1,info$event)
info$event <- as.numeric(as.character(sub("Last follow-up",0,a)))
head(info)

coxph(Surv(time, event)~clinical_stage, data=info)
fit <- survfit(Surv(time, event)~clinical_stage, data=info)
ggsurvplot(fit, conf.int=TRUE, pval=TRUE)
```


- **Perform differential gene expression analysis** using [limma](https://bioconductor.org/packages/release/bioc/html/limma.html)
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("limma")  

library("limma")
library("ggrepel")

###differential expression use limma###
exprSet<-dat[,row.names(info)]
head(exprSet)

par(cex = 0.7)
n.sample=ncol(exprSet)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
boxplot(exprSet, col = cols,main="expression value",las=2)


group<-info$clinical_stage
design <- model.matrix(~0+group)
colnames(design)=levels(factor(group))
rownames(design)=colnames(exprSet)
contrast.matrix<-makeContrasts(paste0(unique(group),collapse = "-"),levels = design)
fit <- lmFit(exprSet,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
tempOutput <- topTable(fit2, coef=1, n=Inf)
nrDEG <- na.omit(tempOutput) 

nrDEG_t<-as.data.frame(nrDEG)
colnames(nrDEG_t)<-c("logFC","AveExpr","t","P.Value","padj","B")
head(nrDEG_t)

```


- **VOLCONA PLOT**  
```
#draw picture
library(ggplot2)
DEG=nrDEG_t

colnames(DEG)
plot(DEG$logFC,-log10(DEG$P.Value))
#logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs( logFC)) )
logFC_cutoff=1  #set cut off
DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
table(DEG$change)
this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',])
)
g = ggplot(data=DEG, aes(,x=logFC, y=-log10(P.Value),color=change)) + geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+ xlab("log2 fold change") + ylab("-log10 p-value") + 
  ggtitle(this_tile) +  theme(plot.title = element_text(size=15,hjust = 0.5)) + 
  scale_colour_manual(values = c('blue','black','red')) ## corresponding to the levels(res$change)
print(g)
#ggsave(g,filename = 'volcano.pvalue.png')
```



- **HEATMAP PLOT** 

```
#install.packages("pheatmap")
library(stats)
library(ggplot2)
library(pheatmap)

#extract top 100 differential expressed genes and 
diff=nrDEG_t[order(nrDEG_t[,"logFC"],decreasing=TRUE),][1:100,] 
head(diff)
dim(diff)

#standalization
aa=t(dat[row.names(diff),])
aa=t(scale(aa))

p<-pheatmap(aa,show_rownames=F,show_colnames=F,cluster_cols=T, cluster_rows=T,cex=1, clustering_distance_rows="euclidean", cex=1,clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE)
p

#get samples/features ID
#colnames(aa[,p$tree_col[["order"]]]) 
#rownames(aa[p$tree_row[["order"]],])  

```




*************************
*Q1:How large difference among the results from pearson correlation, spearman collelation and kendall collelation*  
*Q2:which raw information should be extracted from GSE102349_series_matrix.txt file?*  
*Q3: What limma actually do? Why need to do normalization?  
*（additional）Try to understand each figures*
***************************




 
 
 
 
 
