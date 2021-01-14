#Autora: Elisabte Ricart Ferrer
#Data: 14-01-21


#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("annotate")
#BiocManager::install("illuminaHumanv3.db")

#Carregeum els paquets
library(RSQLite)
library(GEOquery)
library(lumi)
library(umap)
library(illuminaio)

#Definim el directori de treball
workingdir<-("C:/Users/elisa/Desktop/TFG/GSE75819")
setwd(workingdir)

###################
#Descarregar dades#
###################

gseid<-"GSE75819" #id experiment
gse <- getGEO(gseid,GSEMatrix=TRUE) #descarreguem mostres normalizades
#phenoData: ens descarreguem les dades normalitzades per obtenir aquest pheno data (covariables approx)
pD <- pData(gse[[1]])

# Descarreguem dades crues (no normalitzades)
GEO<-getGEOSuppFiles(gseid) 
#descarreguem un fitxer comprimit (.tar)
setwd(paste(getwd(),gseid,sep="/")) #accedir a la carpeta GSE65127

#Llegim les dades
GSE75819<-lumiR("GSE75819_non-normalized.txt.gz")


##################
#Data exploration#
##################
GSE75819 #30 mostres i 48803 sondes


##################################
#######Quality control############
##################################

#boxplot
colos<-rainbow(30) #agafem colors per les mostres
pdf("boxplot_no_normalitzat_GSE75819.pdf")
boxplot(GSE75819,col=colos,las=3,cex.axis=0.6) 
dev.off()

# density plots (histograma)
pdf("histograma_no_normalitzat_GSE75819.pdf")
hist(GSE75819,col=colos)
dev.off()

###############
#Normalization#
###############
GSE75819.norm<-lumiN(GSE75819,"quantile")

##boxplot
pdf("boxplot_norm_GSE75819.pdf")
boxplot(GSE75819.norm,col=colos,las=3,cex.axis=0.6) 
dev.off()
#density plot (histograma)
pdf("histograma_normalitzat_GSE75819.pdf")
hist(GSE75819.norm,col=colos)
dev.off()

#transformaci? en base log2
exprs(GSE75819.norm) <- log2(exprs(GSE75819.norm))

#################################
#Sample aggregation: hclust & PCA
################################
x<-exprs(GSE75819.norm) #assignar la mtriu d'expressi? normal. a x

#Cluster jerarquic
pdf("cluster_jerarquic_GSE75819.pdf")
clust.cor.ward<- hclust(dist(t(x)),method="ward.D2")
plot(clust.cor.ward, main="hierarchical clustering(euclidian)", hang=-1,cex=0.8)
dev.off()

##Principal Component Analysis
library(scatterplot3d)
#Calculem les principal components
summary(pca.filt <- prcomp(t(x), scale=T )) 
#Fem el gràfic en 3D de les PC
pdf("PCA_norm_GSE75819.pdf")
pca3d<-scatterplot3d(x=pca.filt$x[,1],y=pca.filt$x[,2],z=pca.filt$x[,3],  xlab='PC1 (2.31%)', ylab='PC2 (2.23%)', zlab='PC3 (2.14%)', main='PCA', pch=16,color=colos,col.grid="lightblue")
text(pca3d$xyz.convert(pca.filt$x+0.5), labels=rownames(pca.filt$x), cex=0.6)
dev.off()

##############
#####DEG######
##############
library(limma)
#creem el vector amb la condicio a comparar
cond<-as.factor(as.vector(pD$source_name_ch1))

#creem una matriu que assigni CTRL o VITIL a cada mostra(1/0)
design<-model.matrix(~0+cond)
#posem nom de les mostres a les files
rownames(design)<-sampleNames(GSE75819.norm) 
#eliminem el tros "cond" del nom de les columnes
colnames(design)<-c("Lesional","Normal") 

#apliquem el model lineal
fit<-lmFit(GSE75819.norm,design) 
#contrast matrix ens permet comprar entre els nivell de les catagories
contrast.matrix<-makeContrasts(Lesional-Normal,levels=design) ###Lesional- Nolesional
#apliquem la nostra constrat matrix al nostre model lineal 
fit2<-contrasts.fit(fit,contrast.matrix)
fite<-eBayes(fit2)

#corregim els valors per multiple testing (FDR)
top.table<-topTable(fite,coef=1,number=Inf,adjust="fdr")
head(top.table)

#seleccionem els signifcatius (p-valor < 0.05)
results.adj.p0.05<-top.table[top.table$adj.P.Val<0.05,]
dim(results.adj.p0.05)
#9417 sondes diferencialment expressades 

##############
#Annotation
##############
#Carreguem els paquets per anotar
library(annotate)
library(illuminaHumanv3.db)

#obtenim el simbol, el nom del gen i el chromosoma de cada sonda
sym<-Filter(Negate(is.na), (mget(rownames(results.adj.p0.05), env=illuminaHumanv3SYMBOL)))
name<-mget(names(sym), env=illuminaHumanv3GENENAME)
chr<-mget(names(sym), env=illuminaHumanv3CHR)

#obtenim els valors d'intensitats, el logFC, el p-valor i el p-valor ajustat de cada mostra
dat <- exprs(GSE75819.norm)[names(sym),] 
logFC <- results.adj.p0.05[names(sym),]$logFC
pval <- results.adj.p0.05[names(sym),]$P.Value
adj.pval<-results.adj.p0.05[names(sym),]$adj.P.Val

#Creem un html amb els resultats 
affyids<-names(sym)
genelist <- list(names(sym))
filename <- "Results_GSE75819.html"
title <- "Differentially expressed VITIL vs CTRL"
othernames <- list(sym,name,chr,round(logFC, 1), round(pval, 4), round(adj.pval, 4), round(dat, 2)) 
head <- c("Probe ID", "Symbol", "Gene Name", "Chr","logFC", "p-value","adj.p-value",sampleNames(GSE75819.norm))
repository <- list("affy")
htmlpage(genelist, filename, title, othernames, head,repository = repository)

#Creem una taula txt amb els resultats
sign_table <- (cbind(unlist(sym),unlist(name),unlist(chr),unlist(logFC),unlist(pval),unlist(adj.pval)))
sign_table <- as.data.frame(sign_table)
colnames(sign_table) <- c("sym","name","chr","logFC", "pval", "adj.p-val")
write.table(file = "GSE75819_genes.txt", x = sign_table, sep = "\t", quote = F, row.names = T)


##############
#Results plots
##############
#volcano plot
pdf("volcano_GSE75819.pdf")
volcanoplot(fite,coef=1,highlight=10,names=fite$genes$NAME,main="Volcano plot", xlim = c(-6, 6))
dev.off()

#heatmap
pdf("heatmap75819.pdf")
heatmap(dat)
dev.off()
