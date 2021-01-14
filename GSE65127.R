#Autora: Elisabte Ricart Ferrer
#Data: 14-01-21

#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("affy")
#BiocManager::install("GEOquery")
#BiocManager::install("annotate")
#BiocManager::install("hgu133plus2.db")

#Carregar paquets
library(affy)
library(GEOquery)


#Definir el directori de treball
workingDir<-("C:/Users/elisa/Desktop/TFG/GSE75819")
setwd(workingDir)

###################
#Descarregar dades#
###################

gseid<-"GSE65127"#ID experiment
gse <- getGEO(gseid,GSEMatrix=TRUE) #descarreguem mostres normalizades
#phenoData: ens descarreguem les dades normalitzades per obtenir aquest pheno data (covariables approx)
pD <- pData(gse[[1]])

#Descarreguem dades crues (no normalitzades)
GEO<-getGEOSuppFiles(gseid) 
#descarreguem un fitxer comprimit (.tar)

setwd(paste(workingDir,gseid,sep="/")) #accedir a la carpeta GSE65127
untar("GSE65127_RAW.tar", exdir = getwd()) #descomprimit arxiu

#llegir fitxers amb les intensitats
celFiles = c(list.files( pattern = "_LST.CEL.gz"),list.files( pattern = "_NLST.CEL.gz"))
viti <- ReadAffy(filenames=as.vector(celFiles))# les 20 motres 
#escursar noms de les mostres
sampleNames(viti)<-c(substr(sampleNames(viti)[1:10],start = 17, stop = 23),substr(sampleNames(viti)[11:20],start = 17, stop = 24))

##################
#Data exploration#
##################
viti #20 mostres 54675 sondes
annotation(viti)
#paquet per annotar chip: "hgu133plus2"

##################################
#Quality control for 3'IVT arrays#
##################################

#box plots
colos<-rainbow(20) #agafem 20 colors per les mostres
pdf("boxplot_no_normalitzat_GSE65127.pdf")
boxplot(viti,col=colos,las=3,cex.axis=0.6) 
dev.off()

# density plots (histograma)
pdf("histograma_no_normalitzat_GSE65127.pdf")
hist(viti,col=colos)
legend("topright",sampleNames(viti),fill=colos,cex=0.6)
dev.off()

###############
#Normalization#
###############
#RMA
viti.rma<-rma(viti)

# box plots
pdf("boxplot_normalitzat_GSE65127.pdf")
boxplot(x = exprs(viti.rma),col=colos,las=3,cex.axis=0.6) 
dev.off()
# density plots (histograma)
pdf("histograma_normalitzat_GSE65127.pdf")
hist(viti.rma,col=colos)
dev.off()

#################################
#Sample aggregation: hclust & PCA
################################
x<-exprs(viti.rma) #assignar la matriu numèrica a objecte x

#Cluster jerarquic
pdf("clust_jerarquic_euclidian.pdf")
clust.cor.ward<- hclust(dist(t(x)),method="ward.D2")
plot(clust.cor.ward, main="hierarchical clustering(euclidian)", hang=-1,cex=0.8)
dev.off()

#Principal Component Analysis
library(scatterplot3d)
pdf("PCA_GSE65127.pdf")
#Calculem Principal Components
summary(pca.filt <- prcomp(t(x), scale=T )) 
#Fem el gràfic 3D de les PCs
pca3d<-scatterplot3d(x=pca.filt$x[,1],y=pca.filt$x[,2],z=pca.filt$x[,3],  xlab='PC1 (15.61%)', ylab='PC2 (11.68%)', zlab='PC3 (8.97%)', main='PCA', pch=16,color=colos,col.grid="lightblue")
text(pca3d$xyz.convert(pca.filt$x+0.5), labels=rownames(pca.filt$x), cex=0.6)
dev.off()

##############
#####DEG######
##############
library(limma)
#creem el vector amb la condicio a comparar
cond<-as.factor(c(rep("Lesional",10), rep("NonLesional",10)))

#creem una matriu que assigni CTRL o VITIL a cada mostra(1/0)
design<-model.matrix(~0+cond)
#posem nom de les mostres a les files
rownames(design)<-sampleNames(viti.rma) 
#eliminem el tros "cond" del nom de les columnes
colnames(design)<-gsub("cond","",colnames(design)) 

#apliquem el model lineal
fit<-lmFit(viti.rma,design) 
#contrast matrix ens permet comprar entre els nivell de les catagories
contrast.matrix<-makeContrasts(Lesional - NonLesional,levels=design) ###Lesional- Nolesional
#apliquem la nostra constrat matrix al nostre model lineal 
fit2<-contrasts.fit(fit,contrast.matrix)
fite<-eBayes(fit2)

#corregim els valors per multiple testing (FDR)
top.table<-topTable(fite,coef=1,number=Inf,adjust="fdr")
head(top.table)

#seleccionem els signifcatius (p-valor < 0.05)
results.adj.p0.05<-top.table[top.table$adj.P.Val<0.05,]
dim(results.adj.p0.05)
#135 sondes diferencialment expressades

##############
#Annotation
##############
#Carreguem els paquets per anotar
library(annotate)
library(hgu133plus2.db) 

#obtenim el simbol, el nom del gen i el chromosoma de cada sonda
sym<-Filter(Negate(is.na), (mget(rownames(results.adj.p0.05), env=hgu133plus2SYMBOL)))
name<-mget(names(sym), env=hgu133plus2GENENAME)
chr<-mget(names(sym), env=hgu133plus2CHR)

#obtenim els valors d'intensitats, el logFC, el p-valor i el p-valor ajustat de cada mostra
dat <- exprs(viti.rma)[names(sym),] 
logFC <- results.adj.p0.05[names(sym),]$logFC
pval <- results.adj.p0.05[names(sym),]$P.Value
adj.pval<-results.adj.p0.05[names(sym),]$adj.P.Val

#Creem un html amb els resultats 
affyids<-names(sym)
genelist <- list(names(sym))
filename <- "Results_GSE65127.html"
title <- "Differentially expressed VITIL vs CTRL"
othernames <- list(sym,name,chr,round(logFC, 1), round(pval, 4), round(adj.pval, 4), round(dat, 2)) 
head <- c("Probe ID", "Symbol", "Gene Name", "Chr","logFC", "p-value","adj.p-value",sampleNames(viti.rma))
repository <- list("affy")
htmlpage(genelist, filename, title, othernames, head, repository = repository)

#Creem una taula txt amb els resultats
sign_table <- as.data.frame((cbind(unlist(sym),unlist(name),unlist(chr),unlist(logFC),unlist(pval),unlist(adj.pval))))
colnames(sign_table) <- c("sym","name","chr","logFC", "pval", "adj.p-val")
length(unique(sign_table$sym))
write.table(file = "GSE65127_genes.txt", x = sign_table, sep = "\t", quote = F, row.names = T)

##############
#Results plots
##############
#volcano plot
pdf("volcano_GSE65127.pdf")
volcanoplot(fite,coef=1,highlight=10,names=fite$genes$NAME,main="Volcano plot", xlim = c(-6, 6))
dev.off()

#heatmap
pdf("heatmap65127.pdf")
heatmap(dat)
dev.off()
