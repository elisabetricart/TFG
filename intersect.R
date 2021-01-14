#Autora: Elisabte Ricart Ferrer
#Data: 14-01-21

#Carreguem el paquet de dades
library(VennDetail)

#Definim el directori de treball
setwd("C:/Users/elisa/Desktop/TFG/GSE75819")

#Llegim les taules de resultats
GSE75819<-read.table(file = "GSE75819/GSE75819_genes.txt",header = T,sep = "\t")
GSE65127<- read.table(file = "GSE65127/GSE65127_genes.txt", header=T, sep="\t")

#treiem la intersecció entre les dues llistes de gens
int<-intersect(GSE65127$sym,GSE75819$sym)

#fem el Venn Diagram
ven <- venndetail(list(GSE65127=GSE65127$sym,GSE75819=GSE75819$sym))
pdf("VennDetail.pdf")
plot(ven)
dev.off()
