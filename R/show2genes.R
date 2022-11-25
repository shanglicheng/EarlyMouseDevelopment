{
col.g1<-c(103,39,112)/255              #colour associated with gene 1
col.g2<-c(col2rgb("darkorange"))/255   #colour associated with gene 1
n<-c(1,1,1) #no gene expressed

f<-c(0,0.961,0.961) #both genes expressed


A=col.g1[1]-n[1]
B=col.g2[1]-n[1]
C=n[1]

D=col.g1[2]-n[2]
E=col.g2[2]-n[2]
Ff=n[2]

G=col.g1[3]-n[3]
H=col.g2[3]-n[3]
I=n[3]

alpha<-f[1]-A-B-C
beta<-f[2]-D-E-Ff
gamma<-f[3]-G-H-I
}

#############
source("~/script/p0.R")
xList = myP3()
read = xList$reads
rpkm = xList$rpkms
sample_inf_order = xList$sample_inf_order
tsne_corr = xList$tSNE_coor

####
#difmap_corr = read.table(file = "~/P_3/PseudoSpace/Epi/DiffusionMapCoor.txt", header = T, row.names = 1, sep = "\t")
#tsne_corr = difmap_corr
####

subCluster = "1"
read_sub = read[,grep(subCluster, sample_inf_order$Lineage)]
rpkm_sub = rpkm[,grep(subCluster, sample_inf_order$Lineage)]
tsne_corr_sub = tsne_corr[grep(subCluster, sample_inf_order$Lineage),]
sample_inf_order = sample_inf_order[grep(subCluster, sample_inf_order$Lineage),]

match(sample_inf_order$LibraryName, colnames(read_sub))

gene1 = "Wnt3"
gene2 = "Klf2"#gene2Vector[i]#"Wnt3"

data.1.norm = as.vector(t(log2(rpkm_sub[gene1,]+1)))
data.2.norm = as.vector(t(log2((rpkm_sub[gene2,]+1))))

data.1.norm = (data.1.norm - min(data.1.norm))/(max(data.1.norm) - min(data.1.norm))
data.2.norm = (data.2.norm - min(data.2.norm))/(max(data.2.norm) - min(data.2.norm))

data.norm<-cbind(data.1.norm, data.2.norm)#data.1.norm and data.2.norm are expression levels of gene 1 and gene 2 (normalised to vary between 0 and 1)

#
#data.norm[,2] = rep(0, nrow(data.norm))

z1<-apply(data.norm, 1, function(x) max(c(0,A*x[1]+B*x[2]+C+alpha*x[1]*x[2])))
z2<-apply(data.norm, 1, function(x) max(c(0,D*x[1]+E*x[2]+Ff+beta*x[1]*x[2])))
z3<-apply(data.norm, 1, function(x) max(c(0,G*x[1]+H*x[2]+I+gamma*x[1]*x[2])))

#z1,z2,z3 are the RGB coordinates associated to each cell.

shapeTag = gsub("^5$", replacement = "1", sample_inf_order$EmbryonicDay)
shapeTag = gsub("^5.5$", replacement = "4", shapeTag)
shapeTag = gsub("^6$", replacement = "7", shapeTag)
shapeTag = gsub("^6.5$", replacement = "8", shapeTag)
shapeTag = as.numeric(shapeTag)

with(tsne_corr_sub, plot(x = V1, y = V2, pch = shapeTag, col = rgb(z1,z2,z3), lwd = 0.5))
