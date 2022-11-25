library(destiny)
gene_read <- read.table(file = "~/P_3/read.txt",header = T, row.names = 1, sep = "\t")
gene_rpkm <- read.table(file = "~/P_3/rpkm.txt",header = T, row.names = 1, sep = "\t")
cluster_inf <- read.table(file = "~/P_3/tSNE_G3000_P100_Log10RPKM/ManuallyCluster.txt", header = T, row.names = 1, sep = "\t")
sample_inf_order <- read.table(file = "~/P_3/PseudoSpace/sample_inf_order.txt",header = T,row.names = 1,sep = "\t")
row.names(sample_inf_order) <- sample_inf_order$Library.name

cluster_num <- unique(cluster_inf$x)
cluster_num

i=1
read <- gene_read[,grep(cluster_num[i],cluster_inf$x)]
rpkm <- gene_rpkm[,grep(cluster_num[i],cluster_inf$x)]
meta.data <- sample_inf_order[grep(cluster_num[i],cluster_inf$x),]

plot(c(1:768), match((meta.data$LibraryName), colnames(read)))

source("~/script/p0.R")

gene_order = myMVG(reads = read, geneNumForPlot = 500, cv2Cutoff = 0.3, outIndex = "~/P_3/epi3/")

read = read[gene_order,]
rpkm = rpkm[gene_order,]

read = read[-grep("ERCC",rownames(read)),]
rpkm = rpkm[-grep("ERCC",rownames(rpkm)),]

#
read_MVG <- read[1:1000,]
sf.genes<-estimateSizeFactorsForMatrix(read_MVG) #estimate size factor
norm.data<-t( t(read_MVG) / sf.genes) 
dim(norm.data)
#####################
matrix <- log10(1 + norm.data)
dim(matrix)

test <- rowMeans(matrix)
genes_filter <- names(test)[test > log10(10+1)]
matrix <- matrix[na.omit(genes_filter),]
dim(matrix)
matrix = as.matrix(matrix)
sigmas <- find_sigmas(t(matrix), verbose = F)
dif.map.blue <- DiffusionMap(t(matrix), distance="cosine",sigma = optimal_sigma(sigmas))
plot(dif.map.blue, 1:2, pch=21,col="black", bg=as.character(meta.data$Embryonic.day))


df = dif.map.blue@eigenvectors[,1:3]
row.names(df) = colnames(matrix)
#df = read.table(file = "~/P_3/PseudoSpace/dfmap.txt")
write.table(x = df, file = "~/P_3/epi3d/df.new.txt", sep = "\t")
set.seed(10000)
sub2 = kmeans(x = df, centers = 3)
xxx = data.frame(x = sub2$cluster, row.names = colnames(matrix))
tmp = xxx$x
tmp = gsub(pattern = 1, replacement = "transition", x = tmp)
tmp = gsub(pattern = 2, replacement = "anterior", x = tmp)
tmp = gsub(pattern = 3, replacement = "posterior", x = tmp)
xxx$x = tmp


library(ggplot2)
g = ggplot(data = as.data.frame(df), mapping = aes(x = DC1, y = DC2, color = as.factor(sub2$cluster), shape = as.factor(meta.data$EmbryonicDay)))
g = g + geom_point()
#g = g + geom_text(aes(label=rownames(df_plot)),colour = "black", size = 0.5)
g = g + scale_color_manual(values = c("#999999","#E69F00","#56B4E9"))
g = g + scale_shape_manual(values = c(1,4,7,8))
g = g + theme_bw()
g = g + theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              legend.title=element_blank())
g = g + geom_text(aes(label=rownames(df)),colour = "black", size = 0.2)
#g = g + geom_line(df_cor_od_average, mapping = aes(x = aDC1, y = aDC2))
g


col = gsub(pattern = 1, replacement = "#999999", sub2$cluster)
col = gsub(pattern = 2, replacement = "#E69F00", col)
col = gsub(pattern = 3, replacement = "#56B4E9", col)

table(meta.data$EmbryonicDay)
shap = meta.data$EmbryonicDay
shap = gsub(pattern = "^5$", replacement = 1, x = shap)
shap = gsub(pattern = "^5.5$", replacement = 4, x = shap)
shap = gsub(pattern = "^6$", replacement = 7, x = shap)
shap = gsub(pattern = "^6.5$", replacement = 8, x = shap)
shap = as.integer(shap)


scatterplot3d::scatterplot3d(df[,1], df[,2], df[,3], color = col, phi = shap)

library(rgl)
library(plot3D)
within(
  data = as.data.frame(df),
  plot3d(DC1, DC2, DC3, col=col, size=5, pch = shap),
  points3d(DC1, DC2, DC3, col=col, size=6, pch = shap)
  #M <- par3d(DC1, DC2, DC3, col=col, size=5, pch = shap) 
)

movie3d( spin3d(), duration=10, dir="~/P_3/epi3d/movies/", clean=FALSE ) 
