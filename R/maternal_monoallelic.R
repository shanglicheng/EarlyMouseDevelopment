#sample.female.sheet <- read.table(file = "~/P_3/XCI/CASTCells.txt", header = T, row.names = NULL, sep = "\t")
#sample_inf_order <- read.table(file = "~/P_3/XCI/sample_inf_order_CAST.txt", header = T, row.names = 1, sep = "\t")
#lineage_clust <- read.table(file = "~/P_3/tSNE_G3000_Log10RPKM/ManuallyCluster.txt")
#sample_inf_order <- data.frame(sample_inf_order,Lineage = lineage_clust$x)
#sample.female.sheet$Embryo_ID
#sample_inf_order_cast <- sample_inf_order[c(grep(11,sample_inf_order$Embryo),
#                                            grep(12,sample_inf_order$Embryo),
#                                            grep(23,sample_inf_order$Embryo),
#                                            grep(24,sample_inf_order$Embryo),
#                                            grep(27,sample_inf_order$Embryo),
#                                            grep(28,sample_inf_order$Embryo)),]
#write.table(sample_inf_order_cast, file = "~/P_3/XCI/sample_inf_order_cast_selected.txt", sep = "\t")

sample_inf_order_cast <- read.table(file = "~/P_3/XCI/sample_inf_order_cast_selected.txt", header = T, row.names = 1, sep = "\t")

#############
#############
source("~/P_3/XCI/read_rpkmfiles.r")

t1=read_rpkmfiles("~/P_3/XCI/refseq_mm9_allelehits.txt")
snp.alle.hits_1640=t1$counts
row.names(snp.alle.hits_1640)=t1$genes
colnames(snp.alle.hits_1640)=t1$samples
dim(snp.alle.hits_1640)

write.table(x = colnames(snp.alle.hits_1640), file = "~/P_3/XCI/snp.alle.hits_1640.txt", sep = "\t")

#note that this SNP analysis pipeline is based on mm9 database, and you can get mm9 gene annotation 
#from  /mnt/kauffman/sandberglab/pipeline3.0/additionaldata/gene-bed-files/mm9refgene.bed
mm9.bed=read.csv("~/P_3/XCI/mm9refgene.bed", sep="\t", header = F)
mm9.bed=data.frame(symbol=as.vector(sub(":.*", "", mm9.bed[,4])) , mm9.bed)
colnames(mm9.bed)[2:4]=c("chr", "startp", "endp")
mm9.bed[1:2, 1:5]
nrow(unique(mm9.bed[,1:2]))

# genes on chr Y and X
mm9.chrY=mm9.bed[mm9.bed[,"chr"]=="chrY",]
# order the genes on chr
mm9.chrY.symbol=as.vector(unique(mm9.chrY$symbol))

mm9.chrX=mm9.bed[mm9.bed[,"chr"]=="chrX",]
mm9.chrX.symbol.uniq <- as.vector(unique(mm9.chrX$symbol))
mm9.chrX.unique <- mm9.chrX[match(mm9.chrX.symbol.uniq, mm9.chrX$symbol),]
mm9.chrX.unique.order= mm9.chrX.unique[order(mm9.chrX.unique$startp),]


# to select all chrX genes, and select your cells groups for analysis

sample.female.select_1640 <- as.character(sample_inf_order_cast$LibraryName)
dim(snp.alle.hits_1640)

snp.alle.hits_1640.df <- as.data.frame(snp.alle.hits_1640)

str(mm9.chrX.unique.order)
mm9.chrX.unique.order$symbol <- as.character(mm9.chrX.unique.order$symbol)
str(mm9.chrX.unique.order)

allehits.chrX_1640 <- snp.alle.hits_1640.df[mm9.chrX.unique.order$symbol,]
allehits.chrX_1640 <- allehits.chrX_1640[-grep("NA.*",rownames(allehits.chrX_1640)),]

#
#female and male chromosome
Chr_female <- paste(sample_inf_order_cast$LibraryName,"_",sample_inf_order_cast$FemaleChr,"only", sep = "")
Chr_male <- paste(sample_inf_order_cast$LibraryName,"_",sample_inf_order_cast$MaleChr,"only", sep = "")

#expression of female and male chromosome 
Chr_female_alle_hit <- allehits.chrX_1640[,Chr_female]
Chr_male_alle_hit <- allehits.chrX_1640[,Chr_male]

#order by Epi(1), Exe(2), VE(3)
Chr_female_alle_hit <- cbind(cbind(Chr_female_alle_hit[,grep("1",sample_inf_order_cast$Lineage)],
                                   Chr_female_alle_hit[,grep("2",sample_inf_order_cast$Lineage)]),
                             Chr_female_alle_hit[,grep("3",sample_inf_order_cast$Lineage)])
Chr_male_alle_hit <- cbind(cbind(Chr_male_alle_hit[,grep("1",sample_inf_order_cast$Lineage)],
                                   Chr_male_alle_hit[,grep("2",sample_inf_order_cast$Lineage)]),
                             Chr_male_alle_hit[,grep("3",sample_inf_order_cast$Lineage)])
sample_inf_order_cast_order <- rbind(rbind(sample_inf_order_cast[grep("1",sample_inf_order_cast$Lineage),],
                                           sample_inf_order_cast[grep("2",sample_inf_order_cast$Lineage),]),
                                     sample_inf_order_cast[grep("3",sample_inf_order_cast$Lineage),])
#pvalue between F and M
pvalue <- NULL
for (i in (1:ncol(Chr_female_alle_hit))) {
  tmp <- wilcox.test(Chr_female_alle_hit[,i],Chr_male_alle_hit[,i],paired = TRUE, conf.level = 0.95)
  pvalue[i] <- tmp$p.value
}
#cell expression
Chr_female_alle_hit_sum <- apply(Chr_female_alle_hit,2,sum)
Chr_male_alle_hit_sum <- apply(Chr_male_alle_hit,2,sum)

#maternal monoallelic expression
MME <- Chr_female_alle_hit_sum/(Chr_female_alle_hit_sum+Chr_male_alle_hit_sum)
#
for (i in (1:length(MME))) {
  if (pvalue[i] >= 1) {
    MME[i] = 0.5
  }
}
table(sample_inf_order_cast_order$Lineage)
nrow(sample_inf_order_cast_order)
MME_df <- data.frame(y = c(c(1:197),seq(198,243,2),c(244:416)), 
                     x = MME*100, 
                     cast_inf = sample_inf_order_cast_order$FemaleChr, 
                     pvalue = pvalue,
                     lineageinf = factor(sample_inf_order_cast_order$Lineage, levels = c("1","2","3")),
                     EmbryonicDay = factor(sample_inf_order_cast_order$EmbryonicDay),
                     LineageCAST = paste(sample_inf_order_cast_order$Lineage,sample_inf_order_cast_order$FemaleChr,sep = ""))

str(MME_df)
library(ggplot2)
g <- ggplot(data = MME_df, mapping = aes(x = x, y = y))
g <- g + geom_vline(aes(xintercept=0), colour="black", linetype="dotted", size=1)
g <- g + geom_vline(aes(xintercept=20), colour="gray", linetype="dotted")
g <- g + geom_vline(aes(xintercept=50), colour="black", linetype="dotted", size=1)
g <- g + geom_vline(aes(xintercept=80), colour="gray", linetype="dotted")
g <- g + geom_vline(aes(xintercept=100), colour="black", linetype="dotted", size=1)

g <- g + geom_hline(aes(yintercept=197.5), colour="black", linetype="dotted")
g <- g + geom_hline(aes(yintercept=243.5), colour="black", linetype="dotted")
g <- g + geom_hline(aes(yintercept=416.5), colour="black", linetype="dotted")

g <- g + geom_point(data = MME_df,mapping = aes(colour = EmbryonicDay, shape = LineageCAST))
g <- g + theme_bw()
g <- g + theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               axis.line = element_line(colour = "black"))
g <- g + theme(panel.border = element_blank())
g <- g + scale_x_continuous(name = "Allelic Expression",limits = c(0,100), breaks = c(0,20,50,80,100))
g <- g + scale_y_continuous(limits = c(-0,417))
g
