# 差异分析 --------------------------------------------------------------------

#load library
library(tidyr)
library(dplyr)
library(stringi)
library(stringr)
library(limma)
library(edgeR)
library(DESeq2)  

#####################################################################################2.3.疾病与正常组间差异表达基因、


gene_exp1 <- read.table('D:\\项目\\YQGY-10202-9\\原始数据\\TCGA\\HiSeqV2',
                        sep = '\t',check.names = F,header = T,quote = "",row.names = 1)
# 
# exprSet <-data
exprSet <-gene_exp1
# D:\项目\YQNN-10201-9\原始数据\TCGA
# exprSet <- exp1

# 
# 
# #肿瘤样本保留01和02去进行TIDE分析
# exprSet <- exprSet[,c(1,which(substr(colnames(exprSet),14,15)=='11'))]
# rownames(exprSet) <- exprSet$sample
# gene_exp1 <- exprSet
# write.table(gene_exp1,file="TCGA.symbol.txt",sep="\t",row.names = F,quote = F)
# 
# 
# # #正常样本保留11，肿瘤样本保留01和02
exprSet <- exprSet[,c(which(substr(colnames(exprSet),14,15)=='11'),
                      which(substr(colnames(exprSet),14,15)=='01'))]

nNormal <- ncol(exprSet[,which(substr(colnames(exprSet),14,15)>10)])
nTumor <- ncol(exprSet[,which(substr(colnames(exprSet),14,15)<10)])

gene_exp1 <-exprSet

# # exprSet <- t(exprSet)
# 
# # 
# gene_exp1 <- aggregate(.~Tag,gene_exp1,max)
# # 
# # 
# rownames(gene_exp1) <- gene_exp1[,1]
# # gene_exp1 <- gene_exp1[-1,-1]
# # # 
# 
# # 
# 
# # load("mRNAmatrix.symbol.rda")
# orde_tes  <- read.table('FISSION.sample.txt',
#                         sep = '\t',check.names = F,header = T)
# 
# name <- orde_tes$Accession
# gene_exp1 <- gene_exp1[,name]
#
# gene  <- read.table('LSCC-sample.txt',
#                     sep = '\t',check.names = F,header = T)
# #
# gene_exp1 <- gene_exp1[,gene$sampleID]
# write.table(gene_exp1,file="GSE16011.mRNAmatrix.txt",sep="\t",row.names = T,quote = F)
# 
# # Tag <- rownames(gene_exp1)
# 
# gene_exp1 <-cbind(Tag,gene_exp1)
# 
save(gene_exp1,file='LSCC.mRNAmatrix.symbol.rda')


# write.table(gene_exp1,file="TCGA.gene_exp1.gene.txt",sep="\t",row.names = T,quote = F)


# load('GSE12720.mRNAmatrix.symbol.rda')

# rownames(gene_exp1) <- gene_exp1[,1]

# 
load("LSCC.mRNAmatrix.symbol.rda")

write.table(gene_exp1, 'LSCC.mRNAmatrix.symbol.txt', sep = '\t', col.names = NA, quote = FALSE)


# 求log
# dat1 <-apply(gene_exp1,2,function(x){log2(x+1)})

# 逆转log
dat1 <-round(apply(gene_exp1,2,function(x){2^x-1}))
gene_exp1 <- dat1

# DESeq2 


# 加载包
library(DESeq2)   
library(pheatmap)  # 用于作热图的包
library(ggplot2)   # 用于作图的包
# 读入数据，注意设置工作路径
# countData <- as.matrix(read.csv("gene_count_matrix.csv",row.names="gene_id"))
# countData <-  read.table('F:\\项目\\YCZK-129\\原始数据\\mRNA.txt',
#                          sep = '\t',check.names = F,header = T,row.names = 1,quote = "")


countData <-gene_exp1

# countData <- countData[rowMeans(countData)>1,]              # 去除表达量过低的基因




condition <- factor(c(rep("control",nNormal),rep("Keys",nTumor)))

colData <- data.frame(row.names=colnames(countData), condition)


dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)

dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE) 
#将结果用result()函数来获取
res <- results(dds1)

summary(res) 


# res格式转化：用data.frame转化为表格形式
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
# 依次按照pvalue值log2FoldChange值进行排序
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

write.table(res1, './01.limma/DESeq2.mRNA.txt', sep = '\t', col.names = NA, quote = FALSE)


# 获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
res1_up<- res1[which(res1$log2FoldChange >= 2 & res1$padj < 0.05),]      # 表达量显著上升的基因
res1_down<- res1[which(res1$log2FoldChange <= -2 & res1$padj < 0.05),]    # 表达量显著下降的基因
res1_total <- rbind(res1_up,res1_down)


write.table(res1_total, './01.limma/DESeq2.mRNA.out.txt', sep = '\t', col.names = NA, quote = FALSE)



allDiff <- res1
diffSig <- res1_total


# mirna_lncrna_dabase$geneName
# #
# # #
# 
# gene_exp1 <- combat_edata
# gene  <- read.table('F:\\项目\\YQXX0268\\原始数据\\Glutamine.txt',
#                     sep = '\t',check.names = F,header = T)
# 
# gene_exp1 <- gene_exp1[gene$Gene,]
# # gene_exp1 <- gene_exp1[mirna_lncrna_dabase$geneName,]
# # df[is.na(gene_exp1)] <- 0
# # 
# write.table(gene_exp1,file="TCGA.gene.exp.txt",sep="\t",row.names = T,quote = F)
# 
# gene  <- read.table('Methylation.ssgsea.test.txt',
#                     sep = '\t',check.names = F,header = T)
# 
# gene_exp1 <- gene_exp1[,gene$sample]

# # 
# # 
# nNormal <- ncol(gene_exp1[,-1][,which(substr(colnames(gene_exp1[,-1]),14,15)>10)])
# nTumor <- ncol(gene_exp1[,-1][,which(substr(colnames(gene_exp1[,-1]),14,15)<10)])
# 
# gene_exp2 <- as.data.frame(gene_exp1[,-1])
# load("GSE16011.mRNAmatrix.symbol.rda")

# gene_exp1 <- gene_exp1[,c(18:47)]
# gene_exp1 <- gene_exp1[,c(1:17,32:47)]
# 



gene_exp2 <- as.data.frame(gene_exp1)

# nNormal <- 88
# nTumor <- 308
# # 
# nNormal <- 8
# nTumor <- 31
# # 
# nNormal <- 233
# nTumor <- 138
# 
nNormal <- 8
nTumor <- 276
# colnames(gene_exp2)
nNormal <- 142
nTumor <- 142

ID <- colnames(gene_exp2)
dim(gene_exp2)
Type <- c(rep('Control',nNormal),rep('AS',nTumor))

# Type <- c(rep('Ovarian_cancer',nTumor),rep('Control',nNormal))
ID <- cbind(ID,Type)
colnames(ID) <- c("id","Type")
ID <- as.data.frame(ID)

gene_exp3 <- t(gene_exp2)
id <- rownames(gene_exp3)
gene_exp3 <- cbind(id,gene_exp3)

# head(exp1)[1:5,1:5]

riskTime <- merge(ID,gene_exp3,by = 'id')
rownames(riskTime) <- riskTime$id
riskTime <- riskTime[ID$id,]

# riskTime <- riskTime[order(riskTime$Type),F]

group_list= riskTime$Type %>% factor(.,levels = c("Control","AS"),ordered = F)
pvalue <-0.05
logFoldChange <- 0.5 #调节差异值



design = model.matrix(~0+group_list, data=group_list)
colnames(design) <- levels(group_list)
rownames(design) <-rownames(riskTime)


gene_exp1[is.na(gene_exp1)] <- 0
dat1 <-gene_exp1

# dat1 <- scale(dat1,center = T,scale = T)
# dat1 <- scale(dat1)
# dat1 <- apply(gene_exp1[,-1],2,function(x){log2(x+1)})
# str(dat1)

# df[is.na(dat1)] <- 0
# dat1 <-apply(gene_exp1,2,function(x){log2(x+1)})
max(dat1)
min(dat1)
# dat1 <-dat1+1



dat1 <- dat1[,rownames(design)]

contrast.matrix <- makeContrasts(AS-Control, levels = design)#Case比Control
# 这里是死亡的比幸存的

fit <- lmFit(dat1, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
options(digits = 4)#输出结果的小数点后保留4位


allDiff=topTable(fit2,coef=1,number=Inf)
allDiff <- na.omit(allDiff)
write.table(cbind(Symbol=rownames(allDiff),allDiff),file="01.limma/lncRNA.limma.txt",sep="\t",row.names = F,quote = F)
diffSig = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC>logFoldChange | allDiff$logFC<(-logFoldChange))),]
# allDiff$adj.P.Val
# diffSig <- diffSig[-14,]
write.table(cbind(Symbol=rownames(diffSig),diffSig),file="01.limma/lncRNA.limma.OUT.txt",sep="\t",row.names = F,quote = F)

# write.table(diffexp,file="diffexp.txt",sep="\t",row.names = T,quote = F)


# head(gene_exp1)[1:5,1:5]
# 
# 
# exp(4)
# power(4)
# log()


# 密度复热图 --------------------------------------------------------------------


##加载包

library(circlize)
library(dendextend)
library('ComplexHeatmap')
library('circlize')
library("RColorBrewer")
library(pheatmap)

# gene_exp1 <-apply(gene_exp1,2,function(x){log2(x+1)})
diff <- rownames(diffSig)

# cir1 <- gene_exp1[diff,-1]

cir1 <- gene_exp1[diff,]

# cir1 <- mat1[1:50,]




# scale对其归一化处理
cir1 = scale(cir1, center = TRUE, scale = TRUE)
# 自定义颜色
mycols <- colorRamp2(breaks = c(-2, 0, 2),
                     colors = c("#269ccf", "white", "#f0962f"))

# 配置分组详细信息
ha1 = HeatmapAnnotation(group = c(rep("Control", 12), rep("LSCC", 116)))
ha2 = HeatmapAnnotation(expression = anno_points(rnorm(128)))
densityHeatmap(cir1, top_annotation = ha1, bottom_annotation = ha2,col =c("#269ccf","white","#f6d04d","#f0962f"))

# 绘图

pdf(file = '01.limma/多重热图.pdf',width = 8,height = 8)
densityHeatmap(cir1, top_annotation = ha1, bottom_annotation = ha2,col =c("#269ccf","white","#f6d04d","#f0962f")) %v%
  # HeatmapAnnotation(Expression = anno_barplot(1:63)) %v%
  # Heatmap(cir1, name = "Type",height = unit(12, "cm"),
  #         show_row_names = F,show_column_names = F, col = mycols,cluster_columns=T,cluster_rows=T, km  = 4) %v%
  
  ComplexHeatmap::pheatmap(cir1,
                           # cellwidth = 6,cellheight = 0.5,
                           # method="spearman", #计算gene或sample之间的相关性的方法，可选"pearson" (default), "kendall", or "spearman"
                           scale="row", #为基因做scale
                           cluster_rows=T,#为基因做聚类##不要给样本做聚类，否则样本名顺序会变
                           cluster_cols=F,#为sample做聚类
                           color = mycols,###表征表达量高低的颜色域
                           show_colnames=F,show_rownames =F,
                           # annotation_col = annotation_col,
                           treeheight_row = "1",treeheight_col = "1",#不画树
                           border_color = "NA",
                           heatmap_height = unit(11, "cm") #图形高度 
  )
dev.off()





# 渐变火山图 -------------------------------------------------------------------


# devtools::install_github("BioSenior/ggVolcano")

library(ggVolcano)

# 
# #读取上述输出的差异倍数计算结果
# gene_diff <- read.delim('edgeR.txt', row.names = 1, sep = '\t', check.names = FALSE)
# #首先对表格排个序，按 FDR 值升序排序，相同 FDR 值下继续按 log2FC 降序排序
# gene_diff <- gene_diff[order(gene_diff$FDR, gene_diff$logFC, decreasing = c(FALSE, TRUE)), ]
# 
# pvalue <-0.05
# logFoldChange <- 1 #调节差异值
# #log2FC≥1 & FDR<0.01 标识 up，代表显著上调的基因
# #log2FC≤-1 & FDR<0.01 标识 down，代表显著下调的基因
# #其余标识 none，代表非差异的基因
# gene_diff[which(gene_diff$logFC >= 1 & gene_diff$FDR < 0.05),'sig'] <- 'up'
# gene_diff[which(gene_diff$logFC <= -1 & gene_diff$FDR < 0.05),'sig'] <- 'down'
# gene_diff[which(abs(gene_diff$logFC) <= 1 | gene_diff$FDR >= 0.05),'sig'] <- 'none'
# 
# #输出选择的差异基因总表
# gene_diff_select <- subset(gene_diff, sig %in% c('up', 'down'))
# write.table(gene_diff_select, file = 'edgeR.out.txt', sep = '\t', col.names = NA, quote = FALSE)
# 
# 
# 
# allDiff <- gene_diff
# diffSig <- gene_diff_select

# allDiff[1,5] <- 3.904e-18

# 如果你的数据没有名为 "regulate "的列，你可以使用add_regulate函数来添加，使用方法如下。
# use the function -- add_regulate to add a regulate column 
# to the DEG result data. 
allDiff <- add_regulate(allDiff, log2FC_name = "log2FoldChange",
                        fdr_name = "padj",log2FC = 2, fdr = 0.05)
allDiff$name <- rownames(allDiff)


# plot 如果用了add_regulate x、y都不用修改
pdf(file="01.limma/差异火山图.pdf",width=6,height=6)
gradual_volcano(allDiff, x = "log2FoldChange", y = "padj",legend_position = "UL",
                label = "name", label_number = 10, output = FALSE,
                log2FC_cut = 2,FDR_cut = 0.05)
dev.off()

?gradual_volcano()





# 渐变韦恩图 -------------------------------------------------------------------


# install.packages("ggVennDiagram")

library(ggVennDiagram)
library(ggplot2)
library(patchwork)


# 读取韦恩图文件
gene_diff <- read.delim('ggVennDiagram1.txt', sep = '\t', check.names = FALSE)


# genes <- paste0("gene",1:1000)
# set.seed(20210302)


gene_list <- list(LASSO = gene_diff$LASSO,
                  # yellow = gene_diff$yellow,
                  # DEGs = gene_diff$DEGs,
                  Boruta = gene_diff$Boruta)
# List of 4

#最简单的用法
ggVennDiagram(gene_list)

# 分类名 category name，可通过category.names参数修改
ggVennDiagram(gene_list, category.names = c("AA","BB","CC"))


# 选取渐变色
RColorBrewer::display.brewer.all()


pdf(file = 'Veen1.pdf',width = 10,height = 8)
# label参数：c("both", "count", "percent", "none") 四选一，默认为第一个；
ggVennDiagram(gene_list, edge_lty = "solid", edge_size = 2, label = "count",scale_color="black"
              , label_size  = 8,set_size = 10)+
  # scale_x_continuous(expand = expansion(mult = .2))+
  # scale_fill_gradient(low="#269ccf",high = "#f0962f")
  scale_fill_distiller(palette = "RdBu")+ 
  scale_color_manual(values = c("#fff7f7","#fff7f7","#fff7f7","#fff7f7"))+
  scale_x_continuous(expand = expansion(mult = .3))

# scale_fill_distiller(palette = "Reds")
dev.off()

?ggVennDiagram()



ggVennDiagram(gene_list, edge_lty = "solid", edge_size = 2)+
  # scale_x_continuous(expand = expansion(mult = .2))+
  # scale_fill_gradient(low="#269ccf",high = "#f0962f")
  scale_fill_distiller(palette = "RdBu")+ 
  scale_color_manual(values = c("#fff7f7","#fff7f7","#fff7f7","#fff7f7"))


#GO、KEGG富集 ---------------------------------------------------------------



# gene <- diff

library(org.Hs.eg.db) #人类注释数据库
library(clusterProfiler)#进行GO富集和KEGG富集
library(dplyr) #进行数据转换
library(ggplot2)#绘图

gene <- c("DAO","HSPA13","LNX1","LNX1","RGCC")
gene <- c("SEC61G","UNG")

gene <- rownames(diffSig)
# 
gene <- read.table('overgene.txt',
                   sep = '\t',check.names = F,header = T)

gene <- gene[,1]
gene

hg<-bitr(gene,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb="org.Hs.eg.db")
id <- hg[,2]
id

# id <- c("100526794",id)


###GO富集
go1 <- enrichGO(gene,OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 1, qvalueCutoff = 1,keyType = 'SYMBOL')
dim(go1)#查看富集结果

write.csv(go1,file="go.csv")#导出GO富集的结果



#进行结果可视化
barplot(go,showCategory=20,drop=T)#柱状图

dotplot(go,showCategory=20)#气泡图

#通过ggplot2将BP、MF、CC途径的富集结果挑选前8条绘制在一张图上
pdf(file = "GO.gene.pdf", width = 8, height = 6)
barplot(go1, split="ONTOLOGY",font.size=13)+ facet_grid(ONTOLOGY~.,scale="free")
dev.off()
dotplot(go, split="ONTOLOGY",font.size=10)+ facet_grid(ONTOLOGY~.,scale="free")

##KEGG富集
ego <- enrichKEGG(
  gene = id,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 1,
  pAdjustMethod  = "BH",
  qvalueCutoff = 1,
  use_internal_data =T
)

write.csv(ego,file="KEGG.csv")#导出KEGG富集的结果
dim(ego)#查看富集结果

pdf(file = "KEGG.gene.pdf", width = 8, height = 6)
barplot(ego,font.size=15,showCategory=6,color = "pvalue")	# 画气泡图
barplot(ego,font.size=15,showCategory=8)	# 画气泡图
dev.off()
dotplot

browseKEGG(ego,'mmu01100')	# 显示通路图











##富集铉图

library(GOplot)
library(openxlsx)

d1 = read.xlsx("EC-david-GO.xlsx")
d2 = read.xlsx("EC-genelist.xlsx")

circ = circle_dat(d1,d2)

d3 = c("heart development","vasculature development","blood vessel development","tissue morphogenesis","blood vessel morphogenesis")#这个是指定你要做那几个GO的弦图。这个必须于david文件中的名字一样。

# chord = chord_dat(circ,d2,d3)
chord = chord_dat(circ,d2)
head(chord)

pdf("KEGG.chord.pdf",height = 18,width = 18)#准备好一块画布。
GOChord(chord,space = 0.02,gene.order = 'logFC',gene.space = 0.25,gene.size = 8)#进行绘图。
dev.off()#然后必须关闭并保存这块画布。






# 八卦图

#安装GOplot包
# BiocManager::install("GOplot")
#加载GOplot包
library(GOplot)

#加载数据
# data(EC)

# Go <- EC$david
# gene<- EC$genelist
# 
# EC$david <- Go
# EC$genelist <- gene
# 
# 


gene <- read.table('limmun.txt',
                   sep = '\t',check.names = F,header = T)
Go <- read.table('kegg.txt',
                 sep = '\t',check.names = F,header = T)

?GOCircle()
?circle_dat()
#创建一个pdf文件,用来保存circleplot
pdf(file="KEGG_circle_plot.pdf",width=14,height=8)
#创建circ对象，第一个参数为GO富集分析结果，第二个参数为差异表达分析结果
circ <- circle_dat(Go, gene)
# circ <- circ[18:29,]
#绘制circleplot图
GOCircle(circ,nsub=6,label.size = 6)
dev.off()


GOCircle(circ)


# 最好看的富集图用微生信
http://www.bioinformatics.com.cn/plot_basic_go_pathway_circlize_plot_140
输入文件 富集炫图.txt


# Treemap图
# install.packages("treemap")
library(treemap)
library(RColorBrewer)
ExampleWOS<-data.frame(group=c("USA","USA","China","France","USA","Ireland","USA","France","UK","USA"),
                       subgroup=c("UNIVERSITY OF CALIFORNIA SYSTEM","HARVARD UNIVERSITY","CHINESE ACADEMY OF SCIENCES",
                                  "INSTITUT NATIONAL DE LA SANTE ET DE LA RECHERCHE MEDICALE INSERM",
                                  "UNIVERSITY OF NORTH CAROLINA","UNIVERSITY COLLEGE CORK",
                                  "UNIVERSITY OF CALIFORNIA SAN DIEGO","INRAE","UNIVERSITY OF LONDON",
                                  "UNIVERSITY OF TEXAS SYSTEM"),
                       
                       value=c(4.770,3.083,1.821,1.758,1.578,1.495,1.495,1.466,1.389,1.374)
)
treemap(ExampleWOS,index=c("group","subgroup"),vSize="value",type="index",palette = "Set1")



data <- read.csv("KEGG.WGCNA.txt", sep = "\t", header = T,row.names = 1)
head(data)


pdf(file="KEGG_WGCNA.pdf",width=14,height=8)

treemap(data,
        index="Description", #指定分组的列
        vSize="Count", #指定面积大小的列
        vColor="p.adjust", #指定颜色深浅的列
        fontsize.labels=c(10, 10), #设置标签字体大小
        align.labels=list(c("center", "center"), c("left", "top")), #设置标签对齐的方式
        border.col="black", #设置边框的颜色  
        border.lwds=c(2,2),#设置边框的线条的宽度
)



treemap(data,
        index="Description", #指定分组的列
        vSize="Count", #指定面积大小的列
        vColor="p.adjust", #指定颜色深浅的列
        type="value", #指定颜色填充数据的类型
        palette='-Blues',
        fontsize.labels=c(10, 10), #设置标签字体大小
        align.labels=list(c("center", "center"), c("left", "top")), #设置标签对齐的方式
        border.col="black", #设置边框的颜色  
        border.lwds=c(2,2)#设置边框的线条的宽度
)

dev.off()
# WGCNA -------------------------------------------------------------------



#####################################################################################2.2.WGCNA寻找与免疫细胞相关模块

###########################################2、WGCNA


# .libPaths("/data2/project-YangHM/library/")


###########################################2、WGCNA

#.libPaths()
#.libPaths(c("/data2/Rlibrary/","/data2/project-YangHM/R/",.libPaths()))

setwd("/data/project-YangHM/project/YQGY-10202-9")

#load library

library(WGCNA)

library(stringr)

#根据上一步得到的常规表达矩阵和分组信息
# gene_exp <- read.table('GSE75011.maxtire.txt',
#                        sep = '\t',check.names = F,header = T)
# #
# 
# gene_exp1 <- gene_exp
# #
# colnames(gene_exp)
# # # 
# # gene_exp1 <- gene_exp[which(rowSums(gene_exp[,-1]) > 0),]
# rownames(gene_exp1) <- gene_exp1$Gene
# gene_exp1 <- aggregate(.~gene_exp,gene_exp1,max)

# rownames(gene_exp1) <- gene_exp1$Gene.Symbol
# save(gene_exp1,file='GSE75011.exprSet.rda')
load('LSCC.mRNAmatrix.symbol.rda')

# dat1 <- apply(gene_exp1,2,function(x){log2(x+1)})

# sample <- read.table('WGCNA.mRNAsi.txt',
#                      sep = '\t',check.names = F,header = T)
# 

dat1 <- gene_exp1


dataExpr = t(dat1)



gsg = goodSamplesGenes(dataExpr, verbose = 3);
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
dim(dataExpr)


#再过滤一遍
# dataExpr <- dataExpr[,colSums(dataExpr) <10*nrow(dataExpr)*0.9]#count :removing all features that have a count of less than say 10 in more than 90% of the samples
#或者大于10%的样本都表达
# dataExpr <- dataExpr[,colSums(dataExpr>0)>nrow(dataExpr)*0.1]
dim(dataExpr)
# save(dataExpr,file = 'lncRNA.WGCNA.exprSet.rda')
# load(file = 'lncRNA.WGCNA.exprSet.rda')

dim(dataExpr)
##dataExpr <- dataExpr[-c(4,11,138),]###去除离群值
sampleTree = hclust(dist(dataExpr), method = "complete")
pdf(file = "Fig1.Sample_cluster.pdf", width = 6.5, height = 5)
par(mar = c(0,4,2,0),xpd=F)
save(sampleTree,file = 'Sample_cluster.Rdata')
# load("Sample_cluster.Rdata")
plot(sampleTree, main = "Sample clustering to detect outliers", cex=0.3,sub="", xlab="",cex.main=1.5)
# abline(h = 91, col = "red")
dev.off()









#制作样本性状
# sample <- read.table('WGCNA.sample.txt',
#                      sep = '\t',check.names = F,header = T)

# sample$Type1 <- sample$score
# ##sample$Type1 <- c(rep('0',76),rep('1',78))
# names(sample)

# sample$Astrocytes
# sample$`B-cells`
# sample$CLP
# sample$Eosinophils
# sample$Erythrocytes
# sample$`Memory_B-cells`
# sample$MEP
# sample$`naive_B-cells`
# sample$Neurons
# sample$`pro_B-cells`

# 
# dataTraits <- data.frame(Regeneration=sample$Astrocytes,sample$`B-cells`,sample$CLP,sample$Eosinophils,sample$Erythrocytes,sample$`Memory_B-cells`,sample$MEP,sample$`naive_B-cells`,sample$Neurons,sample$`pro_B-cells`)
# rownames(dataTraits) <- sample$sample
# colnames(dataTraits) <- c("Astrocytes","B_cells","CLP","Eosinophils","Erythrocytes","Memory_B_cells","MEP","naive_B_cells","Neurons","pro_B_cells")
# head(dataTraits)
# 
# 
# sampleTree2 = hclust(dist(dataExpr), method = "complete")
# traitColors = numbers2colors(as.numeric(dataTraits$"Astrocytes",dataTraits$"B_cells"), signed = FALSE)
# pdf("Fig2.cluster.dendrogram.pdf",width=6.5,height=5)
# plotDendroAndColors(sampleTree2, traitColors,
#                     groupLabels = colnames(dataTraits),cex.dendroLabels=0.3,
#                     main = "Sample dendrogram and trait heatmap")
# dev.off()
# save(sampleTree2,file = 'cluster.dendrogram.Rdata')





# 
# #制作样本性状
sample <- read.table('WGCNA.type_1.txt',
                     sep = '\t',check.names = F,header = T)

dataTraits <- data.frame(Regeneration=sample$"Type")
rownames(dataTraits) <- sample$Sample
colnames(dataTraits) <- c("Type")
head(dataTraits)

# sample$GSVA
# 
# dataTraits <- data.frame(Regeneration=sample$"GSVA",sample$"UC")
# rownames(dataTraits) <- sample$sample
# colnames(dataTraits) <- c("Control","UC")
# head(dataTraits)


sampleTree2 = hclust(dist(dataExpr), method = "complete")
# traitColors = numbers2colors(as.numeric(dataTraits$"UC",dataTraits$"Control"), signed = FALSE)
traitColors = numbers2colors(as.numeric(dataTraits$"Type"), signed = FALSE)
pdf("Fig2.cluster.dendrogram.pdf",width=6.5,height=5)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = colnames(dataTraits),cex.dendroLabels=0.3,
                    main = "Sample dendrogram and trait heatmap")
dev.off()
# save(sampleTree2,file = 'cluster.dendrogram.Rdata')

# as.numeric(dataTraits$"aDC",dataTraits$"Macrophages",dataTraits$"Mast cells",dataTraits$"Neutrophils",dataTraits$"NK CD56bright cells",dataTraits$"NK cells",dataTraits$"T cells")




###power值散点图
enableWGCNAThreads()   #多线程工作
powers =seq(from = 1, to=20, by=1)  #幂指数范围1:20
sft = pickSoftThreshold(dataExpr, powerVector = powers, verbose = 5)
save(sft,file='sft.Rdata')
load('sft.Rdata')
pdf('Fig2.Soft.Threshold.pdf',width = 10,height = 6)

par(mfrow = c(1,2))
cex1 = 0.85
###拟合指数与power值散点图
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.9,col="red") #可以修改
###平均连通性与power值散点图
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()




###邻接矩阵转换
softPower =sft$powerEstimate #最佳power值
softPower

softPower =9
# #服务器跑
#softPower
adjacency = adjacency(dataExpr, power = softPower)

# 
# ###TOM矩阵,服务器跑
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
save(dissTOM,file='dissTOM.Rdata')
# load('dissTOM.Rdata')
###基因聚类，服务器跑
geneTree = hclust(as.dist(dissTOM), method = "average")
save(geneTree,file='geneTree.Rdata')





load('geneTree.Rdata')
pdf(file="Gene.cluster.pdf",width=6,height=5)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()





###动态剪切模块识识别
#服务器跑

minModuleSize =  100  #模块基因数目
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
save(dynamicMods,file='dynamicMods.Rdata')
load('dynamicMods.Rdata')
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
pdf(file="Fig4.DynamicTree.pdf",width=6.5,height=4.5)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()





###相似模块聚类
MEList = moduleEigengenes(dataExpr, colors = dynamicColors)
MEs = MEList$eigengenes
#如果不合并，运行下两行
# moduleColors = MEList$validColors
# table(moduleColors)
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
pdf('Fig5.merge.cluster.pdf',width = 6.5,height = 5)
par(mar=c(2,3,3,3))
#png('Fig4.merge.cluster.png',width = 600,height = 400)
plot(METree, main = "Clustering of module eigengenes",cex=1,
     xlab = "", sub = "")
MEDissThres = 0.35  #剪切高度可修改
abline(h=MEDissThres, col = "red")
dev.off()
###相似模块合并
merge = mergeCloseModules(dataExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
pdf(file="Fig6.DynamicTree.pdf",width=6.5,height=5)
plotDendroAndColors(geneTree, cbind(dynamicColors,mergedColors),c("Dynamic Tree Cut","Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
moduleColors = mergedColors
save(moduleColors,file = 'moduleColors.Rdata')
table(moduleColors)
colorOrder = c("MEmidnightblue", standardColors(30))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs





###模块与性状数据热图
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
moduleTraitCor = cor(MEs, dataTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
write.table(moduleTraitPvalue,'moduleTraitPvalue.txt',sep = '\t',quote = F,row.names = T)
write.table(moduleTraitCor,'moduleTraitCor.txt',sep = '\t',quote = F,row.names = T)
pdf(file="Fig7.Module_trait.pdf",width=8,height=8)
par(mar = c(5, 10, 3, 3))


#必须至少是两列，如果是一列就如下处理
dataTraits1 <- cbind(dataTraits,dataTraits)
# labeledHeatmap(Matrix = cbind(moduleTraitCor,moduleTraitCor),
#                xLabels = colnames(dataTraits1),
#                yLabels = names(MEs),
#                ySymbols = names(MEs),
#                colorLabels = FALSE,
#                colors = blueWhiteRed(500),
#                textMatrix = cbind(textMatrix,textMatrix),
#                setStdMargins = FALSE,
#                cex.text =1,
#                main = paste("Module-trait relationships"))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(dataTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(500),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text =1,
               main = paste("Module-trait relationships"))
dev.off()











###计算MM和GS值
#modNames = substring(names(MEs), 3)
#modNames <- c('brown','greenyellow','yellow','cyan','grey60','pink')
modNames <- unique(moduleColors)
geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
#names(geneModuleMembership) = paste("MM", modNames, sep="")
#names(MMPvalue) = paste("p.MM", modNames, sep="")
traitNames=colnames(dataTraits)
geneTraitSignificance = as.data.frame(cor(dataExpr, dataTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")




###输出每个模块的基因
for (mod in 1:nrow(table(moduleColors)))
{  
  modules = names(table(moduleColors))[mod]
  probes = colnames(dataExpr)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0('WGCNA.',modules,".txt"),sep=",",row.names=F,col.names=F,quote=F)
}

###输出GS_MM数据
probes = colnames(dataExpr)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "GS_MM.csv",sep=",",row.names=F)








#基因与性状的相关性
pdf('Fig7.Gene.trait.MEyellow.pdf',width = 15,height = 17)
# par(mfrow=c(1,2),mar=c(6,5,5,3))
par(mar=c(6,5,5,3))
for (i in colnames(geneTraitSignificance)){
  moduleGenes = moduleColors=='yellow'
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, 'MEyellow']),
                     abs(geneTraitSignificance[moduleGenes,i]),
                     xlab = ("Module Membership in midnightblue module"),
                     ylab = paste("Gene significance"),
                     main = paste0('Trait ',substr(i,4,nchar(i)),"\n"),
                     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, col = 'yellow')
  abline(v=0.6,h=0.2,col="red")
}
dev.off()

###输出每个模块的基因
for (mod in 1:nrow(table(moduleColors)))
{  
  modules = names(table(moduleColors))[mod]
  probes = colnames(dataExpr)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0('WGCNA.',modules,".txt"),sep=",",row.names=F,col.names=F,quote=F)
}


#多重机器学习 ------------------------------------------------------------------


dir.create ("04.LASSO")
#############随机森林

library(randomForest)
library(ROCR)
library(genefilter)
library(Hmisc)

#####获取差异基因表达量

load('LSCC.mRNAmatrix.symbol.rda')

dim(gene_exp1)

# gene_exp1 <-apply(gene_exp1,2,function(x){log2(x+1)})

gene <- read.table('overgene.txt',
                   sep = '\t',check.names = F,header = T)


gene_exp <- gene_exp1[gene$Symbol,]

gene_exp_EDG <- t(gene_exp)


trait <- read.table('LSCC.sample.txt',
                    sep = '\t',check.names = F,header = T)

colnames(gene_exp)

trait2 <- cbind(trait,gene_exp_EDG)
Key <- as.factor(trait2$Type)
trait2 <- cbind(Key,trait2)
str(trait2)

trait2 <- as.data.frame(trait2[,-c(2,3)])

# colnames(trait2)[9] <- "ABHD14A"
# colnames(trait2)[11] <- "BSCL2"


set.seed(100)
rf=randomForest(Key ~ ., trait2, importance = TRUE, ntree = 400, proximity=TRUE)
pdf("rf.pdf",width = 5,height = 5)
plot(rf)
dev.off()
#查看变量的重要性
rf

importance <- importance(x=rf)

#绘制变量的重要性图
pdf("varImpPlot.pdf",width = 8,height = 5)
varImpPlot(rf)
dev.off()
write.table(importance,'importance.txt',sep="\t",quote=F,col.names=T)




# Boruta特征选择鉴定关键分类变量
# install.packages("Boruta")
library(Boruta)
set.seed(2)

boruta <- Boruta(Key ~ ., trait2, pValue=0.01, mcAdj=T, 
                 maxRuns=300)

boruta
boruta$finalDecision

boruta$impSource


# plot(boruta)
pdf("./04.LASSO/boruta.pdf",width = 5,height = 5)
plot(boruta, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(boruta$ImpHistory),function(i)
  boruta$ImpHistory[is.finite(boruta$ImpHistory[,i]),i])
names(lz) <- colnames(boruta$ImpHistory)  
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(boruta$ImpHistory), cex.axis = 0.7)
dev.off()
#蓝色的盒状图对应一个阴影属性的最小、平均和最大Z分数。
#红色、黄色和绿色的盒状图分别代表拒绝、暂定和确认属性的Z分数。

boruta.finalVars <- data.frame(Item=getSelectedAttributes(boruta, withTentative = F), Type="Boruta")
Tentative.boruta <- TentativeRoughFix(boruta)


write.table(Labels,'./04.LASSO/FinalDecision.txt',sep="\t",quote=F,col.names=F)



#############LASSO


library(tidyverse)
library(glmnet)
source('msvmRFE.R')   #文件夹内自带
library(VennDiagram)
library(sigFeature)
library(e1071)
library(caret)
library(randomForest)
#####获取差异基因表达量

load('Input5.GSE96804.unique.all.rda')

dim(gene_exp1)



gene <- read.table('module.gene.txt',
                   sep = '\t',check.names = F,header = F)


gene_exp <- gene_exp1[gene$V1,]

gene_exp_EDG <- t(gene_exp[,-1])


trait <- read.table('GSE96804.sample.txt',
                    sep = '\t',check.names = F,header = T)

colnames(gene_exp)

trait2 <- cbind(trait,gene_exp_EDG)
Key <- as.factor(trait2$Type)
trait2 <- cbind(Key,trait2)
str(trait2)

dim(trait2)
trait2 <- as.data.frame(trait2[,-(2:4)])
dim(trait2)

# ###### 分别存储 自变量和因变量
# y <- trait2[,1]
# x <- trait2[,2:14]
# 
# set.seed(8)
# 
# 
# model.info <- data.frame()
# for (seed in 1:100){
set.seed(3)

y <- trait2[,1]
x <- trait2[,2:9]
alpha1_fit <- glmnet(x,y,alpha=1,family="binomial")


pdf("./04.LASSO/lambda.pdf",width = 5,height = 5)
plot(alpha1_fit,xvar="lambda",label=TRUE)
dev.off()


alpha1.fit <- cv.glmnet(data.matrix(x),y,type.measure = "class",alpha=1,family="binomial")
pdf("./04.LASSO/alpha1.fit.pdf",width = 5,height = 5)
plot(alpha1.fit)
abline(v=log(c(alpha1.fit$lambda.min,alpha1.fit$lambda.1se)),lty="dashed")#8.8是文字在横轴的位置，可以AI调整
text(log(alpha1.fit$lambda.min),0.03,cex=1.5,
     labels = paste0('lambda.min = \n',round(alpha1.fit$lambda.min,2)))#8.8是文字在横轴的位置，可以AI调整
text(log(alpha1.fit$lambda.1se),0.08,cex=1.5,
     labels = paste0('lambda.lse = \n',round(alpha1.fit$lambda.1se,2)))
dev.off()
print(alpha1.fit)

coef <- coef(alpha1_fit,s=alpha1.fit$lambda.1se)
index <- which(as.matrix(coef)!= 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
lassoGene


myCoefs <- coef(alpha1.fit, s="lambda.min")
myCoefs


alpha1.fit$lambda.min


myCoefs <- coef(alpha1.fit, s="lambda.min");
lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
lasso_fea <- lasso_fea[-1]
lasso_fea

myCoefs@x

write.table(myCoefs@x,'./04.LASSO/lasso.myCoefs.txt',sep="\t",quote=F,col.names=F)


# model.info <- rbind(model.info,
#                     cbind(data.frame(Seed = seed,MG=paste(lasso_fea,sep = " ",collapse=","))))
# }





#####跑种子

model.info <- data.frame()
for (seed in 1:500){
  set.seed(seed)
  
  y <- trait2[,1]
  x <- trait2[,2:11]
  alpha1_fit <- glmnet(x,y,alpha=1,family="binomial")
  
  
  # pdf("lambda.pdf",width = 5,height = 5)
  # plot(alpha1_fit,xvar="lambda",label=TRUE)
  # dev.off()
  
  
  alpha1.fit <- cv.glmnet(data.matrix(x),y,type.measure = "class",alpha=1,family="binomial")
  # pdf("alpha1.fit.pdf",width = 5,height = 5)
  # plot(alpha1.fit)
  # dev.off()
  print(alpha1.fit)
  
  coef <- coef(alpha1_fit,s=alpha1.fit$lambda.1se)
  index <- which(as.matrix(coef)!= 0)
  actCoef <- coef[index]
  lassoGene=row.names(coef)[index]
  lassoGene
  
  
  myCoefs <- coef(alpha1.fit, s="lambda.min")
  myCoefs
  
  
  alpha1.fit$lambda.min
  
  
  myCoefs <- coef(alpha1.fit, s="lambda.min");
  lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
  lasso_fea <- lasso_fea[-1]
  lasso_fea
  
  
  model.info <- rbind(model.info,
                      cbind(data.frame(Seed = seed,MG=paste(lasso_fea,sep = " ",collapse=","))))
}
write.table(model.info, file="./04.LASSO/model.info.lass28.txt",sep="\t",quote=F,row.names=F)




lasso_fea

##############SVM


library(forecast)
library(e1071)
source(msvmRFE.R)
# input <- trait2[,c("JUN","CDKN1A","HMGB2","DDIT3","DDIT4")]
input <- trait2

source("F:\\workspace\\msvmRFE.R")

#采用五折交叉验证 (k-fold crossValidation）
svmRFE(input, k = 5, halve.above = 100) #分割数据，分配随机数

nfold = 8
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=5, halve.above=100) #特征选择


top.features = WriteFeatures(results, input, save=F) #查看主要变量
head(top.features)

#把SVM-REF找到的特征保存到文件
write.csv(top.features,"./04.LASSO/feature_svm.csv")


# 运行时间主要取决于选择变量的个数，一般的电脑还是不要选择太多变量
# 选前5个变量进行SVM模型构建，体验一下
featsweep = lapply(1:8, FeatSweep.wrap, results, input) #13个变量
save(featsweep,file = "featsweep.1.RData")
featsweep

# 选前300个变量进行SVM模型构建，然后导入已经运行好的结果
#featsweep = lapply(1:300, FeatSweep.wrap, results, input) #300个变量
#save(featsweep,file = "featsweep.RData")
load("featsweep.1.RData")

no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

#dev.new(width=4, height=4, bg='white')
pdf("B_svm-error.pdf",width = 5,height = 5)
PlotErrors(errors, no.info=no.info) #查看错误率

dev.off()





#dev.new(width=4, height=4, bg='white')
pdf("B_svm-accuracy.pdf",width = 5,height = 5)
Plotaccuracy(1-errors,no.info=no.info) #查看准确率
dev.off()

# 图中红色圆圈所在的位置，即错误率最低点
which.min(errors) 


top.features[1:which.min(errors),]

write.csv(top.features[1:which.min(errors),],"feature_svm.tezheng.csv")







#####xgboost算法
# install.packages("xgboost")
library(xgboost)


library(xgboost)
library(readr)
library(stringr)
library(caret)
library(car)


library("xgboost")

library("Matrix")

###xgboost仅适用于数值型向量，因此在训练模型前需要对数据进行相应的转化操作。
trait2 <- as.data.frame(trait2[,-(1:2)])
dim(trait2)
colnames(trait2)

# colnames(trait2)[9] <- "ABHD14A"
# colnames(trait2)[11] <- "BSCL2"


train_matrix <- sparse.model.matrix(Type ~ .-1, data = trait2)
train_label <- as.numeric(trait2$Type)

train_fin <- list(data=train_matrix,label=train_label)
dtrain <- xgb.DMatrix(data = train_fin$data, label = train_fin$label)


set.seed(9)
xgb <- xgboost(data = dtrain,max_depth=6, eta=0.5, 
               
               objective='binary:logistic', nround=25)
importance <- xgb.importance(train_matrix@Dimnames[[2]], model = xgb) 

head(importance)
pdf("xgboost.pdf",width = 5,height = 5)
xgb.ggplot.importance(importance)
dev.off()
importance$Gain
importance$Feature


write.table(importance$Feature,'xgboost1.txt',sep="\t",quote=F,col.names=F)


# 单基因GSEA分析 ---------------------------------------------------------------
dir.create ("06.GSEA")
# 
# install.packages("rlang")
# library(semPlot)
# install.packages("semPlot", dependencies = TRUE)

# 以该基因表达值的中位值来对样本进行分组
library ("DESeq2")
library(dplyr)
library(GSEABase)
library(clusterProfiler)
library(enrichplot)
library(RColorBrewer)

library(tidyr)
library(dplyr)
library(stringi)
library(stringr)
library(limma)
library(edgeR)
library(GseaVis)
# devtools::install_github("junjunlab/GseaVis")

load('LSCC.mRNAmatrix.symbol.rda')


gene_exp <- gene_exp1
# gene_exp <- apply(gene_exp1,2,function(x){log2(x+1)})

min(gene_exp)
max(gene_exp)


gene <- "LNX1"
gene.exp <- as.numeric(gene_exp[gene,])
# gene.exp <- is.na(gene.exp)

label <- if_else(gene.exp < median(as.numeric(gene.exp)),0,1)

# label <- if_else(gene.exp < median(gene.exp), 0, 1)

group.low <- gene_exp[,label == 0]
group.high <- gene_exp[,label == 1]


# dim(group.low)

group <- cbind(group.high,group.low)

# write.csv(group,"group.NAP1L2.csv")

group1 <- t(group)


# group1 <- group1[-1,]

# group1$Type <- c(rep('High',31),rep('Low',30))


Type <- c(rep('High',64),rep('Low',64))

riskTime <-as.data.frame(cbind(Type,group1))
riskTime$Type



group_list= riskTime$Type %>% factor(.,levels = c("High","Low"),ordered = F)


design = model.matrix(~0+group_list, data=group_list)
colnames(design) <- levels(group_list)
rownames(design) <-rownames(riskTime)
design <- as.data.frame(design)

dat1 <-gene_exp
dat1 <- as.data.frame(dat1)


# dat1 <-gene_exp1[,-1]
max(dat1)
min(dat1)
dat1 <- dat1[,rownames(design)]

contrast.matrix <- makeContrasts(High-Low, levels = design)#Case比Control
fit <- lmFit(dat1, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
options(digits = 4)#输出结果的小数点后保留4位


allDiff=topTable(fit2,coef=1,number=Inf)
allDiff <- na.omit(allDiff)

allDiff <- allDiff[order(allDiff$logFC,decreasing=T),]

gene <- allDiff$logFC
names(gene) <-rownames(allDiff)

#2、GO
#gmt 文件

c5 <- read.gmt("F:/workspace/GSEA_gmt/c5.go.v7.4.symbols.gmt")
gsea <- GSEA(gene, TERM2GENE=c5, verbose=FALSE, pvalueCutoff = 0.05); head(gsea)
write.table(gsea,'./06.GSEA/gsea.LNX1.go.txt',sep = '\t',quote = F,row.names = F)
dim(gsea)

pdf('./06.GSEA/GSEA.LNX1.GO.pdf',width = 9,height = 7)
# combine
gseaNb(object = gsea,
       geneSetID = gsea@result$ID[1:10],
       subPlot = 2,
       rmHt = T,
       addPval = T,
       pvalX = 1,pvalY = 1.3,
       termWidth = 50,
       # newGsea = T,
       # legend.position = c(0.7,0.8)
       curveCol = c("#96ceb4","#FFF000","#a39391","#e79686","#d9534f","#ffad60","#57D1C9","#a696c8","#2470a0","#d5a4cf")
)

dev.off()
# 
# pdf('GSEA.CXCR4.GO.pdf',width = 10,height = 10)
# gseaplot2(gsea, geneSetID = 1:10,pvalue_table = FALSE,
#           rel_heights = c(1.5, 0.5, 0.5),color= brewer.pal(10,'Paired'))
# dev.off()

#The ridgeplot will visualize expression distributions of core 
# #enriched genes for GSEA enriched categories. It helps users to interpret up/down-regulated pathways.
# pdf('GSEA.MARCD3.GO.ridge.pdf',width = 13,height = 6)
# ridgeplot(gsea,label_format = 60,showCategory = 10)
# dev.off()

#KEGG
c2 <- read.gmt("F:/workspace/GSEA_gmt/c2.cp.kegg.v7.4.symbols.gmt")
gsea <- GSEA(gene, TERM2GENE=c2, verbose=FALSE, pvalueCutoff = 0.05); head(gsea)
write.table(gsea,'./06.GSEA/gsea.LNX1.kegg.txt',sep = '\t',quote = F,row.names = F)
dim(gsea)

pdf('./06.GSEA/GSEA.LNX1.KEGG.pdf',width = 9,height = 7)
# combine
gseaNb(object = gsea,
       geneSetID = gsea@result$ID[1:10],
       subPlot = 2,
       rmHt = T,
       addPval = T,
       pvalX = 1,pvalY = 1.3,
       termWidth = 50,
       # newGsea = T,
       # legend.position = c(0.7,0.8)
       curveCol = c("#96ceb4","#FFF000","#a39391","#e79686","#d9534f","#ffad60","#57D1C9","#a696c8","#2470a0","#d5a4cf")
)

dev.off()
# 
# pdf('GSEA.CXCR4.KEGG.pdf',width = 10,height = 10)
# gseaplot2(gsea, geneSetID = 1:10,pvalue_table = FALSE,
#           rel_heights = c(1.5, 0.5, 0.5),color= brewer.pal(10,'Paired'))
# dev.off()

#The ridgeplot will visualize expression distributions of co

subPlot = 2,
termWidth = 30)

gseaNb(object = gsea,
       geneSetID = 'GOBP_RESPONSE_TO_TYPE_I_INTERFERON')




# 根据通路名称循环画显著的前4条
terms <- c('GOBP_RESPONSE_TO_TYPE_I_INTERFERON',
           'GOBP_DEFENSE_RESPONSE_TO_BACTERIUM',
           'GOBP_ANTIMICROBIAL_HUMORAL_RESPONSE',
           'GOBP_DEFENSE_RESPONSE_TO_VIRUS')

# plot
lapply(terms, function(x){
  gseaNb(object = gsea,
         geneSetID = x,
         newGsea = T,
         addPoint = F,
         addPval = T,
         pvalX = 0.75,pvalY = 0.75,
         pCol = 'black',
         pHjust = 0,
         newHtCol = c("red","white", "blue"),subPlot = 2,
         termWidth = 30)
}) -> gseaList



pdf('GSEA.CXCR4.GO.pdf',width = 10,height = 10)
# combine
cowplot::plot_grid(plotlist = gseaList,ncol = 2,align = 'hv')
dev.off()




# CIBERSORT ---------------------------------------------------------------



library(GSVA)
library(estimate)
library(dplyr)
library(tidyr)
require(ggplot2)
library(ggpubr)
library(pheatmap)
library(corrplot)
library(limma)

# 
# gene_exp <- read.table('cluter2.txt',
#                        sep = '\t',check.names = F,header = F)
# 
# ad=abbbs[gene_exp$V1,]
# 
# write.csv(ad,file = "cluter2.csv")


load("LSCC.mRNAmatrix.symbol.rda")
# resdt <- apply(gene_exp1,2,function(x){log2(x+1)})

resdt <- gene_exp1
resdt<-na.omit(resdt)
# resdt <- resdt[,meta$id]




write.table(resdt,file="exp.CIBERSORT.txt",sep="\t",row.names=T,quote=T)

res.tpm=normalizeBetweenArrays(resdt)
set.seed(2022)
start_time=Sys.time()
source("F:/workspace/CIBERSORT.R")
options(stringsAsFactors = F)
results=CIBERSORT("F:/workspace/ref.txt",###对比文件
                  "exp.CIBERSORT.txt",####"exp.CIBERSORT.txt",
                  perm=100, QN=T)
end_time=Sys.time()
cat(end_time-start_time)
write.table(results, file="CIBERSORT.Results.txt", sep="\t", row.names=T, col.names=T, quote=F)

CIBERSORT.filter <- results[which(results[,ncol(results)-2] <= 0.6),1:(ncol(results)-3)]
write.table(cbind(rownames(CIBERSORT.filter),CIBERSORT.filter), file="CIBERSORT.filter.csv", sep=",", row.names=F, col.names=T, quote=F)
data=t(results)
write.csv(data,file="CIBERSORT.免疫细胞比例.csv")

meta <-  read.table('LSCC.sample.txt',header = T,sep = '\t',,row.names= 1)
meta$id<-rownames(meta)


data=as.data.frame(results[which(results[,(ncol(results)-2)]<= 1),c(1:(ncol(results)-3))])
data$id <-  rownames(data)

data12 <-merge(meta,data,by.x="id",by.y = "id")
rownames(data12) <- data12$id

data12 <- data12[order(data12$Type,decreasing = F),]

data22 <- data12[,-1]

write.table(data12, file="CIBERSORT.dbox.txt", sep="\t", row.names=T, col.names=T, quote=F)



dbox <- data22
dbox  <- gather(dbox,Immunocyte,Score,2:ncol(dbox))

colnames(dbox) <- c("Type","Immunocyte","Score")

# write.table(dbox, file="CIBERSORT.dboxlast.txt", sep="\t", row.names=T, col.names=T, quote=F)
dbox <-  read.table('CIBERSORT.dboxlast.txt',header = T,sep = '\t',,row.names= 1)



results<-  read.table('CIBERSORT.免疫细胞比例.txt',header = T,sep = '\t',,row.names= 1)



pdf('GSE3912.免疫细胞占比小提琴图.pdf',width = 14,height = 8)
ggplot(dbox, aes(x=Immunocyte, y=Score,fill=Type,)) + 
  geom_violin(trim=T,color="black",scale = "width") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  geom_boxplot(width=0.4,position=position_dodge(0.9))+ #绘制箱线图
  stat_summary(fun= median, geom = "point",
               size = 1.5, color = "white",position=position_dodge(0.95))+#stat_summary() 可以加均值/中位值等 用fun= mean/median加均值/中位值
  scale_fill_manual(values = c( "#004BFB","#F91F10","white"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=30,hjust = 1,colour="black",family="Times",size=15), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
        
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属性
                                 size=16),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=18),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),#不显示网格线
        panel.grid.minor = element_blank()+
          stat_compare_means(size = 4,aes(group=Type),
                             label = 'p.signif'))+  #显示显著性
  ylab("Proportion of cells")+xlab("Immunocyte")+ #设置x轴和y轴的标题
  
  stat_compare_means(size = 8,aes(group=Type),
                     label = 'p.signif',label.y = 0.7)
dev.off()



library(RColorBrewer)

# 画箱线图
pdf('免疫细胞占比箱线图.pdf',width = 14,height = 10)
cdl= c(brewer.pal(n = 12, name = "Set3"),brewer.pal(n = 12, name = "Paired"))


ggboxplot(dbox,x='Immunocyte',y='Score',fill ='Type',###按照类型就是fill ='Type'
          ylab = 'Score',
          xlab ='Immunocyte',palette = "jco",
          #add ='jitter',add.params = list(size=1),
          size =1)+
  rotate_x_text(45)+
  theme(axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=13),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        axis.line = element_line(size=1))+
  stat_compare_means(size = 8,aes(group=Type),
                     label = 'p.signif')+scale_fill_manual(
                       values = c("blue",
                                  "red",
                                  "#29AF73"),
                       guide = guide_legend(reverse = TRUE)
                     )
dev.off()



results <- results[meta$id,]

# 画堆积图
library(ggplot2)
library(tidyverse)
colour <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

my_theme <- function(){
  theme(panel.grid = element_blank(),       # 网格线
        panel.border = element_blank(),     # 面板边框
        legend.position="right",            # legend位置
        legend.text = element_text(size=8), # legend内容大小
        legend.title = element_text(size=8),# legend标题大小
        axis.line = element_line(size=1),   # 坐标轴线
        text = element_text(family="Times"),# 文本字体
        axis.text.y = element_text(size = 8,face='bold',color='black'),# y轴标签样式
        axis.text.x = element_text(size = 8,face='bold',color='black',angle=90,hjust=1),        # x轴标签样式，angle=45 倾斜 45 度
        axis.title = element_text(size=10,face="bold"),  # 轴标题
        plot.title = element_text(hjust=0.7,size=15))    # 距，设置绘图区域距离边的据类，上、右、下、左
}  

p1 <- results[,1:22] %>% reshape2::melt() %>%
  ggplot(aes(x=Var1,y=value,fill=Var2)) +
  geom_bar(stat='identity') +
  # coord_flip()  +#反转坐标系
  scale_fill_manual(values =colour ) +
  theme_bw()+ theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+my_theme()


pdf("Cluster.cibersort.pdf")
p1
dev.off()





#免疫细胞和关键基因相关性####
modelgene <- c("WDR54","KAT2B","NBEAL2","LNX1")
ssgseaOut <-data.frame(results[,1:22],check.names = F)
# ssgsea <- t(ssgseaOut[-1,])
ssgsea <- ssgseaOut

# expr <- gene_exp1

expr <- gene_exp1[modelgene,rownames(ssgsea)]



library(psych)
#批量计算相关性
immuscore <- function(gene){
  
  y <- as.numeric(expr[gene,])
  colnames <- colnames(ssgsea)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- corr.test(as.numeric(ssgsea[,x]), y , method="spearman")
    data.frame(gene=gene,immune_cells=x,cor=dd$r,p.value=dd$p )
  }))
}



#批量计算genelist跟免疫浸润相关性的结果
data <- do.call(rbind,lapply(modelgene,immuscore))
head(data)
#保存到文件
write.csv(data, "gene.ssGSEA.correlation.csv", quote = F, row.names = F)




# 绘制散点图
data[is.na(data)] <- 0
#增加一列，区分p值的大小,使用两个ifelse实现三分类
data$pstar <- ifelse(data$p.value < 0.05,
                     ifelse(data$p.value < 0.01,"**","*"),
                     "")


pdf('gene.ssGSEA.correlation.pdf',width = 9,height = 6)
# 配色
my_palette <- colorRampPalette(c("blue", "white", "red"), alpha=TRUE)(n=399)

ggplot(data, aes(x=immune_cells,y=gene)) +
  geom_point(aes(size=-log10(p.value),color=cor)) +
  scale_color_gradientn('Pearson Correlation', 
                        colors=my_palette) + # 用自定义颜色画点
  theme_bw() +
  geom_text(aes(label=pstar),col ="black",size = 7)+
  theme(panel.grid.minor = element_blank(), #不画网格
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title = element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black")) #边框

dev.off()




#增加一列，区分p值的大小,使用两个ifelse实现三分类
data$pstar <- ifelse(data$p.value < 0.05,
                     ifelse(data$p.value < 0.01,"**","*"),
                     "")

#相关性用颜色的不同来表示，相关性的大小用颜色的深浅来反映；
#有差异的把*号打印在热图上
pdf('gene.ssGSEA.correlation.pdf',width = 8,height = 4)
ggplot(data, aes(immune_cells, gene)) + 
  geom_tile(aes(fill = cor), colour = "black",size=1.5)+
  scale_fill_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c")+
  geom_text(aes(label=pstar),col ="black",size = 7)+
  # geom_text(aes(label=pstar),col ="black",size = 2.5,vjust=1.5)+
  theme_minimal()+# 不要背景
  theme(axis.title.x=element_blank(),#不要title
        axis.ticks.x=element_blank(),#不要x轴
        axis.title.y=element_blank(),#不要y轴
        axis.text.x = element_text(angle = 45, hjust = 1,size = 15),# 调整x轴文字
        axis.text.y = element_text(size = 15))+#调整y轴文字
  #调整legen
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))
dev.off()

# 语义相似性图 ------------------------------------------------------------------


######语义相似性图

library(org.Hs.eg.db)
library(GOSemSim)
library(reshape2)
library(ggplot2)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

id_gsym <- read.table('表达模式一致的基因.txt',
                      sep = '\t',check.names = T,header = T)

id_gsym$ENTREZID <- as.character(id_gsym$ENTREZID)
head(id_gsym)

#用godata()函数来构建相应物种的Molecular Function本体的GO DATA
mf <- godata('org.Hs.eg.db', ont="MF", computeIC = FALSE)
#用godata()函数来构建相应物种的Cellular Component本体的GO DATA
cc <- godata('org.Hs.eg.db', ont="CC", computeIC = FALSE)

#用mgeneSim来计算MF本体，基因之间的语义相似度，结果为一个行列相同的矩阵
simmf <- mgeneSim(id_gsym$ENTREZID, semData = mf, measure = "Wang", drop = NULL, combine = "BMA")
#用mgeneSim来计算CC本体，基因之间的语义相似度，结果为一个行列相同的矩阵
simcc <- mgeneSim(id_gsym$ENTREZID, semData = cc, measure = "Wang", drop = NULL, combine = "BMA")

#计算基因在MF本体和CC本体下的几何平均值，一个打分值同时包括基因的分子功能和细胞定位两个信息
fsim <- sqrt(simmf * simcc)


#将基因的名字由ENTREZID改为gene symbol，方便看懂。
colnames(fsim) = id_gsym$SYMBOL
rownames(fsim) = id_gsym$SYMBOL



#将基因自己和自己的相似度设为NA，方便接下来去掉。
for (i in 1:ncol(fsim)){
  fsim[i,i] <- NA
}

y <- melt(fsim) #把宽格式数据转化成长格式，其实就是把正方形矩阵转成三列
y <- y[!is.na(y$value),] #删掉带NA的行

# 把每两个基因之间的相似度保存到文件，只需要保存第一列基因名和第三列数值
write.csv(y[,c(1,3)], "GOSemSim.csv", row.names = F)

y <- read.csv("GOSemSim.csv")

head(y)

#计算每个基因跟其他基因相似度的平均值
y.mean <- aggregate(.~Var1,y,mean) 
m <- y.mean$value
names(m) <- y.mean$Var1
#按平均值给基因名排序，便于画图
y$Var1 <- factor(y$Var1, levels=names(sort(m)))

f <- function(y) {
  r <- quantile(y, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  r[3] <- mean(y)
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

pdf('语义相似性.pdf',width = 6,height = 6)
p1 <- ggplot(y, aes(Var1, value, fill = factor(Var1))) + 
  scale_fill_brewer(palette="Set3") + #配色
  guides(fill=FALSE) + #不显示图例
  
  stat_summary(fun.data= f, geom='boxplot') + 
  geom_hline(aes(yintercept=0.5), linetype="dashed") + #画一条虚线
  
  coord_flip() + # x、y坐标轴互换
  xlab("") + ylab("") + 
  theme(axis.text.x = element_text(family = "Arial", size = 16, face = "bold"),
        axis.text.y = element_text(family = "Arial", size = 16, face = "bold")) + 
  theme_bw() + 
  theme(panel.border=element_rect(size=1)) #边框粗细 
p1

dev.off()




#八、独立预后####


#load library
library(tidyr)
library(dplyr)
library(stringi)
library(stringr)
library(survival)
library(survminer)
library(glmnet)
library(survivalROC)
library(rms)


#把单因素和多因素中的基因换成临床性状即可
#需要注意的是，临床性状可以是字符串，也可以是数字
#若要把临床性状变成数字，则该临床性状必须是等级越高越严重，例如Stage，M分期，T分期，N分期等
#如果临床性状是亚型且亚型之间不存在严重程度的递进关系，则要保留其原始数据，不能改为数字

##input为独立预后的输入文件，每一个样本的所有使用的性状必须是完整的，否则多因素会报错
##通过与OSinput合并，保证每一个样本都具备基因表达量和临床性状信息
##我这里使用的示例数据是处理好之后的数据，这里展示处理过程
input <- read.table('Input8.survivalInput3.txt',header = T,sep = '\t',check.names = F,fill = NA)
# input <- merge(OSinput[,c('Patient','fullid')],input,by='Patient')
# rownames(input) <- input$Patient
#合并临床性状与风险评分
#这里所用的风险评分是单因素+LASSO，seed为2，最佳阈值分高低风险组
#此处未考虑该模型的有效性，仅将该模型风险评分用来示例独立预后分析
merge1 <- read.table('risk.txt',sep = '\t',header = T,check.names = F)
# vali <- read.table('test.risk.txt',sep = '\t',header = T,check.names = F)
merge1 <- rbind(train,vali)
#去掉多余的DFStime和DFStatus以及基因表达量，具体删除的列自己定
#以下用风险评分做独立预后的示例
#我这里选择把风险评分和高低风险组都保留在变量中以便随时切换
merge1 <- merge1[,-c(4:6,8)]
riskTime <- merge(merge1,input,by='fullid')

write.table(riskTime,file="riskTime.xls",sep="\t",row.names=F,quote=F)



riskTime <- read.table('riskTime.txt',header = T,sep = '\t',check.names = F,fill = NA)

#让临床信息排在后面，无关信息排在前面，方便单因素的循环
# indep_rt <- riskTime[,c(1,5,2:4,7:ncol(riskTime))]
indep_rt <- riskTime
# indep_rt <- indep_rt[,-5]
outTab2 = data.frame()


#独立预后-单因素分析
#此处用风险评分riskScore进行多因素分析，从第7列开始逐个因素进行分析
#如果风险评分在多因素分析中pvalue大于0.05，换成高低风险组（即样本是高风险组或低风险组的属性）
for(i in colnames(indep_rt[,4:ncol(indep_rt)])){
  cox <- coxph(Surv(futime, fustatus) ~ indep_rt[,i], data = indep_rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab2=rbind(outTab2,
                cbind(id=i,
                      HR=coxSummary$conf.int[,"exp(coef)"],
                      HR.95L=coxSummary$conf.int[,"lower .95"],
                      HR.95H=coxSummary$conf.int[,"upper .95"],
                      pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
outTab2[,c(2:5)] <- apply(outTab2[,c(2:5)],2,function(x){as.numeric(x)})
outTab2 <- outTab2[order(outTab2$pvalue),]
write.table(outTab2,file="indep.uniCox.xls",sep="\t",row.names=F,quote=F)




#独立预后-多因素分析
indep.uni <- outTab2[outTab2$pvalue<=0.05,]$id


indep_rt1 <- indep_rt[,c('futime','fustatus',unique(indep.uni))]
multiCox1=coxph(Surv(futime, fustatus) ~ ., data = indep_rt1)
multiCox1=step(multiCox1,direction = "both")
multiCoxSum1=summary(multiCox1)
outTab3=data.frame()
outTab3=cbind(
  HR=multiCoxSum1$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum1$conf.int[,"lower .95"],
  HR.95H=multiCoxSum1$conf.int[,"upper .95"],
  pvalue=multiCoxSum1$coefficients[,"Pr(>|z|)"])
outTab3=cbind(id=row.names(outTab3),outTab3)
outTab3=data.frame(outTab3)
if (length(row.names(outTab3)) <=1){
  next
}
outTab3[,c(2:5)] <- apply(outTab3[,c(2:5)],2,function(x){as.numeric(x)})
outTab3 <- outTab3[order(outTab3$pvalue,decreasing = F),]
indep.final <- outTab3[outTab3$pvalue<1,]$id
write.table(outTab3,'indep.mulCox.xls',row.names = F,col.names = T,sep = '\t',quote = F)

###或者这个绘制森林图

###单因素森林图
outTab3 <- read.table('indep.uniCox.xls',header = T,sep = '\t',check.names = F)
tabletext <- outTab3

#tabletext <-subset(outTab3, outTab3$pvalue < 0.05)


tabletext[,2:ncol(tabletext)] <- apply(tabletext[,2:ncol(tabletext)],2,
                                       function(x){format(as.numeric(x),digits=4)})
tabletext <- rbind(c("Variable","HR","lower 95%CI","upper 95%CI","pvalue"),tabletext)
pdf('Fig.uniForest.pdf',width = 12,height = 7)
forestplot(tabletext,
           mean=c(NA, outTab3$HR),#HR值
           lower=c(NA, outTab3$HR.95L), #95%置信区间下限
           upper=c(NA, outTab3$HR.95H),#95%置信区间上限
           graph.pos=6,#图在表中的列位置
           graphwidth = unit(.25,"npc"),#图在表中的宽度比
           fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),#box颜色
           boxsize=0.4,#box大小固定
           lwd.ci=3,
           ci.vertices.height = 0.2,ci.vertices=T,#不显示区间
           zero=1,#zero线横坐标
           lwd.zero=2,#zero线宽
           xticks = c(0,1,5,10,15,20,25,30),#横坐标刻度根据需要可随意设置
           lwd.xaxis=2,
           xlab="HR",
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),#第一行顶部加黑实线
                           "2" = gpar(lwd=1, col="grey50", lty=2),#第二行顶部加灰色虚线
                           "15" = gpar(lwd=2, col="black")),#最后一行底部加黑线，""中数字为nrow(tabletext) + 1
           txt_gp=fpTxtGp(label=gpar(cex=1.2),#各种字体大小设置
                          ticks=gpar(cex=1.2),
                          xlab=gpar(cex=1.2),
                          title=gpar(cex=1.2)),
           lineheight = unit(1.5,"cm"),#固定行高
           colgap = unit(0.5,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
           
)
dev.off()

###多因素森林图
outTab3 <- read.table('indep.mulCox.xls',header = T,sep = '\t',check.names = F)

tabletext <-subset(outTab3, outTab3$pvalue < 1)


tabletext[,2:ncol(tabletext)] <- apply(tabletext[,2:ncol(tabletext)],2,
                                       function(x){format(as.numeric(x),digits=4)})
tabletext <- rbind(c("Variable","HR","lower 95%CI","upper 95%CI","pvalue"),tabletext)
pdf('Fig.mulForest.pdf',width = 10,height = 4)
forestplot(tabletext,
           mean=c(NA, outTab3[outTab3$pvalue <1,]$HR),#HR值
           lower=c(NA, outTab3[outTab3$pvalue <1,]$HR.95L), #95%置信区间下限
           upper=c(NA, outTab3[outTab3$pvalue <1,]$HR.95H),#95%置信区间上限
           graph.pos=6,#图在表中的列位置
           graphwidth = unit(.25,"npc"),#图在表中的宽度比
           fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),#box颜色
           boxsize=0.4,#box大小固定
           lwd.ci=3,
           ci.vertices.height = 0.2,ci.vertices=T,#不显示区间
           zero=1,#zero线横坐标
           lwd.zero=2,#zero线宽
           xticks = c(0,1,5,10,15,20),#横坐标刻度根据需要可随意设置
           lwd.xaxis=2,
           xlab="HR",
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),#第一行顶部加黑实线
                           "2" = gpar(lwd=1, col="grey50", lty=2),#第二行顶部加灰色虚线
                           "4" = gpar(lwd=2, col="black")),#最后一行底部加黑线，""中数字为nrow(tabletext) + 1
           txt_gp=fpTxtGp(label=gpar(cex=1.2),#各种字体大小设置
                          ticks=gpar(cex=1.2),
                          xlab=gpar(cex=1.2),
                          title=gpar(cex=1.2)),
           lineheight = unit(2.3,"cm"),#固定行高
           colgap = unit(0.5,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
           
)
dev.off()



nomod <- read.table('汇总1.txt',header = T,sep = '\t',check.names = F,fill = NA)
# nomod <- indep_rt
nomod <- nomod[,c(1:4,6:9)]

# nomod <- indep_rt[,c(1:6,9,11)]

# 绘制打分曲线
rocCol <- rainbow(7)
aucText <- c()

pdf("比较.第五年.cliROC.pdf",width = 8,height =8)

par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime =nomod$futime,status = nomod$fustatus,marker = nomod$riskScore,predict.time = 365*5,method = "KM" )###三年预测期
plot(roc$FP,roc$TP,type = "l",xlim = c(0,1),ylim = c(0,1),col=rocCol[1],
     xlab = "False positive rate",ylab = "True positive rate",
     lwd=3,cex.main=1.3,cex.lab=1.2,cex.axis=1.2,font=1.2,thresholds="best",print.thres="best")
aucText <- c(aucText,paste0("riskScore for 5−Year (AUC=",sprintf("%.2f",roc$AUC),")"))
abline(0,1)

# 
# 
# plot(roc$FP,roc$TP,type = "l")
# plot(rr,thresholds="best",print.thres="best")
# 

# 绘制其余ROC
j <- 1
for(i in colnames(nomod[,5:ncol(nomod)])){
  roc=survivalROC(Stime =nomod$futime,status = nomod$fustatus,marker = nomod[,i],predict.time =365*5,method = "KM" )###三年预测期
  j <- j+1
  lines(roc$FP,roc$TP,type = "l",xlim = c(0,1),ylim = c(0,1),col=rocCol[j],lwd=3)
  aucText <- c(aucText,paste0(i," for 5−Year (AUC=",sprintf("%.2f",roc$AUC),")"))
}
legend("bottomright",aucText,lwd=2,bty="n",col=rocCol)
dev.off()





#十二、森林图及风险曲线，确定种子后再画####
#单因素森林图，找到合适的种子再画####

library("survival")
library("survminer")
library("forestplot")

# tabletext <- outTab

tabletext <- outTab[outTab$pvalue <0.1,]


tabletext[,2:ncol(tabletext)] <- apply(tabletext[,2:ncol(tabletext)],2,
                                       function(x){format(as.numeric(x),digits=4)})
tabletext <- rbind(c("Variable","HR","lower 95%CI","upper 95%CI","pvalue"),tabletext)
pdf('Fig3.unicox.pdf',width = 10,height = 5)
forestplot(tabletext, mean=c(NA, outTab[outTab$pvalue <0.1,]$HR),#HR值
           lower=c(NA, outTab[outTab$pvalue <0.1,]$HR.95L), #95%置信区间下限
           upper=c(NA, outTab[outTab$pvalue <0.1,]$HR.95H),#95%置信区间上限
           graph.pos=6,#图在表中的列位置
           graphwidth = unit(.25,"npc"),#图在表中的宽度比
           fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),#box颜色
           boxsize=0.4,#box大小固定
           lwd.ci=3,
           ci.vertices.height = 0.2,ci.vertices=T,#不显示区间
           zero=1,#zero线横坐标
           lwd.zero=2,#zero线宽
           xticks = c(0,1,2,3,4),#横坐标刻度根据需要可随意设置
           lwd.xaxis=2,
           xlab="HR",
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),#第一行顶部加黑实线
                           "2" = gpar(lwd=1, col="grey50", lty=2),#第二行顶部加灰色虚线
                           "7" = gpar(lwd=2, col="black")),#最后一行底部加黑线，""中数字为nrow(tabletext) + 1
           txt_gp=fpTxtGp(label=gpar(cex=1.2),#各种字体大小设置
                          ticks=gpar(cex=1.3),
                          xlab=gpar(cex=1.2),
                          title=gpar(cex=1.2)),
           
           lineheight = unit(1.8,"cm"),#固定行高
           colgap = unit(0.5,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F)

dev.off()





rt=read.table("cox.result_1.txt",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
#定义颜色
forestCol="red"
clrs=fpColors(box=forestCol, line="darkblue", summary="royalblue")
#定义图片文字
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )
#绘制森林图
pdf('sub.multiCox_2.pdf',width=9, height=4, onefile=FALSE)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=5,
           boxsize=0.2,
           title="Overall Survival",
           xlab="Hazard ratio",
           txt_gp=fpTxtGp(ticks=gpar(cex=1.1), xlab=gpar(cex = 1.25))
)
dev.off()






############绘制森林图函数
bioForest=function(coxFile=null, forestFile=null, forestCol=null){
  #读取输入文件
  rt=read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  data=as.matrix(rt)
  HR=data[,1:3]
  hr=sprintf("%.3f",HR[,"HR"])
  hrLow=sprintf("%.3f",HR[,"HR.95L"])
  hrHigh=sprintf("%.3f",HR[,"HR.95H"])
  pVal=data[,"pvalue"]
  pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
  #定义颜色
  clrs=fpColors(box=forestCol, line="darkblue", summary="royalblue")
  #定义图片文字
  tabletext <- 
    list(c(NA, rownames(HR)),
         append("pvalue", pVal),
         append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )
  #绘制森林图
  pdf(file=forestFile, width=9, height=6, onefile=FALSE)
  forestplot(tabletext, 
             rbind(rep(NA, 3), HR),
             col=clrs,
             graphwidth=unit(50, "mm"),
             xlog=T,
             lwd.ci=4,
             boxsize=0.6,
             title="Overall survival",
             xlab="Hazard ratio",
             txt_gp=fpTxtGp(ticks=gpar(cex=1.1), xlab=gpar(cex = 1.25))
  )
  dev.off()
}

#调用函数,绘制森林图
bioForest(coxFile="cox.result.txt", forestFile="forest.pdf", forestCol="red")






# 多因素森林图####
tabletext <- outTab1
tabletext[,2:ncol(tabletext)] <- apply(tabletext[,2:ncol(tabletext)],2,
                                       function(x){format(as.numeric(x),digits=4)})
tabletext <- rbind(c("Variable",'coef',"HR","lower 95%CI","upper 95%CI","pvalue"),tabletext)
pdf('Fig3.multicox.pdf',width = 10,height =7)

forestplot(tabletext,
           mean=c(NA,outTab1$HR),#HR值
           lower=c(NA,outTab1$HR.95L), #95%置信区间下限
           upper=c(NA,outTab1$HR.95H),#95%置信区间上限
           graph.pos=7,#图在表中的列位置
           graphwidth = unit(.25,"npc"),#图在表中的宽度比
           fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),#box颜色
           boxsize=0.3,#box大小固定
           lwd.ci=3,
           ci.vertices.height = 0.2,ci.vertices=T,#不显示区间
           zero=1,#zero线横坐标
           lwd.zero=2,#zero线宽
           xticks = c(0,0.5,1,1.5,2),#横坐标刻度根据需要可随意设置
           lwd.xaxis=2,
           xlab="HR",
           hrzl_lines=list("1" = gpar(lwd=2, col="black"),#第一行顶部加黑实线
                           "2" = gpar(lwd=1, col="grey50", lty=2),#第二行顶部加灰色虚线
                           "7" = gpar(lwd=2, col="black")),#最后一行底部加黑线，""中数字为nrow(tabletext) + 1
           txt_gp=fpTxtGp(label=gpar(cex=1.2),#各种字体大小设置
                          ticks=gpar(cex=1.3),
                          xlab=gpar(cex=1.2),
                          title=gpar(cex=1.2)),
           lineheight = unit(2,"cm"),#固定行高
           colgap = unit(0.5,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
dev.off()

# 批量找miRNA、lncRNA关系 -------------------------------------------------------

library(BiocManager)
#读取及基本信息

library("multiMiR")

gene_exp1 <- read.table('lncRNA.txt',
                        sep = '\t',check.names = F,header = F)



# x <- c("Cd40","Pax8","Tslp","Timp1","Rorc"
#        ,"Il1rl1","Gata6","Serpine1")
# 
# 
# x <- c("Fcnaos","Lockd","Lockd")
# 


x <- gene_exp1[,1]

gene2mir <- get_multimir(org     = 'hsa', #小鼠mmu 人hsa
                         target  = x,
                         table   = 'validated',
                         summary = TRUE,
                         predicted.cutoff.type = 'n',
                         predicted.cutoff      = 500000)

table(gene2mir@data$database)

ez = gene2mir@data[gene2mir@data$database=="tarbase",];dim(ez)

write.table(gene2mir@data,file="circRNA.txt",sep="\t",row.names = T,quote = F)






gene_exp1 <- read.table('06.ceRNA/miRNA-mRNA.txt',
                        sep = '\t',check.names = F,header = T)
x <- gene_exp1[,1]
# x <- gene_exp1$mature_mirna_acc

gene2mir <- get_multimir(org     = 'hsa',
                         mirna = x,
                         table = "validated",
                         summary = TRUE, 
                         predicted.cutoff.type = 'n',
                         predicted.cutoff      = 500000)
table(gene2mir@data$database)

ez = gene2mir@data[gene2mir@data$database=="tarbase",];dim(ez)

write.table(gene2mir@data,file="miRNA-lncRNA.txt",sep="\t",row.names = T,quote = F)









# ENCORI数据库获取 -------------------------------------------------------------

load("F:\\workspace\\starBase.Rdata")

gene_exp1 <- read.table('09.ceRNA/miRNA-mRNA.txt',
                        sep = '\t',check.names = F,header = T)



data12 <-merge(gene_exp1,mirna_lncrna_dabase,by.x="mirnaid",by.y = "miRNAname")

write.table(data12,file="miRNA-lncRNA.ENCORI12.txt",sep="\t",row.names = T,quote = F)
























