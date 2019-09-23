# WGCNA 思考题

---

<table><tr><td bgcolor=#D1EEEE>
title: "WGCNA 思考题"
author: "dawen"
date: "2019/8/30"
</td></tr></table>

```{r warning=FALSE, message=FALSE}
rm(list = ls())
library(tidyr)
library(tidyverse)
library(stringr)
library(VennDiagram)
library(WGCNA)
library(knitr)
load(file = "../02_data/TCGA_TNBC_pair_annotation_result.Rdata")
## step1_构建WGCNA所需表达矩阵
# --------------------------------------------------------------------------------

exprSet_lncRNA = as.data.frame(t(exprSet_lncRNA))
exprSet_mRNA = as.data.frame(t(exprSet_mRNA))
exprSet_mRNA$sample = rownames(exprSet_mRNA)
exprSet_lncRNA$sample = rownames(exprSet_lncRNA)
tmp = merge(exprSet_mRNA,exprSet_lncRNA, by = "sample")

rownames(tmp) = tmp$sample
tmp = tmp[,-1]

# -----------------------------------------------------------------------------------
group_list=factor(ifelse(as.numeric(substr(rownames(exprSet_mRNA),14,15)) < 10,'tumor','normal'))
table(group_list)
datTraits <- data.frame(row.names=rownames(exprSet_mRNA),
                        subtype=group_list)

rt = tmp
datExpr = rt
nSamples = nrow(datExpr)


# step2_找到合适的beta值
#---------------------------------------------------------------------------------------
powers1=c(seq(1,10,by=1),seq(12,30,by=2))
RpowerTable=pickSoftThreshold(datExpr, powerVector=powers1)[[2]]
cex1=0.7
#pdf(file="softThresholding.pdf")
par(mfrow=c(1,2))
plot(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n")
text(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2], labels=powers1,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
plot(RpowerTable[,1], RpowerTable[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n")
text(RpowerTable[,1], RpowerTable[,5], labels=powers1, cex=cex1,col="red")
beta1=14


# step3_对基因聚类成模块并且可视化
#---------------------------------------------------------------------------------------
net = blockwiseModules(
  datExpr,
  power = 14,
  maxBlockSize = 6000,
  TOMType = "unsigned", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "AS-green-FPKM-TOM",
  verbose = 3
)
table(net$colors) 
mergedColors = labels2colors(net$colors)
table(mergedColors)
# 将目标基因聚成几个模块
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


# step4_将样本信息添加进去,分析样本性状与模块的关系
# --------------------------------------------------------------------------
design=model.matrix(~0+ datTraits$subtype)
design = as.data.frame(design)
colnames(design)=levels(datTraits$subtype)
moduleColors <- labels2colors(net$colors)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, design , use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
```

---

<font face="黑体" color=black size=5>WGCNA主要是能够将相似的基因聚类成一个module，然后探讨不同module对于对于样本性状的影响，从上图可以看出brown,turquoise与normal的性状正相关，与另外两个模块则相关性降低。肿瘤组织的形状则刚好相反，差异基因是tumor vs.normal 正常推论来说上调基因应该与brown,turquoise交集少，下调基因则反之。</font>

---



###  step7_提取每个模块的基因

```{r warning=FALSE, message=FALSE}
# --------------------------------------------------------------------------------------
# Select module 
Select_module <- function(x){
  
  # Select module probes
  probes = colnames(datExpr) ## 我们例子里面的probe就是基因名
  inModule = (mergedColors==x);
  blue_modProbes = probes[inModule];
  
}
# brown_modProbes
# -------------------------------------------------------------------------------------
brown_modProbes = Select_module("brown")
length(brown_modProbes)
head(brown_modProbes)


# blue_modProbes
# --------------------------------------------------------------------------------------
blue_modProbes = Select_module("blue")
length(blue_modProbes)
head(blue_modProbes)

# grey_modProbes
# --------------------------------------------------------------------------------------
grey_modProbes = Select_module("grey")
length(grey_modProbes)
head(grey_modProbes)

# turquoise_modProbes
# --------------------------------------------------------------------------------------
turquoise_modProbes = Select_module("turquoise")
class(turquoise_modProbes)
length(turquoise_modProbes)


# step8_从差异基因里选出上调基因和下调基因
# --------------------------------------------------------------------------------------
load(file = "../02_data/TCGA_TNBC_pair_annotation_result.Rdata")
logFC_cutoff = 2
DEG = mRNA_exprSet
DEG = na.omit(DEG)

DEG$change = as.factor( ifelse( DEG$padj< 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
                                ifelse( DEG$log2FoldChange > logFC_cutoff , 'UP', 'DOWN' ), 'NOT' ) )
table(DEG$change)
up_gene = DEG[DEG$change == "UP",]
down_gene = DEG[DEG$change == "DOWN",]
# 进一步处理方便做venn图
# ---------------------------------------------------------------------------------
up_gene = as.data.frame(up_gene$gene_id %>%
                          str_split("\\|",simplify = T))
names(up_gene) = c("symbol", "probe_id","gene_type")

down_gene = as.data.frame(down_gene$gene_id %>%
                            str_split("\\|",simplify = T))
names(down_gene) = c("symbol", "probe_id","gene_type")
```



### step9_差异基因与各个颜色模块基因的交集

```{r warning=FALSE, message=FALSE}
library("clusterProfiler")
library("org.Hs.eg.db")
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(enrichplot)
#library(plyr)
library(ggrepel)
# 提取不同模块与差异基因交集的大部分基因
# ----------------------------------------------------------------------------------------
# brown_modProbes
brown_down_gene = down_gene[down_gene$symbol %in% brown_modProbes,]
dim(brown_down_gene)


# turquoise
turquoise_down_gene = down_gene[down_gene$symbol %in% turquoise_modProbes,]
dim(turquoise_down_gene)


#grey
grey_up_gene = up_gene[up_gene$symbol %in% grey_modProbes,]
dim(grey_up_gene )


#blue
blue_up_gene = up_gene[up_gene$symbol %in% blue_modProbes,]
dim(blue_up_gene )
# 用venn 图可视化
# ----------------------------------------------------------------------------------
### brown_modProbes
grid.newpage()
venn.plot <- venn.diagram(list(up_gene = up_gene$symbol, #画图
                               brown_modProbes = brown_modProbes), NULL, 
                          fill = c("#E31A1C","#E7B800"), 
                          alpha = c(0.5,0.5), cex = 4, cat.fontface=3, 
                          category.names = c("up_gene", "brown_modProbes"), 
                          main = "Overlap")
grid.draw(venn.plot)

grid.newpage()
venn.plot <- venn.diagram(list(down_gene = down_gene$symbol, #画图
                               brown_modProbes = brown_modProbes), NULL, 
                          fill = c("#E31A1C","#E7B800"), 
                          alpha = c(0.5,0.5), cex = 4, cat.fontface=3, 
                          category.names = c("down_gene", "brown_modProbes"), 
                          main = "Overlap")
grid.draw(venn.plot)
```

```{r warning=FALSE, message=FALSE}
# 对提取的基因进行注释
# brown
# -----------------------------------------------------------------------------------------
gene.df <- bitr(brown_down_gene$symbol, fromType="SYMBOL",
                toType="ENTREZID", 
                OrgDb = "org.Hs.eg.db")
go <- enrichGO(gene = gene.df$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all")
barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free")


ego_KEGG <- enrichKEGG(gene = gene.df$ENTREZID, organism = "hsa", 
                       keyType = "kegg",
                       pvalueCutoff = 0.9,
                       pAdjustMethod = "BH",
                       minGSSize = 10, maxGSSize = 500, 
                       qvalueCutoff = 0.9,
                       use_internal_data = FALSE)
dotplot(ego_KEGG, showCategory = 20,font.size = 8)
```





---

<font face="黑体" color=black size=5> brown大多数基因与下调基因重合</font>

---



```{r warning=FALSE, message=FALSE}
### turquoise_modProbes
# --------------------------------------------------------------------------------
grid.newpage()
venn.plot <- venn.diagram(list(up_gene = up_gene$symbol, #画图
                               turquoise_modProbes = turquoise_modProbes), NULL, 
                          fill = c("#E31A1C","#E7B800"), 
                          alpha = c(0.5,0.5), cex = 4, cat.fontface=3, 
                          category.names = c("up_gene", "turquoise_modProbes"), 
                          main = "Overlap")
grid.draw(venn.plot)

grid.newpage()
venn.plot <- venn.diagram(list(down_gene = down_gene$symbol, #画图
                               turquoise_modProbes = turquoise_modProbes), NULL, 
                          fill = c("#E31A1C","#E7B800"), 
                          alpha = c(0.5,0.5), cex = 4, cat.fontface=3, 
                          category.names = c("down_gene", "turquoise_modProbes"), 
                          main = "Overlap")
grid.draw(venn.plot)
```


```{r warning=FALSE, message=FALSE}
# turquoise
# -----------------------------------------------------------------------------------------
gene.df <- bitr(turquoise_down_gene$symbol, fromType="SYMBOL",
                toType="ENTREZID", 
                OrgDb = "org.Hs.eg.db")
go <- enrichGO(gene = gene.df$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all")
barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free")


ego_KEGG <- enrichKEGG(gene = gene.df$ENTREZID, organism = "hsa", 
                       keyType = "kegg",
                       pvalueCutoff = 0.9,
                       pAdjustMethod = "BH",
                       minGSSize = 10, maxGSSize = 500, 
                       qvalueCutoff = 0.9,
                       use_internal_data = FALSE)
dotplot(ego_KEGG, showCategory = 20,font.size = 8)
```



----

<font face="黑体" color=black size=5> turquoise大多数基因与下调相关</font>

---



```{r warning=FALSE, message=FALSE}
### blue_modProbes
# --------------------------------------------------------------------------------
grid.newpage()
venn.plot <- venn.diagram(list(up_gene = up_gene$symbol, #画图
                               blue_modProbes = blue_modProbes), NULL, 
                          fill = c("#E31A1C","#E7B800"), 
                          alpha = c(0.5,0.5), cex = 4, cat.fontface=3, 
                          category.names = c("up_gene", "blue_modProbes"), 
                          main = "Overlap")
grid.draw(venn.plot)

grid.newpage()
venn.plot <- venn.diagram(list(down_gene = down_gene$symbol, #画图
                               blue_modProbes = blue_modProbes), NULL, 
                          fill = c("#E31A1C","#E7B800"), 
                          alpha = c(0.5,0.5), cex = 4, cat.fontface=3, 
                          category.names = c("down_gene", "blue_modProbes"), 
                          main = "Overlap")
grid.draw(venn.plot)

```

```{r warning=FALSE, message=FALSE}
# blue
# -----------------------------------------------------------------------------------------
gene.df <- bitr(blue_up_gene$symbol, fromType="SYMBOL",
                toType="ENTREZID", 
                OrgDb = "org.Hs.eg.db")
go <- enrichGO(gene = gene.df$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all")
barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free")


ego_KEGG <- enrichKEGG(gene = gene.df$ENTREZID, organism = "hsa", 
                       keyType = "kegg",
                       pvalueCutoff = 0.9,
                       pAdjustMethod = "BH",
                       minGSSize = 10, maxGSSize = 500, 
                       qvalueCutoff = 0.9,
                       use_internal_data = FALSE)
dotplot(ego_KEGG, showCategory = 20,font.size = 8)
```



---

<font face="黑体" color=black size=5> blue大多数与上调基因重合，与tumor性状相关</font>

---



```{r warning=FALSE, message=FALSE}
### grey_modProbes
# --------------------------------------------------------------------------------
grid.newpage()
venn.plot <- venn.diagram(list(up_gene = up_gene$symbol, #画图
                               grey_modProbes = grey_modProbes), NULL, 
                          fill = c("#E31A1C","#E7B800"), 
                          alpha = c(0.5,0.5), cex = 4, cat.fontface=3, 
                          category.names = c("up_gene", "grey_modProbes"), 
                          main = "Overlap")
grid.draw(venn.plot)

grid.newpage()
venn.plot <- venn.diagram(list(down_gene = down_gene$symbol, #画图
                               grey_modProbes = grey_modProbes), NULL, 
                          fill = c("#E31A1C","#E7B800"), 
                          alpha = c(0.5,0.5), cex = 4, cat.fontface=3, 
                          category.names = c("down_gene", "grey_modProbes"), 
                          main = "Overlap")
grid.draw(venn.plot)
```


```{r warning=FALSE, message=FALSE}

# grey
# -----------------------------------------------------------------------------------------
gene.df <- bitr(grey_up_gene$symbol, fromType="SYMBOL",
                toType="ENTREZID", 
                OrgDb = "org.Hs.eg.db")
go <- enrichGO(gene = gene.df$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all")
barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free")


ego_KEGG <- enrichKEGG(gene = gene.df$ENTREZID, organism = "hsa", 
                       keyType = "kegg",
                       pvalueCutoff = 0.9,
                       pAdjustMethod = "BH",
                       minGSSize = 10, maxGSSize = 500, 
                       qvalueCutoff = 0.9,
                       use_internal_data = FALSE)
dotplot(ego_KEGG, showCategory = 20,font.size = 8)
```

---

<font face="黑体" color=black size=5> grey大多数与上调基因相关</font>

---

---

<table><tr><td bgcolor=#C0FF3E><font face="黑体" color=black size=5>综上所诉，WGCNA分析通过对基因的聚类分成不同的模块可以识别不同性状所特有的基因模块，可以结合差异分析进一步筛选，也可以单独进行，尤其适合某一性状具有多种分类的情况，此时可用wgcna来寻找不同性状的特定的一类基因。</font></td></tr></table>
---

