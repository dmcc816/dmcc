
rm(list=ls())
library(clusterProfiler)
library(org.Rn.eg.db)
#Hs人的数据库，大鼠org.Rn.eg.db，小鼠 org.Mm.eg.db
library(dplyr)
library(DOSE)
library(ggplot2)
library(tidyr)
library(forcats)
# 路径设置
data <- read.csv("ZCF_vs_Vehicle.csv") # 数据读取
#将Symbol转换为EntrezID
entrez <-bitr(
  data$gene,
  fromType = 'SYMBOL',
  toType = 'ENTREZID',
  OrgDb = 'org.Rn.eg.db'
)

head(entrez)
id = unique(as.character(entrez[,2]))


######## KEGG富集分析 #########
kk <- enrichKEGG(gene = id,
                 organism = "rno",#小鼠mmu,大鼠rno
                 keyType = 'kegg',#entrez/kegg
                 pvalueCutoff = 1,#可以自定义cutoff
                 pAdjustMethod = "fdr",
                 qvalueCutoff = 1)#可以自定义cutoff
# geneID转换
kk=setReadable(kk, OrgDb = org.Rn.eg.db , keyType="ENTREZID")

#提取KEGG富集结果表格作图
KEGG_result<-kk@result
KEGG_result <- separate(data=KEGG_result, 
                        col=GeneRatio,
                        into = c("GR1","GR2"),
                        sep = "/")
KEGG_result <- mutate(KEGG_result, 
                      GeneRatio = (as.numeric(GR1)/as.numeric(GR2)))
KEGG_result <- KEGG_result %>% arrange(desc(Count))#将term顺序按照Count进行排序
plot_data <- KEGG_result[1:10,]#提取前15个pathway

# 定义您感兴趣的 ID 列表
specific_ids <- c("hsa04022", "hsa04024", "hsa04062", "hsa04660", "hsa04657", 
                  "hsa04064", "hsa04630", "hsa04659", "hsa04060", "hsa04623")  # 替换为实际的 ID 值
# 根据 ID 列提取行
plot_data <- KEGG_result[KEGG_result$ID %in% specific_ids, ]

##### KEGG柱状图
ggplot(plot_data,aes(x=Count,y=fct_reorder(Description, Count)))+#Count可以改为specific_ids排序
  geom_col(aes(fill=p.adjust))+
  facet_grid(scale='free_y',space = 'free_y')+
  #修改柱状图颜色
  scale_fill_gradient(low='#E27371',high = '#5D82A7')+
  #标题修改
  labs(title='KEGG Enrichment',
       y='Patway',
       x='Count')+
  guides(fill=guide_colorbar(reverse = T))+
  theme_bw()#主题
# 将 Description 转换为因子，并去掉排序，恢复默认顺序


##### KEGG富集分析气泡图
ggplot(plot_data,aes(x=GeneRatio,y=fct_reorder(Description, GeneRatio)))+
  geom_point(aes(size=GeneRatio,fill=pvalue),
             shape=21,
             color='black')+
  facet_grid(scale='free_y',space = 'free_y')+
  scale_fill_gradient(low='#E27371',high = '#5D82A7')+   #修改气泡图颜色
  scale_size_continuous(range = c(2,6)) +  # 这里设置气泡大小的范围
  #标题修改
  labs(title='KEGG Enrichment',
       y='Pathway',
       x='Gene Ratio')+
  guides(fill=guide_colorbar(reverse = T,order=1))+
  theme_bw()#主题

#保存结果
write.csv(KEGG_result,"KEGG_result.csv")




######## GO富集分析  ########
go_enrich<-clusterProfiler::enrichGO(gene = id,
                                     ont = 'all',#可选'BP','CC','MF' or 'all'
                                     keyType = "ENTREZID",
                                     OrgDb = org.Hs.eg.db,
                                     pAdjustMethod = "BH",#p值矫正方法
                                     pvalueCutoff = 1,
                                     qvalueCutoff = 1)
# geneID转换
go_enrich=setReadable(go_enrich, OrgDb = org.Hs.eg.db , keyType="ENTREZID")

#提取go富集结果表格 保存
go_result<-go_enrich@result

# 使用 str_wrap 函数处理描述标签，使其在指定字符数后换行
#go_geo@result$Description <- stringr::str_wrap(go_geo$Description, width = 60)
go_result <- separate(data=go_result, 
                      col=GeneRatio,
                      into = c("GR1","GR2"),
                      sep = "/")
go_result <- mutate(go_result, 
                    GeneRatio = (as.numeric(GR1)/as.numeric(GR2)))
result_BP<-go_result%>%filter(ONTOLOGY=='BP')
result_CC<-go_result%>%filter(ONTOLOGY=='CC')
result_MF<-go_result%>%filter(ONTOLOGY=='MF')
#取前5行
BP<-result_BP[1:5,]
CC<-result_CC[1:5,]
MF<-result_MF[1:5,]
all<-rbind(BP,CC,MF)

##### GO富集分析气泡图
ggplot(all,aes(x=GeneRatio,y=fct_reorder(Description, GeneRatio)))+#将term顺序按照GeneRatio进行排序
  #分面
  geom_point(aes(size=GeneRatio,fill=p.adjust),
             shape=21,
             color='black')+
  facet_grid(ONTOLOGY~.,
             scale='free_y',
             space = 'free_y'
  )+
  #修改气泡图颜色
  scale_fill_gradient(low='#E27371',high = '#5D82A7')+
  #标题修改
  labs(title='GO Enrichment',
       y='Pathway',
       x='Gene Ratio')+
  guides(fill=guide_colorbar(reverse = T,order=1))+
  theme_bw()#主题


#### GO富集分析柱状图
ggplot(all,aes(x=Count,
               y=fct_reorder(Description, Count)))+#将term顺序按照Count进行排序
  #分面
  geom_col(aes(fill=p.adjust))+
  facet_grid(ONTOLOGY~.,
             scale='free_y',space = 'free_y'
  )+
  #修改柱状图颜色
  scale_fill_gradient(low='#E27371',high = '#5D82A7')+
  #标题修改
  labs(title='GO Enrichment',
       y='Pathway',
       x='Count')+
  guides(fill=guide_colorbar(reverse = T))+
  theme_bw()#主题

#保存结果
write.csv(go_result,"go_result.csv")



######## GSEA富集分析########
library(ggplot2) #画图使用
library(clusterProfiler) #GSEA富集/数据读取使用
library(GSEABase) #GSEA富集使用
library(dplyr) #数据处理使用
library(data.table) #数据读取使用

# 注意GSEA富集分析需要所有基因，不必筛选差异基因，需要两列数据，gene、log2FC
setwd("F:\\AXB\\lession_test\\4_Funtion\\GSEA")
#读取 gmt文件 
geneset<-read.gmt("Human_KEGG.gmt")
#读取marker基因文件
data=read.csv("维恩富集分析.csv")

# 使用FDR作为排序依据作为备用，排序方式会影响显著性，log2FC结果不理想可以用这个
# 将对FDR值取对数后的结果乘以基因表达差异方向的正负号，从而获得一个综合的统计量"stat"作为排序依据。
data <- data %>% dplyr::mutate(stat = -log10(FDR) * sign(logFC))
geneList <- data$stat #获取GeneList
names(geneList) <- data$gene_id #使用转换好的ID，对GeneList命名
geneList <- sort(geneList, decreasing = T) #从高到低排序


# 直接使用log2FC作为排序依据
geneList <-  data$gene #获取GeneList
names(geneList) <- data$gene  #使用转换好的ID，对GeneList命名
geneList <- sort(geneList, decreasing = T) #从高到低排序
geneList[1:10]#检查一下
GSEA_enrichment <- GSEA(geneList, #排序后的gene
                        TERM2GENE = geneset, #基因集
                        pvalueCutoff = 0.05, #P值阈值
                        minGSSize = 10, #最小基因数量
                        maxGSSize = 5000, #最大基因数量
                        eps = 0, #P值边界
                        pAdjustMethod = "BH") #校正P值的计算方法

result <- data.frame(GSEA_enrichment)# 转化为数据框并保存结果
write.csv(result,"GSEA_result.csv")
dim(GSEA_enrichment@result)#查看富集通路数目
dotplot(GSEA_enrichment,showCategory=15,color="p.adjust")#查看前15通路

## 将通路分为激活和抑制两个部分
dotplot( GSEA_enrichment,split =  ".sign")+facet_grid(~.sign)+
  theme(plot.title = element_text(size = 10,color="black",hjust = 0.5),
        axis.title = element_text(size = 10,color ="black"), 
        axis.text = element_text(size= 10,color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1 ),
        legend.position = "top",
        legend.text = element_text(size= 10),
        legend.title= element_text(size= 10))
# 单通路展示
gseaplot2(GSEA_enrichment,"Herpes simplex virus 1 infection",color="red",pvalue_table = T)
# 多通路展示
gseaplot2(GSEA_enrichment,c("Antigen processing and presentation",
                            "Herpes simplex virus 1 infection"),
          color=c("red","blue"),pvalue_table = T)

# 山脊图
ridgeplot(GSEA_enrichment,
          showCategory = 10,#展示前10个通路
          fill = "p.adjust", #填充色 "pvalue", "p.adjust", "qvalue"
          core_enrichment = TRUE,#是否只使用 core_enriched gene
          label_format = 30,#通路字符超过30就换行
          orderBy = "NES",
          decreasing = T
)+
  theme(axis.text.y = element_text(size=8))


######## GSVA富集分析 #######

#####获取通路基因#####
library("KEGGREST")
library("EnrichmentBrowser") #这个包里面的一些函数会调用KEGGREST里面的函数

keggGet('K04687') 
gs<-keggGet('K04687')


#获取通路中gene信息 
gs[[1]]$GENE 
#查找所有基因 
genes<-unlist(lapply(gs[[1]]$GENE,function(x) strsplit(x,';'))) 
genelist <- genes[1:length(genes)%%3 ==2] 
genelist <- data.frame(genelist)  
#把结果写入表格中 
write.csv(genelist, "K04687.csv", row.names=FALSE)

###去批次
# 设置批次信息
batch <- pheno$batch # 批次
# 设置生物学分类，告诉函数不要把生物学差异整没了 
pheno$cancer <- factor(pheno$cancer, levels = c("Normal", "Cancer", "Biopsy"))
mod <- model.matrix(~as.factor(cancer), data=pheno)

expr_combat <- ComBat(dat = exprSet, batch = batch, mod = mod,par.prior=TRUE, ref.batch=1)
