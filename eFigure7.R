# eFigure 7: Gene Ontology (GO) biological process terms enrichment

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggsci)
result = function(do,showcat){
  ke.result = as.data.frame(do@compareClusterResult)
  ke.result$GeneRatio = sapply(ke.result$GeneRatio, function(x) eval(parse(text=x)))
  ke.result$geneSet = round(as.numeric(ke.result$Count / ke.result$GeneRatio))
  bg=strsplit(ke.result$BgRatio[1],split="/")[[1]][[2]]
  ke.result$BgRatio = sub(paste0("/",bg),"",ke.result$BgRatio)
  ke.result$BgRatio = as.numeric(ke.result$BgRatio)
  ke.result$Cluster=paste0(ke.result[,2]," (",ke.result$geneSet,")")
  orde=unique(paste(ke.result$Cluster))
  ke.result$Cluster = factor(ke.result$Cluster,levels = orde)
  ke.result$Description=paste0(ke.result$Description," (",ke.result$BgRatio,")")
  ord2=paste(unique(paste(ke.result$Description)))
  ke.result$Description = factor(ke.result$Description,levels = ord2)
  cluster = unique(ke.result[,2])
  cat=NULL
  for(i in 1:length(cluster)){
    cat= union(cat,ke.result$Description[which(ke.result[,2]==cluster[i])][1:showcat])
  }
  ke.result = ke.result[which(ke.result$Description %in% cat),]
  if(length(which(is.na(ke.result$qvalue)))!=0){ke.result = ke.result[-which(is.na(ke.result$qvalue)),]}
  return(ke.result)
}

exp_num = readRDS("exp_num.rds")
allv2g=readRDS("allv2g.rds")


genes=merge(exp_num[,c('Name','gene_name','Scaled\nexp.level','chr','system')],
            allv2g[,c('gene_name','gene_id')],by="gene_name")
genes=unique(genes)

ge=bitr(unique(genes$gene_name), fromType = "SYMBOL",
           toType =  "ENTREZID",
           OrgDb = org.Hs.eg.db)
colnames(ge)[1]=c("gene_name")
genes=merge(genes,ge,by="gene_name",all.x=TRUE)
genes=genes[order(genes$`Scaled\nexp.level`,decreasing = TRUE),]

# Panel A
dgn=compareCluster(gene_name~Name,data=unique(genes[,c('gene_name','Name')]),
                   fun="enrichGO",keyType="SYMBOL",OrgDb=org.Hs.eg.db,ont = "BP")
dgn.result=result(dgn,5)
dgn.result = merge(dgn.result,allstat[,c('Name','system')],by="Name")
cluster = unique(dgn.result[,c('Name','Cluster')])
cluster$Name = factor(cluster$Name,levels = levels(exp_num$Name))
cluster = cluster[order(cluster$Name),]
dgn.result$Cluster = factor(dgn.result$Cluster,levels = paste(cluster$Cluster))


genes_in_path = unique(unlist(strsplit(dgn.result$geneID,split="/")))
genes_in_path = data.frame(matrix(0,ncol = length(unique(dgn.result$Description)),nrow = length(genes_in_path)))
rownames(genes_in_path) = unique(unlist(strsplit(dgn.result$geneID,split="/")))
colnames(genes_in_path)=unique(dgn.result$Description)
for(i in 1:ncol(genes_in_path)){
  genes_in_path[unique(unlist(strsplit(dgn.result$geneID[which(dgn.result$Description == colnames(genes_in_path)[i])],split="/"))),i]=1
}
ord=colSums(genes_in_path)
names(ord)=colnames(genes_in_path)
ord=sort(ord,decreasing = T)
genes_in_path = genes_in_path[,names(ord)]

ord2=rowSums(genes_in_path)
names(ord2) = rownames(genes_in_path)
ord2=sort(ord2,decreasing = T)
genes_in_path = genes_in_path[names(ord2),]

g=ggplot(dgn.result)+
  geom_point(aes(y=Description,x=Cluster,color=qvalue,size=factor(Count)),alpha=.7)+theme_bw()+
  scale_color_gradient(low = "darkblue",high = "red")+
  scale_y_discrete(limits=rev(names(ord)))+
  scale_size_manual(name="# genes",values = c(3,4,5,6,7))+
  ylab("")+xlab("GO BP")+
  facet_grid(cols = vars(system),drop = TRUE,space = "free",scales = "free")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),legend.position = "left")+
  guides(color = guide_colorbar(show.limits = TRUE))

genes_in_path$gene = rownames(genes_in_path)
genes_in_path = melt(genes_in_path,id.vars = 'gene')
genes_in_path = genes_in_path[-which(genes_in_path$value==0),]
colnames(genes_in_path)[2]="Description"
for(i in 1:nrow(genes_in_path)){
  genes_in_path$value[i] = length(unique(dgn.result$Cluster[intersect(which(dgn.result$Description==paste(genes_in_path$Description[i])),grep(genes_in_path$gene[i],dgn.result$geneID))]))
}
genes_in_path$Description=factor(genes_in_path$Description,levels = rev(names(ord)))
genes_in_path$gene = factor(genes_in_path$gene,levels = names(ord2))
g2=ggplot(genes_in_path)+geom_tile(aes(gene,Description,fill=factor(value)))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
        axis.text.y=element_blank(),axis.title.y = element_blank())
library(egg)
pdf("GO_BP_celltypes_5each_regroup.pdf",width = 14,height = 5.5)
ggarrange(g,g2,ncol=2,widths = c(2,1))
dev.off()

# Panel B
do=compareCluster(gene_name~system,data=unique(genes[,c('gene_name','system')]),
                  fun="enrichGO",keyType="SYMBOL",OrgDb=org.Hs.eg.db,ont = "MF",pvalueCutoff=0.1)
do.result=result(do,5)
g=ggplot(do.result)+
  geom_point(aes(y=Description,x=Cluster,color=qvalue,size=factor(Count)),alpha=.7)+theme_bw()+
  scale_color_gradient(low = "darkblue",high = "red")+
  xlab("")+ylab("GO molecular function")+
  scale_size_manual(name="# genes",values = c(3,4,5,6,7))+
  guides(color = guide_colorbar(show.limits = TRUE))
pdf("GO_MF_systems.pdf",height = 4,width = 6.6)
print(g)
dev.off()