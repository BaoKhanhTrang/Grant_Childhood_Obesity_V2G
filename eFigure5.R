# eFigure 5: The statistics of the 111 implicated genes

library(reshape2)
library(ggplot2)
allv2g=readRDS("allv2g.rds")


cell_gene=dcast(unique(allv2g[,c('gene_name','Name')]),gene_name~Name,fun.aggregate = length)
cell_gene$sum=rowSums(cell_gene[,-1])
cell_gene = cell_gene[,c('gene_name','sum')]
tmp=dcast(unique(allv2g[,c('gene_name','Name')]),gene_name~Name,value.var = "Name")
tmp$Names = apply(tmp[,-1],1,function(x) paste(na.exclude(x),collapse = ","))
cell_gene = merge(cell_gene,tmp[,c('gene_name','Names')],by="gene_name")
tmp=dcast(unique(allv2g[,c('gene_name','system')]),gene_name~system)
tmp$system=apply(tmp[,-1],1,function(x) paste(x[!is.na(x)],collapse=","))
tmp$system[which(tmp$system=="immune,metabolic,neural,other")] = "all"
tmp$system[which(tmp$system=="immune,metabolic,neural")] = "all_main"
cell_gene=merge(cell_gene,tmp,by="gene_name")
cell_gene=cell_gene[order(cell_gene$sum,decreasing = T),]
cell_gene=merge(cell_gene,exp_num[,c('gene_name','locus')],by="gene_name")
tmp=dcast(unique(allv2g[,c('gene_name','Name','system')]),
          gene_name~system,value.var = "Name",fun.aggregate = length)
cell_gene=merge(cell_gene[,c('gene_name','sum','Names','system','locus')],tmp,by="gene_name")
cell_gene = unique(cell_gene)
sort(table(cell_gene$Names[which(cell_gene$sum==1)]))

# Panel A 
tmp=as.data.frame(sort(table(cell_gene$sum)))
tmp = tmp[order(as.numeric(paste(tmp$Var1))),]
tmp$Var1 = factor(tmp$Var1,levels = rev(paste(tmp$Var1)))
pdf("Number of cell types vs number of genes.pdf",width = 3)
ggplot(tmp) + geom_bar(aes(Var1,Freq),stat = 'identity')+theme_bw()+
  ylab("Number of genes")+xlab("Number of cell types")+coord_flip()
dev.off()

# Panel B
tmp=as.data.frame(sort(table(cell_gene$system)))
pdf("Number of systems vs number of genes.pdf",width = 3,height = 3)
ggplot(tmp) + geom_bar(aes(factor(Var1),Freq),stat = 'identity')+theme_bw()+
  ylab("Number of genes")+xlab("System(s)")+coord_flip()
dev.off()
