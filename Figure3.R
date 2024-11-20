library(reshape2)
library(ggplot2)
library(ggsci)
library(paletteer)


allv2g=readRDS("allv2g.rds")
allstat=readRDS("allstat.rds")
metaData=readRDS("metaData.rds")
allExp=readRDS(file="allExp.rds")

pir = resmelt(dcast(data=unique(allv2g[which(allv2g$OCRtype=="PIR"),c('proxy_rsid','gene_name','Celltype')]),gene_name~Celltype,fun.aggregate = length))
colnames(pir)[3] = "PIR"
pir = pir[-which(pir$PIR==0),]
prom = melt(dcast(data=unique(allv2g[which(allv2g$OCRtype=="Prom"),c('proxy_rsid','gene_name','Celltype')]),gene_name~Celltype,fun.aggregate = length))
colnames(prom)[3] = "Prom"
prom = prom[-which(prom$Prom==0),]
colnames(prom)[2]=colnames(pir)[2]="Celltype"
prom$pair=paste(prom$gene_name,prom$Celltype,sep = "_")
pir$pair=paste(pir$gene_name,pir$Celltype,sep = "_")
prom$PIR=0
pir$Prom=0
exp_num = merge(pir,prom,by="pair",all=TRUE)
exp_num$gene_name = exp_num$gene_name.x
exp_num$gene_name[which(is.na(exp_num$gene_name.x))] = exp_num$gene_name.y[which(is.na(exp_num$gene_name.x))]
exp_num$Celltype = apply(exp_num[,c('Celltype.x','Celltype.y')],1,function(x) paste(unique(na.exclude(x))))
exp_num$PIR=rowSums(exp_num[,c('PIR.x','PIR.y')],na.rm = TRUE)
exp_num$Prom=rowSums(exp_num[,c('Prom.x','Prom.y')],na.rm = TRUE)
exp_num = exp_num[,10:13]
rm(pir,prom)


scaled = allExp
scaled[is.na(scaled)]=0
for(i in 2:ncol(scaled)){
  scaled[,i] = scales::rescale(x = log2(unlist(allExp[,i])),to = c(0, 100))
}
scaled = scaled[-which(duplicated(scaled$gene_name)),]
rownames(scaled)=scaled$gene_name
scaled=scaled[which(scaled$gene_name %in% unique(allv2g$gene_name)),]
tmp = melt(scaled,id.vars = 'gene_name')
colnames(tmp)[2]="id"
tmp = merge(tmp,metaData[,c('id','Celltype')],by="id")
tmp = tmp[which(tmp$Celltype %in% allstat$Celltype),-1]
scaled = dcast(gene_name~Celltype,data=tmp,value.var = 'value',fun.aggregate = function(x) mean(x))
tmp = melt(scaled,id.vars = 'gene_name',variable.name = "Celltype")

exp_num = merge(exp_num,tmp,by=c("gene_name","Celltype"),all.x=TRUE)
colnames(exp_num)[ncol(exp_num)]="Scaled\nexp.level"
exp_num$gene_name=factor(exp_num$gene_name,levels = levels(allv2g$gene_name))
exp_num$`Scaled\nexp.level`=ceiling(exp_num$`Scaled\nexp.level`)
exp_num = merge(exp_num,allv2g[,c('gene_name','chr')],by="gene_name")
exp_num=unique(exp_num)
exp_num = merge(exp_num,unique(allstat[,c('Celltype','Name','system')]),by="Celltype")
exp_num=unique(exp_num)
exp_num$chr = factor(exp_num$chr,levels = unique(allv2g$chr))
exp_num$`Scaled\nexp.level`[which(is.na(exp_num$`Scaled\nexp.level`))]=0
exp_num$`Scaled\nexp.level`[which(is.infinite(exp_num$`Scaled\nexp.level`))]=0
exp_num$color=paste(exp_num$PIR,exp_num$Prom,sep=",")
exp_num = exp_num[order(exp_num$PIR,exp_num$Prom,decreasing = F),]
exp_num$shape=2
exp_num$shape[which(exp_num$PIR==0)]=15
exp_num$shape[which(exp_num$Prom==0)]=19
exp_num$shape[which(exp_num$PIR!=0 & exp_num$Prom !=0)]=18
exp_num$shape[which(exp_num$`Scaled\nexp.level`==0)] = 2

mycol=sample(paletteer_d("ggsci::default_igv"),size = length(unique(exp_num$color)),replace = F)
names(mycol)=paste(unique(exp_num$color))
for(i in 1:nrow(exp_num)){
  exp_num$mycolor[i]=mycol[[exp_num$color[i]]] 
}
exp_num=merge(exp_num,unique(allv2g[,c('gene_name','Name','locus')]),by=c("gene_name","Name"))
exp_num$locus = factor(exp_num$locus, levels = levels(proxy$locus))


legend = unique(exp_num[,c('mycolor','color','shape')])                
legend=legend[order(legend$color),]
tmp= unlist(strsplit(legend$color,split=","))
legend$x = as.numeric(paste(tmp[seq(1,length(tmp),by=2)]))
legend$y = as.numeric(paste(tmp[seq(2,length(tmp),by=2)]))
legend = legend[order(legend$x),]
legend$shape[which(legend$x==0)]=15
legend$shape[which(legend$y==0)]=19
legend$shape[which(legend$x!=0 & legend$y !=0)]=18
legend$x = factor(legend$x,levels = rev(unique(legend$x)))
pdf("legend_cell_gene.pdf",height = 3,width = 2)
ggplot(legend) + geom_point(aes(y,x,color=color,fill=color,shape=shape),size=5,alpha=.6)+
  scale_color_manual(name="Number of\nacting SNPs\n#PIR,#Promoter",
                     values = deframe(unique(exp_num[,c('color','mycolor')]))) +
  scale_fill_manual(name="Number of\nacting SNPs\n#PIR,#Promoter",
                     values = deframe(unique(exp_num[,c('color','mycolor')]))) +
  scale_shape_identity()+
  scale_x_continuous(name="#Promoter-SNPs",breaks=unique(legend$y),limits=c(-0.1,5.2))+
  scale_y_discrete(name="#PIR-SNPs")+
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid.minor = element_blank(),
        legend.position = "none")
dev.off()

gene_cell=dcast(unique(allv2g[,c('gene_name','Name')]),Name~gene_name,fun.aggregate = length)
gene_cell$gene_num=apply(gene_cell[,-1],1,function(x) length(which(x!=0)))
gene_cell = gene_cell[,c(1,ncol(gene_cell))]
gene_cell = gene_cell[order(gene_cell$gene_num,decreasing = F),]
gene_cell$Name = factor(gene_cell$Name,levels = paste( gene_cell$Name))
gene_cell=merge(gene_cell,allstat[,c('Name','system')],by="Name")
gene_cell = unique(gene_cell)

g2=ggplot(gene_cell)+geom_bar(aes(gene_num,Name),stat = "identity")+
  facet_grid(rows=vars(system),drop=TRUE,scales = "free", space = "free")+
  scale_x_continuous(expand = c(0,0),name="# implicated\ngenes",breaks = c(0,25,50),limits = c(0,50))+
  theme_bw()+
  theme(panel.grid.minor.x = element_blank(),panel.grid.major.y = element_blank(),
        axis.text.y = element_blank(),axis.title=element_blank())


exp_num$Name = factor(exp_num$Name,levels = levels(gene_cell$Name))
g=ggplot(exp_num)+
  geom_point(aes(gene_name,Name,size=`Scaled\nexp.level`,color=color,shape=shape),alpha=.75)+
  theme_bw()+
  scale_size_binned(breaks = seq(0,80,by=15))+
  scale_color_manual(name="Number of\nacting SNPs\n#PIR,#Promoter",
                     values = deframe(unique(exp_num[,c('color','mycolor')]))) +
  scale_shape_identity() +
  facet_grid(rows=vars(system),cols=vars(locus),drop=TRUE,scales = "free", space = "free")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust=0.5,size = 7),axis.title = element_blank(),
        panel.grid.major = element_line(colour = "grey70",linetype = "dotted",linewidth = 0.1),
        panel.grid.minor = element_blank(),
        strip.text = element_blank(),
        legend.position = "none",
        legend.title = element_text(hjust=0.5),
        legend.margin = margin(t=0,r=-2,b=0,l=0, unit = "lines"),
        legend.justification = "center" )+
  guides(size = guide_bins(show.limits = TRUE))


cell_gene=dcast(unique(allv2g[,c('gene_name','Name','system')]),
                gene_name~system,value.var = "Name",fun.aggregate = length)
cell_gene$sum=rowSums(cell_gene[,-1])
tmp=dcast(unique(allv2g[,c('gene_name','system')]),gene_name~system)
tmp$system=apply(tmp[,-1],1,function(x) paste(x[!is.na(x)],collapse=","))
tmp$system[which(tmp$system=="immune,metabolic,neural,other")] = "all"
tmp$system[which(tmp$system=="immune,metabolic,neural")] = "all_main"
tmp=merge(tmp,cell_gene[,c('gene_name','sum')],by="gene_name")
cell_gene = cell_gene[,-ncol(cell_gene)]
cell_gene = melt(cell_gene,id.vars = "gene_name")
cell_gene = cell_gene[-which(cell_gene$value==0),]
colnames(cell_gene)[2:3]=c("system","cell_num")
cell_gene$gene_name = factor(cell_gene$gene_name,levels = levels(exp_num$gene_name))
tmp$gene_name = factor(tmp$gene_name,levels = levels(exp_num$gene_name))
cell_gene = unique(cell_gene)
cell_gene=merge(cell_gene,exp_num[,c('gene_name','locus')],by="gene_name")
tmp=merge(tmp,exp_num[,c('gene_name','locus')],by="gene_name")
cell_gene = unique(cell_gene)
cell_gene$locus = factor(cell_gene$locus ,levels = levels(proxy$locus))
tmp$locus = factor(tmp$locus ,levels = levels(proxy$locus))


g3=ggplot(cell_gene)+geom_bar(aes(gene_name,cell_num,fill=system),stat = "identity")+
  geom_point(data=tmp[which(tmp$system =="all"),],aes(gene_name,sum+2),color="red",shape="*",size=4)+
  geom_point(data=tmp[which(tmp$system =="all_main"),],aes(gene_name,sum+2),color="blue",shape="*",size=4)+
  scale_y_continuous(expand = c(0,0),breaks = seq(0,50,by=10),limits = c(0,max(tmp$sum)+max(tmp$sum)/5))+
  facet_grid(cols=vars(locus),drop=TRUE,scales = "free", space = "free")+
  theme_bw()+
  scale_fill_manual(name="",values = c("#222222","#a9a9a9","#545454","#d4d4d4"))+
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(), 
        strip.text.x=element_text(angle=90,hjust=0.5),
        axis.text.x = element_blank(),axis.title=element_blank(),
        legend.position = "none")
pdf("CB_EGG.topLD_celltype_vs_genes.pdf",height = 9.6,width = 15)
g3 + plot_spacer()+g+g2+plot_layout(widths = c(13, 1),heights = c(1,6))
dev.off()
