# eFigure 9: modified SPIA analysis
# Panel B
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(patchwork)

dir="pathway/"
res=readRDS(paste0(dir,"SPIA/combine.rds"))
exp_num = readRDS("exp_num.rds")
allExp=readRDS(,"allExp.rds")
allstat=readRDS("allstat.rds")
module_df = readRDS(file="WGCNA_module_df.rds")

eg = bitr(paste(unique(allExp$gene_name)), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
colnames(eg)[1]="gene_name"

res$ID=apply(res,1,function(x) paste(strsplit(strsplit(x['KEGGLINK'],split="\\+")[[1]][[1]],split="\\?")[[1]][[2]]) )
res$ENTREZID=apply(res,1,function(x) paste(unlist(strsplit(x['KEGGLINK'],split="\\+"))[-1],collapse = ", ") )
res$genes = apply(res,1,function(x) paste(eg$SYMBOL[which(eg$ENTREZID %in% unlist(strsplit(x['ENTREZID'],split = ", ")))],collapse = ", ") )
cells=unique(res$Celltypes)

res.sig = res[which(round(res$pGFdr,digits = 2)<=0.05 | round(res$pGFWER,digits = 2) <= 0.05),]
res.sig = res.sig[order(res.sig$Name,res.sig$pGFdr,res.sig$pGFWER,decreasing = F),]
res.sig$num_gene = apply(res.sig,1,function(x) length(grep(",",paste(x['genes'])))+1)

colnames(res.sig)[1]="Term_Description"
res.sig = merge(res.sig,allstat[,c('Name','Celltype','system')],by="Celltype")
res.sig$Name = factor(res.sig$Name,levels = levels(exp_num$Name))
res.sig = res.sig[order(res.sig$pGFWER,decreasing = T),]
res.sig$Term_Description = factor(res.sig$Term_Description,levels = paste(unique(res.sig$Term_Description)))

g1=ggplot(res.sig)+
    geom_point(aes(Name,Term_Description,size=factor(num_gene),color=-log10(pGFWER)))+
    facet_grid(cols = vars(system),drop = T,space = "free",scales = "free")+
    scale_size_manual(name="# genes",values = c(1,2,3))+theme_bw()+
    theme(axis.text.x = element_blank(),axis.title = element_blank(),
    legend.position = "left")+    
    guides(color = guide_colorbar(show.limits = TRUE))

tmp=unique(res.sig[,c('genes','Name','num_gene')])
ge=NULL
for(i in 1:nrow(tmp)){
  gee = data.frame(Name=tmp$Name[i],gene_name=unlist(strsplit(tmp$genes[i],split=", ")))
  ge = rbind(ge,gee)
}
ge = merge(ge,unique(exp_num[,c('gene_name','Name','value','system')]),by=c('gene_name','Name'))
ge$Name = factor(ge$Name,levels(exp_num$Name))
ge$gene_name = factor(ge$gene_name,levels(exp_num$gene_name))
ge=merge(ge,module_df,by="gene_name")

g2=ggplot(ge)+geom_tile(aes(Name,gene_name,fill=value))+theme_bw()+
    scale_fill_gradient2()+
    facet_grid(cols = vars(system),rows = vars(module),scales = "free",space = "free",drop = T)+
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust = 0.5),
        strip.text.x = element_blank(),strip.text.y = element_blank(),
        legend.position = "left"
    )

tmp=unique(res.sig[,c('genes','Term_Description','Name')])
gee=NULL
for(i in 1:nrow(tmp)){
  gee = rbind(gee,data.frame(Term_Description=tmp$Term_Description[i],gene_name=unlist(strsplit(tmp$genes[i],split=", ")),Name=tmp$Name[i]))
}
gee=dcast(gene_name+Term_Description~"num_cell",data=tmp,value.var = "Name")
gee$Term_Description = factor(gee$Term_Description,levels = levels(res.sig$Term_Description))
gee=merge(gee,module_df,by="gene_name")
gee$gene_name = factor(gee$gene_name,levels = rev(levels(exp_num$gene_name)))

g3=ggplot(gee)+geom_tile(aes(gene_name,Term_Description,fill=factor(num_cell)))+
  facet_grid(cols=vars(module),scales = "free",space = "free",drop = T)+
  geom_text(aes(gene_name,Term_Description,label=num_cell),color="white")+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),axis.text.y = element_blank())+
  scale_fill_identity()

pdf("CB_EGG.topLD_SPIA.pdf",width = 12,height = 8)
g1 + g3 +g2 +plot_layout(ncol = 2,heights = c(8,3),widths = c(3,1))
dev.off()