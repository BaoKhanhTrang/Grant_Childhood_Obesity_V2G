# eFigure 8: The pathfindR method
combine = readRDS("pathway/pathfindR/combine.rds")
exp_num = readRDS("exp_num.rds")
module_df = readRDS(file="WGCNA_module_df.rds")


enrich = melt(combine[,c(2,grep("Fold_Enrichment",colnames(combine)))],
              id.vars = "Term_Description",variable.name = "Celltype",value.name = "Fold_Enrichment")
enrich$Celltype = sub("Fold_Enrichment_","",enrich$Celltype)
enrich = enrich[!is.na(enrich$Fold_Enrichment),]
enrich$p=0
enrich$num_gene=0
for(i in 1:nrow(enrich)){
  enrich$p[i] = combine[which(combine$Term_Description==enrich$Term_Description[i]),which(colnames(combine) == paste0("lowest_p_",enrich$Celltype[i]))]
  
  enrich$num_gene[i] = length(unique(unlist(strsplit(unlist(combine[which(combine$Term_Description==enrich$Term_Description[i]),
                          grep(paste0("regulated_",enrich$Celltype[i]),colnames(combine))]),split=", "))))
}
enrich = merge(enrich,allstat[,c('Celltype','Name')],by="Celltype")
enrich$Name = factor(enrich$Name,levels = levels(exp_num$Name))
enrich$Term_Description = factor(enrich$Term_Description,levels = rev(paste(final$Term_Description)))

trio=NULL
for(i in 1:nrow(enrich)){
  ge=combine[which(combine$Term_Description==enrich$Term_Description[i]),
          which(colnames(combine)==paste0("Up_regulated_",enrich$Celltype[i]))]
  up = unique(unlist(strsplit(ge,split=", ")))
  ge=combine[which(combine$Term_Description==enrich$Term_Description[i]),
             which(colnames(combine)==paste0("Down_regulated_",enrich$Celltype[i]))]
  
  down = unique(unlist(strsplit(ge,split=", ")))

  add=data.frame(Term_Description=enrich$Term_Description[i],
                 Celltype=enrich$Celltype[i],
                 gene_name=c(up,down),
                 UpDown=c(rep("up",length(up)),rep("down",length(down))))
  trio = rbind(trio,add)
}
trio = merge(trio,module_df,by="gene_name")  
trio = unique(trio)
module=unique(module_df$module)
enr=dcast(formula = Term_Description+Celltype~"num_gene",
          value.var = "gene_name",fun.aggregate = length,
          data=trio)

enr = merge(enr,enrich[,c("Term_Description","Celltype","p","Name")],by=c("Term_Description","Celltype"))    
enr = merge(enr,allstat[,c('Name','system')],by="Name")
enr$Name = factor(enr$Name,levels = levels(exp_num$Name))
enr$Term_Description = factor(enr$Term_Description,levels = rev(paste(final$Term_Description)))

g1=ggplot(enr)+geom_point(aes(Name,Term_Description,size=factor(num_gene),color=-log10(p)))+
    facet_grid(cols = vars(system),drop = T,space = "free",scales = "free")+
  scale_size_manual(name="# genes",values = c(1,2,3))+theme_bw()+
  theme(axis.text.x = element_blank(),axis.title = element_blank(),
        legend.position = "left")+    
  guides(color = guide_colorbar(show.limits = TRUE))

ge=unique(trio[,c('gene_name','Celltype')])
ge = merge(ge,unique(exp_num[,c('gene_name','Celltype','Name','value','system')]),by=c('gene_name','Celltype'))
ge$Name = factor(ge$Name,levels(exp_num$Name))
ge$gene_name = factor(ge$gene_name,levels(exp_num$gene_name))
ge=merge(ge,module_df,by="gene_name")

g2=ggplot(ge)+geom_tile(aes(Name,gene_name,fill=value))+theme_bw()+
  scale_fill_gradient2()+
  facet_grid(cols = vars(system),rows = vars(module),scales = "free",space = "free",drop = T)+
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust = 0.5),
        strip.text.x = element_blank(),strip.text.y = element_blank(),
        legend.position = "left")

gee=dcast(gene_name+Term_Description~"num_cell",data=trio[,c('gene_name','Celltype','Term_Description')],value.var = "Celltype")
gee$Term_Description = factor(gee$Term_Description,levels = rev(final$Term_Description))
gee=merge(gee,module_df,by="gene_name")
gee$gene_name = factor(gee$gene_name,levels = rev(levels(exp_num$gene_name)))

g3=ggplot(gee)+geom_tile(aes(gene_name,Term_Description,fill=factor(num_cell)))+
  facet_grid(cols=vars(module),scales = "free",space = "free",drop = T)+
  geom_text(aes(gene_name,Term_Description,label=num_cell),color="white")+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),axis.text.y = element_blank())+
  scale_fill_identity()

pdf("CB_EGG.topLD_pathfindR.pdf",width = 10,height = 10)
g1 + g3 +g2 +plot_layout(ncol = 2,heights = c(10,3),widths = c(3,2))
dev.off()