# eFigure 10: The GnRH signaling pathway
library(pathview)

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


tmp = unique(trio[which(trio$Term_Description=="GnRH signaling pathway"),c('gene_name','Celltype')])
tmp = merge(tmp,ge[,c('gene_name','Celltype','Name','value')])
vec=dcast(gene_name~Name,data=tmp,value.var = 'value')
tmp <- AnnotationDbi::mget(paste(vec$gene_name),
        AnnotationDbi::revmap(org.Hs.eg.db::org.Hs.egSYMBOL),
        ifnotfound = NA)

tmp=as.data.frame(t(as.data.frame(tmp)))
colnames(tmp)="EGID"
tmp$gene_name=rownames(tmp)
vec = merge(vec,tmp,by="gene_name")
rownames(vec)=vec$EGID
tmp = data.frame(vec[,-which(colnames(vec) %in% c('gene_name','EGID'))],row.names = vec$EGID)
colnames(tmp)=colnames(vec)[-c(1,ncol(vec))]

pathview(gene.data = tmp,pathway.id = "hsa04912",multi.state=TRUE,
           out.suffix = "_pathfindR",same.layer=F,both.dirs = list(gene=F),
           limit = list(gene=c(0,90)),
           kegg.dir = "pathway/Kegg/")