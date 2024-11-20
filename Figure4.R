library(nVennR)
library(reshape2)
library(dendextend)

coloc=read.table("coloc/Childhood_obesity_eqtl_coloc_results_all_summary_2023-04-21_condPP4.txt",header = T)
tmp=unlist(strsplit(coloc$GeneID.Tissue,split="_"))
coloc$gene_id = apply(coloc,1,function(x) sub("_","",substring(paste(x[3]),first = 0,last = 18)))
coloc$Tissue = apply(coloc,1,function(x) sub(paste0(x[11],"_"),"",x[3]))
coloc=coloc[which(signif(coloc$cond.PP.H4.abf,1)>=0.8),]
allv2g=readRDS("allv2g.rds")

# Panel A ############################################################################################################
pairV2G = unique(allv2g[,c('proxy_rsid','gene_name')])
colnames(coloc)[1:2]=colnames(pairV2G)

simi=merge(unique(coloc),pairV2G,by=c('proxy_rsid','gene_name'))
simi = unique(simi)
simi = simi[which(signif(simi$cond.PP.H4.abf,1)>=0.8),]

genes=list(
  ColocQuiaL = unique(coloc$gene_name),
  V2G = unique(allv2g$gene_name),
  olap = intersect(unique(coloc$gene_name),unique(allv2g$gene_name)),
  Gene_variant_pairs = unique(simi$gene_name)
)

plotVenn(genes,labelRegions = F,opacity = 0.6,fontScale=4,
        setColors=c("#3fc061","#ef8910","#4a89ff","#dc0000"),
        outFile="Genes_coloc_vs_TopLD_v2g.svg"
)


# Panel B ############################################################################################################
###### create circos_genes.txt ####
genes = union(unique(allv2g$gene_name),unique(coloc$gene_name))
library(biomaRt)
library(org.Hs.eg.db)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
pos=getBM(mart = ensembl,attributes = c('chromosome_name','start_position','end_position','hgnc_symbol'),
          filters = 'hgnc_symbol',values = genes)
leftover=setdiff(genes,pos$hgnc_symbol)
prom=read.table("v2g/hg38_genecode.v30_promoter_reference.txt",header=TRUE)
pos2=merge(prom,allv2g[which(allv2g$gene_name %in% leftover),c('gene_id','gene_name')],by="gene_id")
pos2$start = pos2$pro_start + 1500
pos2$end = pos2$pro_end
pos2$start[which(pos2$strand=="-")] = pos2$pro_start[which(pos2$strand=="-")]
pos2$end[which(pos2$strand=="-")] = pos2$pro_end[which(pos2$strand=="-")] - 1500
leftover = unique(pos2$gene_name)
for(i in 1:length(leftover)){
  tmp = data.frame(chromosome_name = unique(pos2$pro_chr[which(pos2$gene_name == leftover[i])]),
                   start_position = min(pos2$start[which(pos2$gene_name == leftover[i])]),
                   end_position = max(pos2$end[which(pos2$gene_name==leftover[i])]),
                   hgnc_symbol=leftover[i])
  pos = rbind(pos,tmp)
}
leftover=setdiff(genes,pos$hgnc_symbol)
pos2 = getBM(mart = ensembl,attributes = c('chromosome_name','start_position','end_position','hgnc_symbol','ensembl_gene_id'),
             filters = 'ensembl_gene_id',values = leftover)

pos2$hgnc_symbol[which(pos2$hgnc_symbol=="")]=pos2$ensembl_gene_id[which(pos2$hgnc_symbol=="")]
pos=rbind(pos,pos2[,1:4])

pos$id="id=0"
genes = setdiff(allv2g$gene_name,coloc$gene_name)
pos$id[which(pos$hgnc_symbol %in% genes)] = "id=1"
genes = setdiff(intersect(unique(coloc$gene_name),unique(allv2g$gene_name)),unique(simi$gene_name))
pos$id[which(pos$hgnc_symbol %in% genes)] = "id=3"
pos$id[which(pos$hgnc_symbol %in% simi$gene_name)] = "id=2"
pos = pos[-which(pos$chromosome_name %in% c("HSCHR16_5_CTG1","HSCHR2_6_CTG7_2")),]
pos$chromosome_name = paste0("hs",gsub("chr","",pos$chromosome_name))
unique(pos$chromosome_name)
pos=pos[order(as.numeric(gsub("hs","",paste(pos$chromosome_name))),as.numeric(pos$start_position),decreasing = F),]

pos$mid=round((pos$end_position+pos$start_position)/2)
pos$start_position = pos$mid - 100
pos$end_position = pos$mid +100
pos=pos[,1:5]

rm(pos2,ensembl,prom,tmp,cs,genes,i,leftover,w)
###### create links for circos #######
snp_gene = unique(coloc[,c('snp','gene_name')])
colnames(snp_gene)[1]="proxy_rsid"
snp_gene = merge(snp_gene,proxy[,c('chr','start','variant_pos','proxy_rsid')],by="proxy_rsid")
snp_gene$chr = gsub("chr","hs",snp_gene$chr)
colnames(pos)[4]="gene_name"
snp_gene = merge(snp_gene,pos[,1:4],by="gene_name")
snp_gene = unique(snp_gene)
which(snp_gene$chromosome_name != snp_gene$chr)
write.table(snp_gene[,c('chr','start','variant_pos','chromosome_name','start_position','end_position')],
            file="circos/eqtl_links_v2.txt",
            col.names = F,row.names = F,sep="\t",quote=F)

pos=pos[-which(pos$gene_name %in% setdiff(pos$gene_name[which(pos$id=="id=0")],snp_gene$gene_name)),]
write.table(pos,file="circos_genes_v2.txt",
            col.names = F,row.names = F,sep="\t",quote = F)

###### create ci_links.txt ####
ocr_gene = unique(allv2g[,c('chr','ocr_start','ocr_end','gene_name')])
ocr_gene = merge(ocr_gene,pos[,1:4],by="gene_name")
ocr_gene$chr = gsub("chr","hs",ocr_gene$chr)
which(ocr_gene$chr != ocr_gene$chromosome_name)
write.table(ocr_gene[,2:7],file="circos/ci_links_v2.txt",
            col.names = F,row.names = F,sep="\t",quote = F)


###### create circos_rsID.txt ####
locus=unique(proxy[which(proxy$index_rsid %in% coloc$lead),c('chr','index_rsid','locus')])
for(i in 1:nrow(locus)){
  locus$start[i] = min(proxy$variant_pos[which(proxy$index_rsid==locus$index_rsid[i])])
  locus$end[i] = max(proxy$variant_pos[which(proxy$index_rsid==locus$index_rsid[i])])
}
locus$chr=gsub("chr","hs",locus$chr)
write.table(locus[,c('chr','start','end')],file="circos/highlights.txt",col.names = F,row.names = F,sep="\t",quote=F)


bed=pos[,1:4]
colnames(bed)=c('chr','start','end','gene_name')
tmp = unique(ocr_gene[,c('chr','ocr_start','ocr_end','gene_name')])
colnames(tmp)=c('chr','start','end','gene_name')
bed=rbind(bed,tmp)
tmp=snp_gene[,c('chr','start','variant_pos','gene_name')]
colnames(tmp)=c('chr','start','end','gene_name')
bed=rbind(bed,tmp)
bed = bed[order(as.numeric(gsub("hs","",bed$chr)),bed$start,bed$gene_name),]
tmp = bed
bed=NULL
for(i in 1:length(unique(tmp$gene_name))){
  a = tmp[which(tmp$gene_name == unique(tmp$gene_name)[i]),]
  a$start = min(a$start)
  a$end = max(a$end)
  a = unique(a)
  bed=rbind(bed,a)
}
for(i in 1:nrow(locus)){
  genes = union(unique(allv2g$gene_name[which(allv2g$index_rsid == locus$index_rsid[i])]),
                unique(coloc$gene_name[which(coloc$lead == locus$index_rsid[i])]))
  
  locus$out_start[i] = min(locus$start[i],pos$start_position[which(pos$gene_name %in% genes)])
  locus$out_end[i] = max(locus$end[i],pos$end_position[which(pos$gene_name %in% genes)])
}
locus = locus[order(as.numeric(gsub("hs","",locus$chr)),decreasing = F),]
write.table(locus[,c('chr','out_start','out_end')],file="circos/circos_regions.bed",
            col.names = F,row.names = F,sep="\t",quote=F
)
# circos plot was generated using Circos v0.68-1 download from https://circos.ca/
# by running 
# circos -conf circos_CO_topLD_coloc_vs_v2g.conf


# Panel C - ADCY3 coloc ############################################################################################################
proxy=readRDS("proxy.rds")
pair=coloc[which(coloc$snp %in% proxy$proxy_rsid[which(proxy$locus=="ADCY3")]),]
pair = pair[which(pair$gene_name=="ADCY3"),]
pair=pair[which(pair$credible95.y),]
colnames(pair)[which(colnames(pair)=="snp")]="proxy_rsid"
pair$pair=paste(pair$proxy_rsid,pair$gene_name,sep = "-")

pairV2G = merge(unique(pair[,c('proxy_rsid','gene_name')]),allv2g,by=c('proxy_rsid','gene_name'))
pairV2G$pair=paste(pairV2G$proxy_rsid,pairV2G$gene_name,sep = "-")
pairV2G = pairV2G[order(pairV2G$Name),]

pv=data.frame(matrix(0,ncol=length(unique(pair$Tissue)),nrow=length(unique(pairV2G$Celltype))))
rownames(pv) = unique(pairV2G$Celltype)
colnames(pv) = sort(unique(pair$Tissue))
for(i in 1:nrow(pv)){
  for(e in 1:ncol(pv)){
    pv[i,e] = length(intersect(pair$pair[which(pair$Tissue == colnames(pv)[e])],
                               pairV2G$pair[which(pairV2G$Celltype == rownames(pv)[i])]))
  }
}

pp = melt(as.matrix(pv))
pp = pp[-which(pp$value==0),]
pp = pp[order(pp$value),]
pv=NULL
for(i in 1:nrow(pp)){
  tmp = data.frame(Celltype=pp$Var1[i],Tissue=pp$Var2[i],
                   proxy_rsid=intersect(pair$proxy_rsid[which(pair$Tissue==pp$Var2[i])],
                                        pairV2G$proxy_rsid[which(pairV2G$Celltype==paste(pp$Var1[i]))]))
  pv = rbind(pv,tmp)
}

pvv=reshape2::dcast(Celltype~Tissue,data=pv,value.var = "proxy_rsid",fun.aggregate = function(x) paste(sort(x),collapse = ","))
rownames(pvv)=pvv$Celltype
pv=reshape2::melt(pvv,id.vars = 'Celltype')
colnames(pv)[2]="Tissue"
colnames(pv)[3]="proxy_rsid"
pv = pv[-which(pv$proxy_rsid==""),]
pv = merge(pv,allstat[,c('Celltype','Name')],by="Celltype")
pv$val = 0
for(i in 1:nrow(pv)){
  pv$val[i] = which(unique(pv$proxy_rsid)==pv$proxy_rsid[i])
}
tmp2=reshape2::dcast(Celltype~Tissue,value.var = "val",data=pv)
rownames(tmp2)=tmp2$Celltype
tmp2=as.matrix(t(tmp2[,-1]))
tmp2[is.na(tmp2)]=0
d <- dist(tmp2, method = "euclidean")
dend <- as.dendrogram(hclust(d, method = "complete"))
dend <- dendextend::seriate_dendrogram(dend, d)
tis_ord=rownames(tmp2)[order.dendrogram(dend)]
d <- dist(t(tmp2), method = "euclidean")
dend <- as.dendrogram(hclust(d, method = "complete"))
dend <- dendextend::seriate_dendrogram(dend, d)
cell_ord=rev(colnames(tmp2)[order.dendrogram(dend)])

pv$Celltype = factor(pv$Celltype,levels = paste(cell_ord))
pv$Tissue = factor(pv$Tissue,levels = paste(tis_ord))
celltypes = allstat[which(allstat$Celltype %in% cell_ord),c('Celltype','Name')]
celltypes$Celltype=factor(celltypes$Celltype,levels = paste(cell_ord))
celltypes = celltypes[order(celltypes$Celltype),]
pv$Name = factor(pv$Name,levels = paste(celltypes$Name))

pdf("ADCY3_coloc2.pdf",height = 3,width = 5)
ggplot(pv)+geom_tile(aes(Name,Tissue,fill=proxy_rsid))+theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  theme(legend.position = "none")
dev.off()