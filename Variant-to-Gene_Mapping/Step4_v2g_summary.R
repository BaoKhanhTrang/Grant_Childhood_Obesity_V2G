allstat=readRDS("allstat.rds")

# create V2g table
allv2g=NULL
for(i in 1:nrow(allstat)){
  tmpp =read.delim(paste0("v2g/CB_EGG.topLD_v2g/CB_EGG.topLD_",allstat$Celltype[i],".prom_v2g.txt"),header=T)
  if(nrow(tmpp)!=0){tmpp$OCRtype="Prom"
  tmpp$source = allstat$Celltype[i]
  }
  tmp =read.delim(paste0("v2g/CB_EGG.topLD_v2g/CB_EGG.topLD_",allstat$Celltype[i],".cis_v2g.txt"),header=T)
  if(nrow(tmp)!=0){
    tmp$source = allstat$Celltype[i]
    tmp$OCRtype="PIR"
  }
  colnames(tmp)[9]=colnames(tmpp)[9]="Celltype"
  allv2g = rbind(allv2g,tmp,tmpp)
}
rm(tmp,tmpp)
prom=read.table("v2g/hg38_genecode.v30_promoter_reference.txt",header=TRUE)
prom=prom[order(as.numeric(paste(gsub("chr","",prom$pro_chr))),as.numeric(prom$pro_start),decreasing = F),]
allv2g$gene_id=factor(allv2g$gene_id,level=unique(prom$gene_id)) 
allv2g=allv2g[order(allv2g$gene_id),]
allv2g$gene_name =factor(allv2g$gene_name,levels = paste(unique(allv2g$gene_name)))
tmp=unlist(strsplit(allv2g$variant_id,split = "\\|"))
allv2g$index_rsid = sub("index_rsid:","",paste(tmp[seq(1,length(tmp),by=2)]))
allv2g$proxy_rsid = sub("proxy_rsid:","",paste(tmp[seq(2,length(tmp),by=2)]))

tmp = reshape2::dcast(data = unique(allv2g[,c('gene_name','Celltype')]),formula = gene_name~Celltype,fun.aggregate = length)
p1=as.data.frame(sapply(tmp,function(x) length(which(x!=0)))) # number of implicated genes
p1$Celltype=rownames(p1)
p1=p1[-1,]
colnames(p1)[1]="Number of\nimplicated genes"
allstat=merge(allstat,p1,by="Celltype",all.x = TRUE)

tmp = dcast(data = unique(allv2g[,c('proxy_rsid','Celltype')]),formula = proxy_rsid~Celltype,fun.aggregate = length)
p1=as.data.frame(sapply(tmp,function(x) length(which(x!=0)))) # number of SNPs
colnames(p1)="Number of\nacting SNPs"
p1$Celltype=rownames(p1)
p1=p1[-1,]
allstat=merge(allstat,p1,by="Celltype",all.x = TRUE)

saveRDS(allv2g,file='allv2g.rds')

# statistics of all proxies
proxy=read.table(("v2g/SNPset/CB_EGG_proxy_topLD.bed"),header=F)
colnames(proxy)[1:3]=c("chr","start","variant_pos")
tmp=unlist(strsplit(proxy$V4,split = "\\|"))
proxy$index_rsid = sub("index_rsid:","",paste(tmp[seq(1,length(tmp),by=2)]))
proxy$proxy_rsid = sub("proxy_rsid:","",paste(tmp[seq(2,length(tmp),by=2)]))
loci=read.table("v2g/SNPset/19lead_SNPs.txt",header = TRUE,sep='\t')
proxy = merge(proxy,loci,by="index_rsid")
proxy$loci = paste(proxy$locus,proxy$index_rsid,sep = ": ")
proxy = proxy[order(proxy$chr,proxy$variant_pos,decreasing = F),]
proxy$loci = factor(proxy$loci,levels = unique(proxy$loci))
locus=proxy[which(proxy$index_rsid==proxy$proxy_rsid),]
locus=locus[order(as.numeric(paste(gsub("chr","",locus$chr))),as.numeric(locus$variant_pos)),]
proxy$locus = factor(proxy$locus,levels = unique(locus$locus))

snp_in_loop = snp_in_openloop = snp_in_ocr = data.frame(matrix(0,nrow = nrow(proxy),ncol = nrow(allstat)))
rownames(snp_in_loop) = rownames(snp_in_openloop) = rownames(snp_in_ocr) = proxy$proxy_rsid
colnames(snp_in_loop) = colnames(snp_in_openloop) = colnames(snp_in_ocr) = allstat$Celltype
for(i in 1:nrow(allstat)){
  inp = read.table(paste0("v2g/CB_EGG.topLD_v2loop/",allstat$Celltype[i],".snps"))
  id=unique(sub("proxy_rsid:","",paste(unlist(strsplit(inp$V1,split = "\\|"))[seq(2,nrow(inp)*2,by=2)])))
  snp_in_loop[which(rownames(snp_in_loop) %in% id),which(colnames(snp_in_loop)==allstat$Celltype[i])] = 1

  inp = read.table(paste0("v2g/CB_EGG.topLD_v2openloops/",allstat$Celltype[i],".snps"))
  id=unique(sub("proxy_rsid:","",paste(unlist(strsplit(inp$V4,split = "\\|"))[seq(2,nrow(inp)*2,by=2)])))
  snp_in_openloop[which(rownames(snp_in_openloop) %in% id),which(colnames(snp_in_openloop)==allstat$Celltype[i])] = 1
  
  inp = read.table(paste0("v2g/CB_EGG.topLD_v2ocr/",allstat$Celltype[i],".snps"))
  id=unique(sub("proxy_rsid:","",paste(unlist(strsplit(inp$V4,split = "\\|"))[seq(2,nrow(inp)*2,by=2)])))
  snp_in_ocr[which(rownames(snp_in_ocr) %in% id),which(colnames(snp_in_ocr)==allstat$Celltype[i])] = 1

}
allstat$num_snp_in_loops = colSums(snp_in_loop)
proxy$num_cell_in_loop = rowSums(snp_in_loop)
allstat$num_snp_in_openloop = colSums(snp_in_openloop)
proxy$num_cell_in_openloop = rowSums(snp_in_openloop)
allstat$num_snp_in_ocr = colSums(snp_in_ocr)
proxy$num_cell_in_ocr = rowSums(snp_in_ocr)

snp_in_cRE = dcast(proxy_rsid~Celltype,data=unique(allv2g[,c('proxy_rsid','Celltype')]),fun.aggregate = length,fill = 0)
rownames(snp_in_cRE)=snp_in_cRE$proxy_rsid
snp_in_cRE = snp_in_cRE[proxy$proxy_rsid,]
rownames(snp_in_cRE)=proxy$proxy_rsid
snp_in_cRE = snp_in_cRE[,-1]
proxy$num_cell_in_cRE = rowSums(snp_in_cRE,na.rm = T)
snp_in_cRE[,setdiff(allstat$Celltype,colnames(snp_in_cRE))]=0
snp_in_cRE = snp_in_cRE[,paste(allstat$Celltype)]
allstat$num_snp_in_cRE = colSums(snp_in_cRE,na.rm = T)
snp_in_cRE[is.na(snp_in_cRE)]=0

saveRDS(allstat,file="allstat.rds")
saveRDS(proxy,file="proxy.rds")
saveRDS(snp_in_loop,file="snp_in_loop.rds")
saveRDS(snp_in_openloop,file="snp_in_openloop.rds")
saveRDS(snp_in_ocr,file="snp_in_ocr.rds")
saveRDS(snp_in_cRE,file="snp_in_cRE.rds")

# create loci for colocalization analysis
leads = unique(allv2g[,c('index_rsid','chr')])
colnames(leads)[1]="proxy_rsid"
leads = merge(leads,proxy[,c('proxy_rsid','variant_pos')],by="proxy_rsid")
colnames(leads)[1]="rsID"
colnames(leads)[3]="pos"
leads$start=0
leads$stop=0
for(i in 1:nrow(leads)){
  start = leads$pos[i]-50000
  end = leads$pos[i]+50000
  blk_start=min(allv2g$variant_pos[which(allv2g$index_rsid == leads$rsID[i])])
  blk_end=max(allv2g$variant_pos[which(allv2g$index_rsid == leads$rsID[i])])
  mid = round(mean(start,end,blk_start,blk_end))
  leads$start[i] = mid - 250000
  leads$stop[i] = mid + 250000
}
leads$chr = sub("chr","",leads$chr)
leads = leads[order(as.numeric(leads$chr),leads$pos,decreasing = F),]
write.table(leads,file="CB_EGG.13sentinels_loci.txt",col.names = F,row.names = F,sep="\t",quote = F)
