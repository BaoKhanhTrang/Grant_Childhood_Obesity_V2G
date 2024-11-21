# eFigure 11: motifbreakR vs ATAC-seq footprint analysis

library(motifbreakR)
library(ggstatsplot)

mb=readRDS("CB_EGG.topLD_94V2G_motifBreak.rds")
proxy = readRDS("proxy.rds")
allv2g = readRDS("allv2g.rds")

# Panel A was generated with pyGenomeTrack chr12:49842317-49966000
# with 25 cell types grouped into 4 system classes

plotMB(results = mb, rsid = "rs7132908", effect = "strong", altAllele = "C")


# Panel B
olap=NULL
for(i in 1:nrow(allstat)){
    #read in overlapped proxies and TF-footprint from ATAC-seq RGT-HINT analysis
  filename = paste0("atac/",allstat$Celltype[i],"/",allstat$Celltype[i],".mpbs.txt") 
  if(file.size(filename)!=0){
    mo = read.delim(filename,header = F,sep="\t")
    mo = mo[,-c(7,8,10,11,14)]
    colnames(mo)=c("chr","mo_start","mo_end","motif","score","strand",
                   "nucleosome_pos","variant_pos","id")
    mo$Name = allstat$Name[i]
    olap= rbind(olap,mo)
  }
}
tmp=unlist(strsplit(olap$id,split = "\\|"))
olap$index_rsid = sub("index_rsid:","",paste(tmp[seq(1,length(tmp),by=2)]))
olap$proxy_rsid = sub("proxy_rsid:","",paste(tmp[seq(2,length(tmp),by=2)]))
olap$TF=""
tmp = unlist(strsplit(olap$motif[grep("HUMAN",olap$motif)],split="_HU"))
olap$TF[grep("HUMAN",olap$motif)] = paste(tmp[seq(1,by=2,length.out=length(grep("HUMAN",olap$motif)))])
tmp = unlist(strsplit(sub("(var.3)","",sub("(var.2)","",olap$motif[-grep("HUMAN",olap$motif)])) ,split="\\."))
olap$TF[-grep("HUMAN",olap$motif)] = paste(tmp[seq(3,length(tmp),by=3)])


proxy$motifBreakR="non-annotated"
proxy$motifBreakR[which(proxy$proxy_rsid %in% unique(mb$SNP_id))]="motifBreakR-annotated"
proxy$MPBS="non-overlapped"
proxy$MPBS[which(proxy$proxy_rsid %in% unique(olap$proxy_rsid))]="TF footprint-overlapped"
rownames(proxy)=proxy$proxy_rsid
df=proxy[,c('motifBreakR','MPBS')]
test <- fisher.test(table(df),)

pdf("CB_EGG.topLD_motifbreak_vs_HINT_olap_fisher.pdf",height = 4,width = 5)
ggbarstats(
  df, MPBS,motifBreakR,
  results.subtitle = FALSE,
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))
  )
)
dev.off()

# Panel C
olap = olap[which(olap$proxy_rsid %in% allv2g$proxy_rsid),]
rm(filename,tmp,i)

fp=unique(olap[,c('proxy_rsid','TF')])
tmp=mb[which(mb$SNP_id %in% fp$proxy_rsid)]
tmp = as.data.frame(tmp,row.names = seq(1,length(tmp)))
tmp = unique(tmp[,c('SNP_id','geneSymbol','effect')])
colnames(tmp)[1:2]=colnames(fp)
simi = merge(fp,tmp)

olap_cell=unique(olap[,c('proxy_rsid','TF','Name')])
olap_cell=merge(olap_cell,simi,by=c('proxy_rsid','TF'))
olap_cell=unique(olap_cell[,1:3])
olap_cell = merge(olap_cell,unique(allv2g[,c('proxy_rsid','Name')]),by=c('proxy_rsid','Name'))
rs=unique(simi$proxy_rsid)
locus=unique(allv2g[which(allv2g$proxy_rsid %in% rs),c('proxy_rsid','locus')])
olap_cell = merge(olap_cell,locus,by="proxy_rsid")


out=as.data.frame(mb,row.names = seq(1,length(mb))) 
out$SNP_id = paste(unlist(out$SNP_id))
outt = unique(out[which(out$effect=="strong"),c('SNP_id','geneSymbol')])
outtt = as.data.frame(table(paste(outt$SNP_id)))
colnames(outtt)=c("proxy_rsid","N.pred.TF")
for(i in 1:nrow(outtt)){
  MPBS = unique(olap$TF[which(olap$proxy_rsid == outtt$proxy_rsid[i])])
  predTF = unique(outt$geneSymbol[which(outt$SNP_id==outtt$proxy_rsid[i])])
  outtt$olap[i] = length(intersect(MPBS,predTF))
  outtt$N.MPBS[i] = length(setdiff(MPBS,predTF))
  outtt$N.pred.TF[i] = length(setdiff(predTF,MPBS))
  cell = unique(allv2g$Name[which(allv2g$proxy_rsid==outtt$proxy_rsid[i])])
  cellFP = unique(olap$Name[which(olap$proxy_rsid == outtt$proxy_rsid[i])])
  outtt$N.cell[i] = length(setdiff(cell,cellFP))
  outtt$N.cellFP[i] = length(setdiff(cellFP,cell))
  outtt$cell.olap[i] = length(intersect(cellFP,cell))
}
pred = reshape2::melt(outtt[,1:4],id.vars = "proxy_rsid")
colnames(pred)[2:3]=c("TF","N.TF")
tmp = reshape2::melt(outtt[,c(1,5:7)],id.vars = "proxy_rsid")
colnames(tmp)[2:3]=c("cell","N.cell")
pred=cbind(pred,tmp[,2:3])
pred$TF = factor(pred$TF,levels = c("N.MPBS","olap","N.pred.TF"))
pred$proxy_rsid = factor(pred$proxy_rsid,levels = proxy$proxy_rsid)
pred$cell = factor(pred$cell,levels = c("N.cellFP","cell.olap","N.cell"))
pred = merge(pred,unique(allv2g[,c('proxy_rsid','locus')]),by="proxy_rsid")
pred$locus = factor(pred$locus,levels = unique(proxy$locus))
g1=ggplot(pred)+geom_bar(aes(proxy_rsid,N.TF,group=factor(TF),fill=TF),stat = 'identity')+
  theme_bw()+
  scale_fill_manual(values = c("#F0A0FF" ,"#0075DC", "#94FFB5" ))+
  scale_y_continuous(expand = c(0,0),name="Number of TFs")+
  facet_grid(cols = vars(locus),drop = TRUE,scales = "free",space = "free")+
  theme(strip.text.x = element_text(angle = 90),panel.spacing = unit(0.05, "lines"),
          axis.title.x = element_blank(), legend.position = "none",
        axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))

pdf("CB_EGG.topLD_motifbreak_vs_HINT.pdf",height = 3.5,width = 11)
print(g1)
dev.off()

# eFigure 12: Genome view of the 7 variants and the corresponding TF motifs  

exam=pred[which(pred$TF=="olap" & pred$N.TF!=0),]
exam = exam[order(exam$proxy_rsid),]

for(i in 1:nrow(exam)){
  oneRSID = mb[mb$SNP_id==exam$proxy_rsid[i]]
  oneRSID = oneRSID[oneRSID$effect=="strong"]
  oneRSID = oneRSID[oneRSID$geneSymbol %in% intersect(outt$geneSymbol[outt$SNP_id==exam$proxy_rsid[i]],
                                                               fp$TF[which(fp$proxy_rsid==exam$proxy_rsid[i])])]
  rs=oneRSID[which(oneRSID$geneSymbol == unique(oneRSID$geneSymbol)[1])]
  wi = sapply(rs$motifPos, `[`, 2) - sapply(rs$motifPos, `[`, 1)
  rs = rs[which(wi == max(wi))]
  rs = rs[which(rs$scoreAlt == max(rs$scoreAlt))]
  tmp=rs
  for(e in 2:length(unique(oneRSID$geneSymbol))){
    rs = oneRSID[which(oneRSID$geneSymbol == unique(oneRSID$geneSymbol)[e])]
    wi = sapply(rs$motifPos, `[`, 2) - sapply(rs$motifPos, `[`, 1)
    rs = rs[which(wi == max(wi))]
    rs = rs[which(rs$scoreAlt == max(rs$scoreAlt))]
    tmp = c(tmp,rs)
  }
  oneRSID = tmp
  oneRSID = calculatePvalue(oneRSID,granularity = 1e-6)
  wid = max(unlist(oneRSID$motifPos)) - min(unlist(oneRSID$motifPos))
  pdf(paste0("CB_EGG.topLD_motifbreak_vs_HINT.",exam$proxy_rsid[i],".pdf"),
      width = 0.45*wid+0.5,height = 3+length(oneRSID))
  plotMB(oneRSID,rsid = exam$proxy_rsid[i],altAllele = unique(oneRSID$ALT))
  dev.off()
}
