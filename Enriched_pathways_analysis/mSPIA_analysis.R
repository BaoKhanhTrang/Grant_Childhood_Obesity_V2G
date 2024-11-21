args = commandArgs(trailingOnly=TRUE)
cell=args[1]

library(clusterProfiler)
library(org.Hs.eg.db)
library(DESeq2)

dir="pathway/"

metaData=readRDS("metaData.rds")
allExp=readRDS(,"allExp.rds")
allstat=readRDS("allstat.rds")

eg = bitr(paste(unique(allExp$gene_name)), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
colnames(eg)[1]="gene_name"
ALL = paste(unique(eg$ENTREZID))
ALL = ALL[!is.na(ALL)]
ALL = ALL[!(ALL=="NA")]

file=list.files(path=paste0(dir,"DEanalysis/"),pattern=paste0(cell,"_overAll.apeglm.rds"))
a=readRDS(paste0(dir,"DEanalysis/",file))
# make sure the fold change is cell type/other, it not, reverse the log
check=strsplit(sub("Wald test p-value: Celltype ","",a@elementMetadata$description[4]),split=" vs ")[[1]][[1]]
if(check=="other"){
        a$tmp = 1/(2^(a$log2FoldChange))
        a$logFC = log2(a$tmp)
}else{
        colnames(a)[2]="logFC"
}
a=as.data.frame(a)
if(any(is.na(a$pvalue))){a=a[-which(is.na(a$pvalue)),]}
a$padj = p.adjust(a$pvalue,method = "BH")

a=merge(a,eg,by="gene_name")
allt = a$stat
names(allt) = paste(a$ENTREZID)
allt = allt[-which(names(allt)=="NA")]

# retrieve V2G of cell type as input
genes = unique(allv2g$gene_name[which(allv2g$Celltype==cell)])

# SPIA
source(paste0(dir,"spia_source.R"))

de = unique(a[which(a$gene_name %in% genes),c('ENTREZID','logFC')])
DE = de$logFC
names(DE) = as.vector(paste(de$ENTREZID))

nDE = IF(mdir=paste0(dir,"Kegg"),pathwaynames = p,DE = paste(de$ENTREZID),org = "hsa")
wig=wi(wfg = wfg,nDE = nDE)
res=spiaDBS(de=DE,all=ALL,allt = allt,betweennesslist = betweennesslist,DEdegree = wig,
    organism = "hsa",plots = TRUE,data.dir = paste0(dir,"Kegg/"),combine = "norminv",plotname = paste0("pathways_for_",cell)
)

saveRDS(res,file=paste0(dir,"SPIA/",cell,".rds")) 

