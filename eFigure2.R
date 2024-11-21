# eFigure 2: Examples of 2 variants mapped to only one gene promoter

library(reshape2)
library(nVennR)

snp_in_cRE = readRDS("snp_in_cRE.rds")
snp_in_loop = readRDS("snp_in_openloop.rds")
snp_in_openloop = readRDS("snp_in_openloop.rds")
snp_in_ocr = readRDS("snp_in_ocr.rds")

allstat=readRDS("allstat.rds")
Celltype= list()
for(i in 1:nrow(allstat)){
  n = which(colnames(snp_in_cRE)==allstat$Celltype[i])
  Celltype[[i]] = list(inLoop=rownames(snp_in_loop)[which(snp_in_loop[,i]!=0)],
                       inOpenLoop=rownames(snp_in_openloop)[which(snp_in_openloop[,i]!=0)],
                       inOCR = rownames(snp_in_ocr)[which(snp_in_ocr[,i]!=0)],
                       inCRE = rownames(snp_in_cRE)[which(snp_in_cRE[,n]!=0)]) 
  names(Celltype)[[i]] = paste(allstat$Celltype[i])
}

# Panel A
i=which(allstat$Celltype=="hESC")
plotVenn(Celltype[[i]],labelRegions = F,opacity = 0.3,fontScale=2,
           setColors=c("#4472C4","#FF0000","#FFC000","#70AD47"),
           outFile=paste0("SNPs within ",allstat$Celltype[i] ," _ by nVennR.svg")
)

# Panel B 
i=which(allstat$Celltype=="ME_YRI")
plotVenn(Celltype[[i]],labelRegions = F,opacity = 0.3,fontScale=2,
           setColors=c("#4472C4","#FF0000","#FFC000","#70AD47"),
           outFile=paste0("SNPs within ",allstat$Celltype[i] ," _ by nVennR.svg")
)

# trackviews were generated using pyGenomeTracks with correspoding cell type ATAC-seq bigwig and chromatin loops BEDPE files
# TFAP2B locus was plot at chr6:50711468-50846616
# GPR1 locus was plot at chr2:206201166-206214091