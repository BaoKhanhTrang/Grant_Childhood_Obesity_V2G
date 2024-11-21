# eFigure 1: Mapping statistics of the childhood obesity-associated loci

library(reshape2)
library(ggplot2)
library(patchwork)
library(nVennR)
library(UpSetR)
library(eulerr)
library(tidyverse)  
proxy=readRDS("proxy.rds")
allv2g=readRDS("allv2g.rds")

SNP = list(total = proxy$proxy_rsid,
            inLoops=rownames(snp_in_loop)[which(rowSums(snp_in_loop)!=0)],
           openLoop = rownames(snp_in_openloop)[which(rowSums(snp_in_openloop)!=0)],
           inOCR = rownames(snp_in_ocr)[which(rowSums(snp_in_ocr)!=0)],
           inCRE = unique(allv2g$proxy_rsid))

# A panel
i=1
snp = list(total = proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],
             inLoops = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inLoops),
             openLoop = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$openLoop),
             inOCR = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inOCR),
             inCRE = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inCRE))
euler_plot <- euler(snp)
g1 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)

i=i+1
snp = list(total = proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],
             inLoops = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inLoops),
             openLoop = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$openLoop),
             inOCR = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inOCR),
             inCRE = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inCRE))
euler_plot <- euler(snp)
g2 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)

i=i+1
snp = list(total = proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],
             inLoops = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inLoops),
             openLoop = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$openLoop),
             inOCR = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inOCR),
             inCRE = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inCRE))
euler_plot <- euler(snp)  g3 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)
g3 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)

i=i+1
snp = list(total = proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],
             inLoops = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inLoops),
             openLoop = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$openLoop),
             inOCR = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inOCR),
             inCRE = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inCRE))
euler_plot <- euler(snp)  g3 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)
g4 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)

i=i+1
snp = list(total = proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],
             inLoops = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inLoops),
             openLoop = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$openLoop),
             inOCR = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inOCR),
             inCRE = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inCRE))
euler_plot <- euler(snp)  g3 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)
g5 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)

i=i+1
snp = list(total = proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],
             inLoops = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inLoops),
             openLoop = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$openLoop),
             inOCR = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inOCR),
             inCRE = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inCRE))
euler_plot <- euler(snp)  g3 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)
g6 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)

i=i+1
snp = list(total = proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],
             inLoops = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inLoops),
             openLoop = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$openLoop),
             inOCR = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inOCR),
             inCRE = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inCRE))
euler_plot <- euler(snp)  g3 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)
g7 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)

i=i+1
snp = list(total = proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],
             inLoops = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inLoops),
             openLoop = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$openLoop),
             inOCR = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inOCR),
             inCRE = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inCRE))
euler_plot <- euler(snp)  g3 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)
g8 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)

i=i+1
snp = list(total = proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],
             inLoops = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inLoops),
             openLoop = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$openLoop),
             inOCR = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inOCR),
             inCRE = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inCRE))
euler_plot <- euler(snp)  g3 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)
g9 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)

i=i+1
snp = list(total = proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],
             inLoops = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inLoops),
             openLoop = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$openLoop),
             inOCR = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inOCR),
             inCRE = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inCRE))
euler_plot <- euler(snp)  g3 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)
g10 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)

i=i+1
snp = list(total = proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],
             inLoops = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inLoops),
             openLoop = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$openLoop),
             inOCR = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inOCR),
             inCRE = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inCRE))
euler_plot <- euler(snp)  g3 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)
g11 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)

i=i+1
snp = list(total = proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],
             inLoops = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inLoops),
             openLoop = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$openLoop),
             inOCR = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inOCR),
             inCRE = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inCRE))
euler_plot <- euler(snp)  g3 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)
g12 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)

i=i+1
snp = list(total = proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],
             inLoops = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inLoops),
             openLoop = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$openLoop),
             inOCR = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inOCR),
             inCRE = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inCRE))
euler_plot <- euler(snp)  g3 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)
g13 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)

i=i+1
snp = list(total = proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],
             inLoops = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inLoops),
             openLoop = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$openLoop),
             inOCR = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inOCR),
             inCRE = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inCRE))
euler_plot <- euler(snp)  g3 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)
g14 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)

i=i+1
snp = list(total = proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],
             inLoops = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inLoops),
             openLoop = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$openLoop),
             inOCR = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inOCR),
             inCRE = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inCRE))
euler_plot <- euler(snp)  g3 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)
g15 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)

i=i+1
snp = list(total = proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],
             inLoops = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inLoops),
             openLoop = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$openLoop),
             inOCR = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inOCR),
             inCRE = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inCRE))
euler_plot <- euler(snp)  g3 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)
g16 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)

i=i+1
snp = list(total = proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],
             inLoops = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inLoops),
             openLoop = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$openLoop),
             inOCR = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inOCR),
             inCRE = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inCRE))
euler_plot <- euler(snp)  g3 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)
g17 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)

i=i+1
snp = list(total = proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],
             inLoops = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inLoops),
             openLoop = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$openLoop),
             inOCR = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inOCR),
             inCRE = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inCRE))
euler_plot <- euler(snp)  g3 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)
g18 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)

i=i+1
snp = list(total = proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],
             inLoops = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inLoops),
             openLoop = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$openLoop),
             inOCR = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inOCR),
             inCRE = intersect(proxy$proxy_rsid[which(proxy$loci==loci$loci[i])],SNP$inCRE))
euler_plot <- euler(snp)  g3 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)
g19 = plot(euler_plot,main = paste(loci$loci[i]),quantities = TRUE,fill=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),alpha=0.8)

pdf("CB_EGG.topLD_SNPs within each locus - euler.pdf",height = 15,width = 15)
gridExtra::grid.arrange(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,
             g11,g12,g13,g14,g15,g16,g17,g18,g19, ncol = 5)
dev.off()

# Panel B 
upse = NULL
for(i in 1:nrow(allstat)){
  id = unique(allv2g$locus[which(allv2g$Name==allstat$Name[i])])
  upse[[i]] = id
}
names(upse)=paste(allstat$Name)
pdf("CB_EGG.topLD_locus within cREs _ intersection between cell types.pdf",height = 7,width = 4)
upset(fromList(upse), order.by = "degree",nsets = 58,nintersects = NA,mb.ratio = c(0.3, 0.7),set_size.show=T,
      mainbar.y.label = "#Loci Intersections", sets.x.label = "Loci per Cell type",keep.order = F)
dev.off()
