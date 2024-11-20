library(nVennR)
library(ggplot2)
library(UpSetR)

allv2g=readRDS("allv2g.rds")
allstat=readRDS("allstat.rds")

# Panel A ############################################################################################################
proxy = readRDS(file="proxy.rds")
snp_in_loop = readRDS(file="snp_in_loop.rds")
snp_in_openloop = readRDS(file="snp_in_openloop.rds")
snp_in_ocr = readRDS(file="snp_in_ocr.rds")

SNP = list(total = proxy$proxy_rsid,
            inLoops=rownames(snp_in_loop)[which(rowSums(snp_in_loop)!=0)],
           openLoop = rownames(snp_in_openloop)[which(rowSums(snp_in_openloop)!=0)],
           inOCR = rownames(snp_in_ocr)[which(rowSums(snp_in_ocr)!=0)],
           inCRE = unique(allv2g$proxy_rsid))

myvenn <- plotVenn(SNP,labelRegions = F,opacity = 0.3,fontScale=2,
                   setColors=c("white","#4472C4","#FF0000","#FFC000","#70AD47"),
                   outFile="CB_EGG.topLD_SNPs within all cell types _ by nVennR.svg")


# Panel B ############################################################################################################
loci=NULL
tmp = dcast(unique(allv2g[,c('Celltype','locus')]),locus~Celltype,fun.aggregate = length)
tmp$sum = rowSums(tmp[,-1])
tmp = tmp[,c('locus','sum')]
tmp$type = "# Mapped celltypes"
loci = rbind(loci,tmp)
tmp=dcast(unique(allv2g[,c('gene_name','locus')]),locus~gene_name,fun.aggregate = length)
tmp$sum=rowSums(tmp[,-1])
tmp = tmp[,c('locus','sum')]
tmp$type="# Target genes"
loci = rbind(loci,tmp)
tmp = dcast(unique(allv2g[,c('proxy_rsid','locus')]),locus~proxy_rsid,fun.aggregate = length)
tmp$sum = rowSums(tmp[,-1])
tmp = tmp[,c('locus','sum')]
tmp$type = "# variants"
loci = rbind(loci,tmp)
tmp=tmp[order(tmp$sum,decreasing = F),]
loci$locus = factor(loci$locus,levels = paste(tmp$locus))
loci$type = factor(loci$type,levels = c("# variants","# Mapped celltypes","# Target genes"))
pdf("CB_EGG.topLD_Locus stat.pdf",width = 6,height = 2.5)
ggplot(loci) + geom_bar(aes(locus,sum),stat = 'identity')+theme_bw()+
  geom_text(aes(locus,sum+2.5,label=sum))+
  facet_grid(cols = vars(type),scales = "free")+
  ylab("")+xlab("locus")+coord_flip()
dev.off()



# Panel C ############################################################################################################
upse = NULL
for(i in 1:nrow(allstat)){
  id = unique(allv2g$proxy_rsid[which(allv2g$Name==allstat$Name[i])])
  upse[[i]] = id
}
names(upse) = paste(allstat$Name)

pdf("CB_EGG.topLD_SNPs within cREs _ intersection between cell types.pdf",height = 7.5,width = 8.5)
upset(fromList(upse), order.by = "degree",nsets = 60,nintersects = NA,mb.ratio = c(0.3, 0.7),set_size.show=T,
      mainbar.y.label = "#SNPs Intersections", sets.x.label = "SNPs per Cell type",keep.order = F,
      queries = list(list(query = intersects,params = list(colnames(snp_in_cRE)[which(snp_in_cRE[names(which(rowSums(snp_in_cRE,na.rm = T)==39)),]!=0)]), color = "red", active = T),
                     list(query = intersects,params = list(colnames(snp_in_cRE)[which(snp_in_cRE['rs7132908',]!=0)]), color = "#1f77b4", active = T),
                     list(query = intersects,params = list(colnames(snp_in_cRE)[which(snp_in_cRE['rs6548240',]!=0)]), color = "#1f77b4", active = T),
                     list(query = intersects,params = list(colnames(snp_in_cRE)[which(snp_in_cRE[names(which(rowSums(snp_in_cRE,na.rm = T)==23)),]!=0)]), color = "#1f77b4", active = T),
                     list(query = intersects,params = list(colnames(snp_in_cRE)[which(snp_in_cRE[names(which(rowSums(snp_in_cRE,na.rm = T)==17)),]!=0)]), color = "#1f77b4", active = T),
                     list(query = intersects,params = list(colnames(snp_in_cRE)[which(snp_in_cRE[names(which(rowSums(snp_in_cRE,na.rm = T)==7)),]!=0)]), color = "#1f77b4", active = T),
                     list(query = intersects,params = list(colnames(snp_in_cRE)[which(snp_in_cRE[names(which(rowSums(snp_in_cRE,na.rm = T)==5)),]!=0)]), color = "#1f77b4", active = T),
                     list(query = intersects,params = list(colnames(snp_in_cRE)[which(snp_in_cRE['rs585944',]!=0)]), color = "#1f77b4", active = T),
                     list(query = intersects,params = list(colnames(snp_in_cRE)[which(snp_in_cRE['rs10865551',]!=0)]), color = "#1f77b4", active = T),
                     list(query = intersects,params = list(colnames(snp_in_cRE)[which(snp_in_cRE['rs36048912',]!=0)]), color = "#1f77b4", active = T),
                     list(query = intersects,params = list("NTERA2"), color = "grey76", active = T),
                     list(query = intersects,params = list("Myotubes"), color = "grey76", active = T),
                     list(query = intersects,params = list("EndoC_BH1"), color = "grey76", active = T),
                     list(query = intersects,params = list("NK"), color = "grey76", active = T),
                     list(query = intersects,params = list("Adipocytes"), color = "grey76", active = T),
                     list(query = intersects,params = list("Treg_act"), color = "grey76", active = T),
                     list(query = intersects,params = list("HNs_hESC"), color = "grey76", active = T),
                     list(query = intersects,params = list("pDC"), color = "grey76", active = T),
                     list(query = intersects,params = list("Osteo_hMSC"), color = "grey76", active = T),
                     list(query = intersects,params = list("PreAdipocytes"), color = "grey76", active = T),
                     list(query = intersects,params = list("Melano_YRI"), color = "grey76", active = T),
                     list(query = intersects,params = list("hESC"), color = "grey76", active = T),
                     list(query = intersects,params = list("Adipo_SGBS"), color = "grey76", active = T),
                     list(query = intersects,params = list("Enteroid"), color = "grey76", active = T),
                     list(query = intersects,params = list("Hepatocytes"), color = "grey76", active = T),
                     list(query = intersects,params = list("hFOB_Diff"), color = "grey76", active = T)
                     
      ))
dev.off()

# a few details of figure 2 were manually added using Power Point