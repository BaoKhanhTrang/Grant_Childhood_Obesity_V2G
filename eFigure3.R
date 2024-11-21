# eFigure 3: The number of genes and cell types each proxy mapped into our cREs

library(reshape2)
library(ggplot2)
allv2g=readRDS("allv2g.rds")
proxy = readRDS("proxy.rds")

# Panel A
snp_gene = dcast(proxy_rsid~gene_name,data=unique(allv2g[,c('proxy_rsid','gene_name')]),fun.aggregate = length)
snp_gene$num_gene = rowSums(snp_gene[,-1])
snp_gene = snp_gene[,c('proxy_rsid','num_gene')]
tmp = dcast(proxy_rsid~Name,data=unique(allv2g[,c('proxy_rsid','Name')]),fun.aggregate = length)
tmp$num_cell = rowSums(tmp[,-1])
snp_gene = merge(snp_gene,tmp[,c('proxy_rsid','num_cell')],by="proxy_rsid")
snp_gene = snp_gene[order(snp_gene$num_gene,snp_gene$num_cell),]
snp_gene = merge(snp_gene,proxy[,c('proxy_rsid','locus')],by="proxy_rsid")
snp_gene=snp_gene[order(snp_gene$num_gene,snp_gene$num_cell),]

pdf("CB_EGG.topLD_SNPs across cells across genes_dotplot.pdf",height = 4,width = 5.5)
ggplot(snp_gene)+geom_point(aes(num_gene,num_cell,color=locus))+
  theme_bw()+
  theme(panel.border = element_blank(),strip.background = element_blank(),
        strip.text = element_blank())
dev.off()

# Panel B
snp_gene = dcast(proxy_rsid~gene_name,data=unique(allv2g[,c('proxy_rsid','gene_name')]),fun.aggregate = length)
snp_gene$num = rowSums(snp_gene[,-1])
snp_gene = snp_gene[,c('proxy_rsid','num')]
snp_gene$type="Number of genes"
tmp = dcast(proxy_rsid~Name,data=unique(allv2g[,c('proxy_rsid','Name')]),fun.aggregate = length)
tmp$num = rowSums(tmp[,-1])
tmp$type="Number of cell types"
snp_gene = rbind(snp_gene,tmp[,c('proxy_rsid','num','type')])
tmp = as.data.frame(table(proxy$locus[which(proxy$proxy_rsid %in% allv2g$proxy_rsid)]))
colnames(tmp)=c('locus','num_snp')

snp_gene = merge(snp_gene,proxy[,c('proxy_rsid','locus')],by="proxy_rsid")
snp_gene = snp_gene[order(snp_gene$num),]
snp_gene$proxy_rsid = factor(snp_gene$proxy_rsid,levels = paste(unique(snp_gene$proxy_rsid)))

pdf("CB_EGG.topLD_SNPs across cells across genes.pdf",height = 6,width = 12)
ggplot(snp_gene) + 
  geom_bar(aes(proxy_rsid,num,fill=locus),stat = 'identity',width = .8)+
  scale_y_continuous(breaks = c(1,unique(round(as.numeric(names(table(snp_gene$num)))/10)*10)),expand = c(0,0))+
  facet_grid(rows = vars(type),scales = "free")+
  theme_bw()+theme(axis.text.x=element_text(angle = 90,hjust = 1,vjust = 0.5),panel.grid.minor = element_blank(),
                   panel.border = element_blank(),strip.background = element_blank(),
                   legend.position = "none",axis.title = element_blank())
dev.off()