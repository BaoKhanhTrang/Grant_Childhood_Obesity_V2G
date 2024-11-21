# Figure 1: Partitioned Linkage Disequilibrium Score Regression analysis for open chromatin regions of all cell types.

library(ggplot2)
library(ggsci)
library(egg)
library(scales)
options(scipen = 100)

# celltypes library
allstat=readRDS("allstat.rds")

# read in results  ####
all=read.table("ldsc/out_atac/summary_CB_EGG_hg38_EUR.munged.atac",header=T)
prom=read.table("ldsc/out_Prom/summary_CB_EGG_hg38_EUR.munged.Prom",header=T)
nar=read.table("ldsc/out_narrow/summary_CB_EGG_hg38_EUR.munged.narrow",header=T)
aft=read.table("ldsc/out_expand/summary_CB_EGG_hg38_EUR.munged.expand",header=T)

all=merge(all,allstat[,c('Celltype','system','Name')],by='Celltype')
all$material="bulk"
all$material[grep("klaus_pancreas_cells",all$Celltype)]="single cell"
prom=merge(prom,allstat[,c('Celltype','system','Platform','Name')],by='Celltype')
prom$material="bulk"
prom$material[grep("klaus_pancreas_cells",prom$Celltype)]="single cell"
nar = merge(nar,allstat[,c('Celltype','system','Platform','Name')],by='Celltype')
nar$material="bulk"
nar$material[grep("klaus_pancreas_cells",nar$Celltype)]="single cell"
aft = merge(aft,allstat[,c('Celltype','system','Platform','Name')],by='Celltype')
aft$material="bulk"
aft$material[grep("klaus_pancreas_cells",aft$Celltype)]="single cell"

allstat$Name = factor(allstat$Name,levels = paste(allstat$Name))
all$Name = factor(all$Name,levels = paste(allstat$Name))
aft$Name = factor(aft$Name,levels = paste(allstat$Name))
nar$Name = factor(nar$Name,levels = paste(allstat$Name))
prom$Name = factor(prom$Name,levels = paste(allstat$Name))

g1=ggplot(allstat)+geom_bar(aes(atac_peaks_N,Name,fill=material),stat = "identity",width = .8)+
  geom_bar(aes(ocr_N,Name,fill=Platform),stat = "identity",width = .7)+
  facet_grid(rows=vars(system),scale = "free",space = "free",drop = TRUE)+
  scale_fill_d3(name="",breaks=c("bulk","single cell","hiC","capture-C"))+
  theme_bw()+
  scale_y_discrete(name="Cell types")+
  scale_x_continuous(expand = c(0,0),name = "Number of open chromatin regions",
                     breaks = seq(0,300000,by=60000),
                     labels = label_number(scale = 1e-3),
                     limits = c(0,305000))+
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),strip.text = element_blank(),
        legend.position = c(0.7,0.91),legend.background = element_blank())

maxp = ceiling(max(max(-log10(all$Enrichment_p)),
                   max(-log10(aft$Enrichment_p)),
                   max(-log10(nar$Enrichment_p)),
                   max(-log10(prom$Enrichment_p)) ))
minp = floor(min(c(min(-log10(all$Enrichment_p)),
                   min(-log10(aft$Enrichment_p)),
                   min(-log10(nar$Enrichment_p)),
                   min(-log10(prom$Enrichment_p)))))
maxbin = round(max(c(max(all$Prop._SNPs),
                     max(aft$Prop._SNPs),
                     max(nar$Prop._SNPs),
                     max(prom$Prop._SNPs))),2)
minbin = round(min(c(min(all$Prop._SNPs),
                     min(aft$Prop._SNPs),
                     min(nar$Prop._SNPs),
                     min(prom$Prop._SNPs))),2)
minenr = min(floor(all$Enrichment-all$Enrichment_std_error), # -12
             floor(aft$Enrichment-aft$Enrichment_std_error),
             floor(nar$Enrichment-nar$Enrichment_std_error),
             floor(prom$Enrichment-prom$Enrichment_std_error)) #17
maxenr = max(ceiling(all$Enrichment+all$Enrichment_std_error),
             ceiling(aft$Enrichment+aft$Enrichment_std_error),
             ceiling(nar$Enrichment+nar$Enrichment_std_error),
             ceiling(prom$Enrichment+prom$Enrichment_std_error))

g2=ggplot(all)+geom_vline(xintercept = 1, linetype='dashed',linewidth=.5)+
  geom_point(aes(Enrichment,factor(Name),size=`Prop._SNPs`,color=-log10(Enrichment_p)),alpha=0.7)+
  geom_errorbarh(aes(y=Name,xmax=Enrichment+Enrichment_std_error,xmin=Enrichment-Enrichment_std_error),color="#1F77B4FF",height=.5,linewidth=.4,position = "dodge")+
  geom_errorbarh(data=all[all$material=="single cell",],aes(y=Name,xmax=Enrichment+Enrichment_std_error,xmin=Enrichment-Enrichment_std_error),color="#D62728FF",height=.5,linewidth=.4,position = "dodge")+
  theme_bw()+
  geom_point(data=all[which(all$Enrichment_p < 0.05),],aes(Enrichment,Name),shape=8,color="white",size=2)+
  scale_size_binned(name = "Proportion of SNP",breaks = seq(0,maxbin,by=0.01),
                    labels = function(x) as.character(round(x,2)),limits = c(0,maxbin))+
  scale_x_continuous(name = "Total OCRs by ATAC-seq",sec.axis = dup_axis(),
                     breaks = c(1,round(seq(minenr,maxenr,length.out=6))))+
  coord_cartesian(xlim=c(minenr,maxenr))+
  scale_color_gradient(name="-log10 p-value",high = "#B2182B",low = "blue",limits = c(minp,maxp))+
  facet_grid(rows = vars(system),scales = "free",space = "free")+
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x.bottom = element_blank(),axis.text.x.top = element_blank(),axis.ticks.x.top = element_blank(),
        panel.grid.minor = element_blank(),legend.position = "none",
        strip.background = element_blank(),strip.text = element_blank())+
  guides(size = guide_bins(show.limits = TRUE))

g2b=ggplot(prom)+geom_vline(xintercept = 1, linetype='dashed',size=.5)+
  geom_point(aes(Enrichment,factor(Name),size=`Prop._SNPs`,color=-log10(Enrichment_p)),alpha=0.7)+
  geom_errorbarh(aes(y=Name,xmax=Enrichment+Enrichment_std_error,xmin=Enrichment-Enrichment_std_error),color="#1F77B4FF",height=.5,size=.4,position = "dodge")+
  geom_errorbarh(data=prom[prom$material=="single cell",],aes(y=Name,xmax=Enrichment+Enrichment_std_error,xmin=Enrichment-Enrichment_std_error),color="#D62728FF",height=.5,linewidth=.4,position = "dodge")+
  theme_bw()+
  geom_point(data=prom[which(prom$Enrichment_p < 0.05),],aes(Enrichment,Name),shape=8,color="white",size=2)+
  scale_size_binned(name = "Proportion of SNP",breaks = seq(0,maxbin,by=0.01),
                    labels = function(x) as.character(round(x,2)),limits = c(0,maxbin))+
  scale_x_continuous(name = "OCRs at promoters",sec.axis = dup_axis(),
                     breaks = c(1,round(seq(minenr,maxenr,length.out=6))))+
  coord_cartesian(xlim=c(minenr,maxenr))+
  scale_color_gradient(name="-log10 p-value",high = "#B2182B",low = "blue",limits = c(minp,maxp))+
  facet_grid(rows = vars(system),scales = "free",space = "free")+
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x.bottom = element_blank(),axis.text.x.top = element_blank(),axis.ticks.x.top = element_blank(),
        panel.grid.minor = element_blank(),legend.position = "none",
        strip.background = element_blank(),strip.text = element_blank())+
  guides(size = guide_bins(show.limits = TRUE))


g3=ggplot(nar)+geom_vline(xintercept = 1, linetype='dashed', linewidth=.5)+
  geom_point(aes(Enrichment,factor(Name),size=`Prop._SNPs`,color=-log10(Enrichment_p)),alpha=0.7)+
  geom_errorbarh(aes(y=Name,xmax=Enrichment+Enrichment_std_error,xmin=Enrichment-Enrichment_std_error),color="#FF7F0EFF",height=.5,size=.4,position = "dodge")+
  geom_errorbarh(data=nar[nar$Platform=="hiC",],aes(y=Name,xmax=Enrichment+Enrichment_std_error,xmin=Enrichment-Enrichment_std_error),color="#2CA02CFF",height=.5,linewidth=.4,position = "dodge")+
  theme_bw()+
  labs(title="Fold enrichment of heritability")+
  geom_point(data=nar[which(nar$Enrichment_p < 0.05),],aes(Enrichment,Name),shape=8,color="white",size=1)+
  scale_size_binned(name = "Proportion of SNP",breaks = seq(0,maxbin,by=0.01),
                    labels = function(x) as.character(round(x,2)),limits = c(0,maxbin))+
  scale_x_continuous(name = "putative cREs",sec.axis = dup_axis(),
                     breaks = c(1,round(seq(minenr,maxenr,length.out=6))))+
  coord_cartesian(xlim=c(minenr,maxenr))+
  scale_color_gradient(name="-log10 p-value",high = "#B2182B",low = "blue",limits = c(minp,maxp))+
  facet_grid(rows = vars(system),scales = "free",space = "free")+
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),plot.caption = element_text(hjust = 0.5,vjust = 0),
        axis.title.x.bottom = element_blank(),axis.text.x.top = element_blank(),axis.ticks.x.top = element_blank(),
        panel.grid.minor = element_blank(),legend.position = "none",
        plot.title = element_text(hjust = 0.5,margin=margin(b=-500)),
        strip.background = element_blank(),strip.text = element_blank())+
  guides(size = guide_bins(show.limits = TRUE))


g4=ggplot(aft)+geom_vline(xintercept = 1, linetype='dashed',size=.5)+
  geom_point(aes(Enrichment,factor(Name),size=`Prop._SNPs`,color=-log10(Enrichment_p)),alpha=0.7)+
  geom_errorbarh(aes(y=Name,xmax=Enrichment+Enrichment_std_error,xmin=Enrichment-Enrichment_std_error),color= "#FF7F0EFF",height=.5,size=.4,position = "dodge")+
  geom_errorbarh(data=aft[aft$Platform=="hiC",],aes(y=Name,xmax=Enrichment+Enrichment_std_error,xmin=Enrichment-Enrichment_std_error),color= "#2CA02CFF",height=.5,linewidth=.4,position = "dodge")+
  theme_bw()+
  geom_point(data=aft[which(aft$Enrichment_p < 0.05),],aes(Enrichment,Name),shape=8,color="white",size=2)+
  scale_size_binned(name = "Proportion of SNP",breaks = seq(0,maxbin,by=0.01),
                    labels = function(x) as.character(round(x,2)),limits = c(0,maxbin))+
  scale_x_continuous(name = "cREs ±500 bases",sec.axis = dup_axis(),
                     breaks = c(1,round(seq(minenr,maxenr,length.out=6))))+
  coord_cartesian(xlim=c(minenr,maxenr))+
  scale_color_gradient(name="-log10 p-value",high = "#B2182B",low = "blue",limits=c(minp,maxp))+
  facet_grid(rows = vars(system),scales = "free",space = "free")+
  theme(axis.title.y = element_blank(),axis.text.y=element_blank(),
        axis.title.x.bottom = element_blank(),axis.text.x.top = element_blank(),axis.ticks.x.top = element_blank(),
        panel.grid.minor = element_blank())+
  guides(size = guide_bins(show.limits = TRUE))


pdf("LDSC_dotplot_newCelltype_mungedSumstat.pdf",width = 11,height = 7)
ggarrange(g1,g2,g2b,g3,g4,ncol = 5,widths = c(4.2,3,3,3,3))
dev.off()