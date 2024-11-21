library(motifbreakR)
library(BSgenome)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(XML)

allstat=readRDS("allstat.rds")
allv2g=readRDS("allv2g.rds")
proxy=readRDS("proxies.rds")

snp=read.delim("v2g/SNPset/topLD_output.txt",header=T)
ind = unique(snp[,c('CHROM','POS1.gb38.','rsID1','MAF1','REF1','ALT1')])
ind = ind[,c('CHROM','POS1.gb38.','rsID1','MAF1','REF1','ALT1','strand')]
pro = unique(snp[,c('CHROM','POS2.gb38.','rsID2','MAF2','REF2','ALT2','Sign')])
colnames(ind)=colnames(pro)=c("chr","pos","SNP_id","MAF","REF","ALT","strand")
snp=rbind(ind,pro)
rownames(snp)=snp$SNP_id
snp$chr = paste0("chr",snp$chr)
snp$start=snp$pos-1
snps.mb = makeGRangesFromDataFrame(snp,keep.extra.columns = TRUE,
                                   end.field = "pos",seqnames.field = "chr",
                                   starts.in.df.are.0based = TRUE,
                                   seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38))
attributes(snps.mb)$genome.package <- attributes(BSgenome.Hsapiens.UCSC.hg38)$pkgname

change.to.search.genome <- function(granges.object, search.genome) {
  sequence <- seqlevels(granges.object)
  ## sequence is in UCSC format and we want NCBI style
  newStyle <- mapSeqlevels(sequence,seqlevelsStyle(search.genome))
  newStyle <- newStyle[complete.cases(newStyle)] # removing NA cases.
  ## rename the seqlevels
  granges.object <- renameSeqlevels(granges.object,newStyle)
  seqlevels(granges.object) <- seqlevelsInUse(granges.object)
  seqinfo(granges.object) <- keepSeqlevels(seqinfo(search.genome),
                                           value = seqlevelsInUse(granges.object))
  return(granges.object)
}
snps.mb = change.to.search.genome(snps.mb,BSgenome.Hsapiens.UCSC.hg38)
snps.mb$REF=DNAStringSet(paste(snps.mb$REF))
snps.mb$ALT=DNAStringSet(paste(snps.mb$ALT))

data(hocomoco)
data(motifbreakR_motif)

mb <- motifbreakR(snpList = snps.mb, filterp = TRUE,
                   pwmList = MotifDb,
                   threshold = 0.005,
                   method = "ic",
                   bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                   BPPARAM = BiocParallel::SerialParam())

saveRDS(mb,file="CB_EGG.topLD_94V2G_motifBreak.rds")