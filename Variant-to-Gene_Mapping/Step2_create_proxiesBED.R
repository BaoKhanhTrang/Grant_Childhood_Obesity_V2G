setwd("v2g/SNPset")
# use LDlinkR to expand proxies

library(LDlinkR)
snps = readLines("leads.txt")

LDproxy_batch(snps,token = Sys.getenv("LDLINK_TOKEN"),genome_build = "grch38_high_coverage",pop="EUR")
add=NULL
for(i in 1:length(snps)){
        file=list.files(".",pattern=snps[i])
        if(length(file)==1){
                tmp=read.table(paste0("./",file),header=T)
                tmp=tmp[which(round(tmp$R2,1) >= 0.8),]
                pos=unlist(strsplit(tmp$Coord,split=":"))
                tmp$chr=paste(pos[seq(1,by=2,length(pos))])
                tmp$end=as.numeric(paste(pos[seq(2,by=2,length(pos))]))
                tmp$start = tmp$end -1
                tmp$index_rsid = snps[i]
                add = rbind(add,tmp)
        }
}

ld=read.delim("topLD_output.txt",header=T,sep="\t")

info=read.delim("topLD_info.txt",header=T,sep="\t")
info$Position=as.numeric(paste(info$Position))
if("None" %in% info$CHROM){
    info=info[-which(info$CHROM=="None"),]
}

proxy=data.frame(chr=paste0("chr",LD$CHROM),start=LD$POS2.gb38. -1, end =LD$POS2.gb38.,index=paste0("|index_rsid:",LD$rsID1,"|proxy_rsid:",LD$rsID2))
proxy = unique(proxy)
tmp=data.frame(chr=paste0("chr",info$CHROM),start=info$Position - 1,end=info$Position,index=paste0("|index_rsid:",info$rsID,"|proxy_rsid:",info$rsID))
tmp=unique(tmp)
proxy=unique(rbind(proxy,tmp))
tmp=data.frame(chr=add$chr,start=add$start,end=add$end,index=paste0("index_rsid:",add$index_rsid,"|proxy_rsid:",add$RS_Number))
proxy=unique(rbind(proxy,tmp))

proxy=proxy[order(as.numeric(gsub("chr","",proxy$chr)),as.numeric(proxy$start)),]
write.table(proxy,file="CB_EGG_proxy_topLD.bed",col.names=F,row.names=F,sep="\t",quote=F)
