dir="pathway/"
setwd(dir)

metaData=readRDS("metaData.rds")
allstat=readRDS("allstat.rds")

cells=intersect(metaData$Celltype,allstat$Celltype)

# for pathfindR
combine=readRDS(paste0("pathfindR/",cells[1],"/",cells[1],"_pathfindR.rds"))
colnames(combine)[3:7]=paste0(colnames(combine)[3:7],"_",cells[1])

for(i in 2:length(cells)){
  file=paste0("pathfindR/",cells[i],"/",cells[i],"_pathfindR.rds")
  if(file.exists(file)){
    df=readRDS(file)
    colnames(df)[3:7]=paste0(colnames(df)[3:7],"_",cells[i])
    combine=merge(combine,df,by=c("ID","Term_Description"),all=TRUE)
  }
}
saveRSD(combine,file="pathfindR/combine.rds")

# for mSPIA
spia=NULL

for(i in 1:length(cells)){
  if(file.exists(pattern = paste0("SPIA/",cells[i],".rds"))){
    tmp=readRDS(paste0("SPIA/",cells[i],".rds"))
    if(nrow(tmp)>0){
      tmp$Celltype = cells[i]
      tmp$trait=trait_ord[t]
      spia = rbind(spia,tmp)
    }
  }
}
saveRDS(spia,file="SPIA/combine.rds")