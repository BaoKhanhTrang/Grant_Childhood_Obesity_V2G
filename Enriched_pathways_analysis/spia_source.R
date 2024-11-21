if(!require(SPIA)){
  BiocManager::install("SPIA")
  library(SPIA)
}
if(!require(KEGGREST)){
  BiocManager::install("KEGGREST")
  library(KEGGREST)
}
if(!require(EnrichmentBrowser)){
  BiocManager::install("EnrichmentBrowser")
  library(EnrichmentBrowser)
}
if(!require(limma)){
  BiocManager::install("limma")
  library(limma)
}
if(!require(splines)){
  install.packages("splines")
  library(splines)
}
if(!require(clusterProfiler)){
  install.packages("clusterProfiler")
  library(clusterProfiler)
}
if(!require(pathview)){
  BiocManager::install("pathview")
  library(pathview)
}

BC<-function(mdir=NULL,pathwaynames=NULL,org=org){
  library(igraph)
  library(KEGGgraph)
  betweennesslist<-NULL
  pathwayID<-NULL
  for(i in 1:length(pathwaynames)){
    mapkpathway<-try(parseKGML(paste0(mdir,"/",pathwaynames[[i]])),TRUE)
    mapkG<- KEGGpathway2Graph(mapkpathway, expandGenes=F)
    g<-igraph.from.graphNEL(mapkG)
    bet<-betweenness(g)
    nodlist<-NULL
    nodeslist<-NULL
    nod<-nodes(mapkpathway)
    for(j in 1:length(bet)){
      nodname<-names(bet[j])
      genename<-nod[[nodname]]@name
      
      for(jj in 1:length(genename)){
        betweenness<-rep(bet[nodname],length(genename))
      }
      nodlist<-t(rbind(genename,betweenness))
      nodeslist<-rbind(nodeslist,nodlist)
    }
    betness<-nodeslist[,2]
    betness<-1+as.numeric(betness)
    names(betness)<-nodeslist[,1]
    name<-names(betness)
    name<-strsplit(as.character(name),paste0(org,":"))
    name<-do.call(rbind,name)
    names(betness)<-name[,2]
    betweennesslist<-c(betweennesslist,list(betness))
    pathwayID<-c(pathwayID,mapkpathway@pathwayInfo@name)
  }	
  names(betweennesslist)<-pathwayID
  return(betweennesslist)
}

SP<-function(mdir=NULL,pathwaynames=NULL,org=org){
  library(igraph)
  library(KEGGgraph)
  nodeslist<-NULL
  pathwayID<-NULL
  for(i in 1:length(pathwaynames)){
    mapkpathway<-try(parseKGML(paste0(mdir,"/",pathwaynames[[i]])),TRUE)
    mapkG<- KEGGpathway2Graph(mapkpathway, expandGenes=T )
    nodes<-nodes(mapkG)
    nodeslist<-c(nodeslist,nodes)
  }
  specific_number<-table(unlist(nodeslist))
  wfg<-specific_number
  wfg<-as.data.frame(wfg)
  nodenames<-wfg$Var1
  nodenames<-strsplit(as.character(nodenames),paste0(org,":"))
  nodenames<-do.call(rbind,nodenames)
  wfg<-wfg[,2]
  names(wfg)<-nodenames[,2]
  return(wfg)
}

IF<-function(mdir=NULL,pathwaynames=NULL,DE=NULL,org=org){
  library(igraph)
  library(KEGGgraph)
  all<-paste0(org,':',DE)
  edge<-NULL
  nodeslist<-NULL
  for(i in 1:length(pathwaynames)){
    mapkpathway<-try(parseKGML(paste0(mdir,"/",pathwaynames[[i]])),TRUE)
    mapkG<-KEGGpathway2Graph(mapkpathway,expandGenes=TRUE)
    node<-nodes(mapkG)
    edL<-edgeData(mapkG)
    circsp<-strsplit(as.character(names(edL)),"\\|")
    geneName<-do.call(rbind,circsp)
    edge<-rbind(edge,geneName)
    nodeslist<-c(nodeslist,node)
  }
  nodeslist<-nodeslist[!duplicated(nodeslist)]
  e<-unique.matrix(edge)
  g<-graph_from_edgelist(e)
  mapk<-igraph.to.graphNEL(g)
  nodes<-nodes(mapk)
  neighborhoods<-ego(g,1,nodes,"out")
  inter<-function(X){
    nde<-length(intersect(names(X),all))
  }
  nDE<-lapply(neighborhoods,inter)
  nDE<-as.numeric(as.matrix(nDE))
  signodes<-setdiff(nodeslist,nodes)
  nDEsn<-rep(0,length(signodes))
  nDE<-c(nDE,nDEsn)
  names(nDE)<-c(nodes,signodes)
  nodenames<-strsplit(as.character(c(nodes,signodes)),paste0(org,":"))
  nodenames<-do.call(rbind,nodenames)
  names(nDE)<-nodenames[,2]
  return(nDE)
}
wi<-function(wfg=NULL,nDE=NULL){
  DEdegree<-nDE[names(wfg)]
  wig<-DEdegree/wfg
  #write.csv(wig,"./GWSPIA_w_8671.csv")
  wig<-1+((wig-min(wig))/(max(wig)-min(wig)))
  return(wig)
}
spiaDBS<-function(de=NULL,all=NULL,allt=NULL,betweennesslist=NULL,DEdegree=NULL,organism=organism,data.dir=NULL,
                  pathids=NULL,nB=2000,plots=FALSE,verbose=TRUE,beta=NULL,combine="fisher",plotname=plotname){
  library(SPIA)

  rel<-c("activation","compound","binding/association","expression","inhibition",
         "activation_phosphorylation","phosphorylation","inhibition_phosphorylation",
         "inhibition_dephosphorylation","dissociation","dephosphorylation",
         "activation_dephosphorylation","state change","activation_indirect effect",
         "inhibition_ubiquination","ubiquination", "expression_indirect effect",
         "inhibition_indirect effect","repression","dissociation_phosphorylation",
         "indirect effect_phosphorylation","activation_binding/association",
         "indirect effect","activation_compound","activation_ubiquination")



  if(is.null(beta)){
    beta=c(1,0,0,1,-1,1,0,-1,-1,0,0,1,0,1,-1,0,1,-1,-1,0,0,1,0,1,1)
    names(beta)<-rel
  }else{

    if(!all(names(beta) %in% rel) | length(names(beta))!=length(rel)){
      stop(paste("beta must be a numeric vector of length",length(rel), "with the following names:", "\n", paste(rel,collapse=",")))
    }
  }


   .myDataEnv <- new.env(parent=emptyenv()) # not exported

   datload<-paste(organism, "SPIA", sep = "")

   if(is.null(data.dir)){
     if(! paste(datload,".RData",sep="") %in% dir(system.file("extdata",package="SPIA"))){
       cat("The KEGG pathway data for your organism is not present in the extdata folder of the SPIA package!!!")
       cat("\n");
       cat("Please generate one first using makeSPIAdata and specify its location using data.dir argument or copy it in the extdata folder of the SPIA package!")
     } else{
       load(file=paste(system.file("extdata",package="SPIA"),paste("/",organism, "SPIA", sep = ""),".RData",sep=""), envir=.myDataEnv)
     }
   }
   if(!is.null(data.dir)){
     if (! paste(datload,".RData",sep="") %in% dir(data.dir)) {
       cat(paste(data.dir, " does not contin a file called ",paste(datload,".RData",sep="")))

     }else{
       load(file=paste(data.dir,paste(datload,".RData",sep=""),sep=""), envir=.myDataEnv)
     }
   }



   datpT=.myDataEnv[["path.info"]] # load SPIA data into datpT object

   if (!is.null(pathids)){
     if( all(pathids%in%names(datpT))){
       datpT=datpT[pathids]
     }else{
       stop( paste("pathids must be a subset of these pathway ids: ",paste(names(datpT),collapse=" "),sep=" "))
     }
   }

  datp<-list(); #normalized weighted directed adjacency matrix
  path.names<-NULL
  hasR<-NULL  # has reaction(s)

  for (jj in 1:length(datpT)){              #compute normalized weighted directed adjacency matrix for each pathway
    sizem<-dim(datpT[[jj]]$activation)[1]   #number of rows(number of nodes in this pathway)
    s<-0;
    con<-0;

    for(bb in 1:length(rel)){              # for each rel
      con=con+datpT[[jj]][[rel[bb]]]*abs(sign(beta[rel[bb]])) # sum all actions for each gene (regardless for direction)
      s=s+datpT[[jj]][[rel[bb]]]*beta[rel[bb]]                # sum actions for each gene
    }
    z=matrix(rep(apply(con,2,sum),dim(con)[1]),dim(con)[1],dim(con)[1],byrow=TRUE); # net total effect of each gene on the other genes in pathway, create new matrix
    z[z==0]<-1;   # gene with no action got a dummy value for division later

    datp[[jj]]<- s/z   #normalize the signed actions of eac gene by total actions of such gene across other gene within this pathway

    path.names<-c(path.names,datpT[[jj]]$title)
    hasR<-c(hasR,datpT[[jj]]$NumberOfReactions>=1)  # pathway has reaction
  }

  names(datp)<-names(datpT)
  names(path.names)<-names(datpT)

  #tor<-lapply(datp,function(d){sum(abs(d))})==0  | is.na(path.names) | hasR  # eliminate pathway that have no action/no name/ have reaction???
  #datp<-datp[!tor]
  #path.names<-path.names[!tor]



  IDsNotP<-names(de)[!names(de)%in%all]
  if(length(IDsNotP)/length(de)>0.01){stop("More than 1% of your de genes have IDs are not present in the reference array!. Are you sure you use the right reference array?")}
  if(!length(IDsNotP)==0){cat("The following IDs are missing from all vector...:\n");
    cat(paste(IDsNotP,collapse=","));
    cat("\nThey were added to your universe...");
    all<-c(all,IDsNotP)}

  if(length(intersect(names(de),all))!=length(de)){stop("de must be a vector of log2 fold changes. The names of de should be included in the refference array!")}

  ph<-pb<-pcomb<-sobs<-pSize<-smPFS<-tA<-tAraw<-KEGGLINK<-NULL;
  set.seed(1)

  if(plots){
    pdf(paste0(dir,"SPIA/",plotname,"_SPIAPerturbationPlots.pdf"))

  }
  #DEdegree<-nDE[names(wfg)]
  #DEdegree<-DEdegree/wfg
  for(i in 1:length(names(datp))){
    path<-names(datp)[i]
    M<-datp[[path]] # get each adjacency matrix
    path<-paste0("path:",organism,path)# fix the path name
    diag(M)<-diag(M)-1; # set self action to -1
    X<-de[rownames(M)] # set of DE gene present in this pathway, if not present value =0
    names(X)<-rownames(M)
    noMy<-sum(!is.na(X)) # number of DE gene in this pathway

    okg<-as.character(intersect(rownames(M),all)) #set of gene from pathway that present in the array
    nokg<-length(rownames(M))-length(okg)  #  number of gene in the pathway but not present in array

    wg<-DEdegree[okg] # get the normalized weight IF/SP of these genes
    tg<-allt[okg]     # get the t value of these genes

    np<-length(okg)
    sobs[i]<-sum(abs(tg)*wg,na.rm = TRUE)/np # "observed average t-Score" t-value multiply by weight, take absolute value, normalize for all considered genes
                                             # modified PADOG can be used to replace the original ORA method used in SPIA.

    ok<-rownames(M)%in%all # location within pathway that have genes of the pathways that are present in array
    pSize[i]<-length(rownames(M))

    if((noMy)>0&(abs(det(M))>1e-7)){ # if there are DE genes in this pathway and determinant of adjacency matrix >1e-7
      gnns<-paste(names(X)[!is.na(X)],collapse="+")
      KEGGLINK[i]<-paste("http://www.genome.jp/dbget-bin/show_pathway?",organism,names(datp)[i],"+",gnns,sep="")
      X[is.na(X)]<-0.0  # set gene in the pathway but have no value (NA) to 0 value

      b<-betweennesslist[[path]][as.character(names(X))] #get the betweenness of that pathway
      b[is.na(b)]<-1   # sub NA value as 1 (no effect)
      X<-b*X # multiply betweenness with differential expression value

      pfs<-solve(M,-X)  # solve linear algebraic equation -X = M*pfs
      smPFS[i]<-sum(pfs-X)
      tAraw[i]<-smPFS[i]


      if(plots){
        #if(interactive()){x11();par(mfrow=c(1,2))}
        par(mfrow=c(1,2))
        plot(X,pfs-X,main=paste("pathway ID=",names(datp)[i],sep=""),
             xlab="Log2 FC",ylab="Perturbation accumulation (Acc)",cex.main=0.8,cex.lab=1.2);abline(h=0,lwd=2,col="darkgrey");
        abline(v=0,lwd=2,col="darkgrey");#abline(0,1,lwd=2,col="darkgrey");
        points(X[abs(X)>0 & X==pfs],pfs[abs(X)>0 & X==pfs]-X[abs(X)>0 & X==pfs],col="blue",pch=19,cex=1.4)
        if(length(which(abs(X)>0 & X==pfs))>0){
          text(X[abs(X)>0 & X==pfs],pfs[abs(X)>0 & X==pfs]-X[abs(X)>0 & X==pfs], labels=names(X[abs(X)>0 & X==pfs]), cex=0.9, font=2,col="blue")
        }
        points(X[abs(X)>0 & X!=pfs],pfs[abs(X)>0 & X!=pfs]-X[abs(X)>0 & X!=pfs],col="red",pch=19,cex=1.4)
        if(length(which(abs(X)>0 & X!=pfs))>0){
          text(X[abs(X)>0 & X!=pfs],pfs[abs(X)>0 & X!=pfs]-X[abs(X)>0 & X!=pfs], labels=names(X[abs(X)>0 & X!=pfs]), cex=0.9, font=2,col="red")
        }
        points(X[abs(X)==0 & X==pfs],pfs[abs(X)==0 & X==pfs]-X[abs(X)==0 & X==pfs],col="black",pch=19,cex=1.4)
        if(length(which(abs(X)==0 & X==pfs))>0){
          text(X[abs(X)==0 & X==pfs],pfs[abs(X)==0 & X==pfs]-X[abs(X)==0 & X==pfs], labels=names(X[abs(X)==0 & X==pfs]), cex=0.9, font=2,col="black")
        }
        points(X[abs(X)==0 & X!=pfs],pfs[abs(X)==0 & X!=pfs]-X[abs(X)==0 & X!=pfs],col="green",pch=19,cex=1.4)
        if(length(which(abs(X)==0 & X!=pfs))>0){
          text(X[abs(X)==0 & X!=pfs],pfs[abs(X)==0 & X!=pfs]-X[abs(X)==0 & X!=pfs], labels=names(X[abs(X)==0 & X!=pfs]), cex=0.9, font=2,col="green")
        }
      }

      stmp<-NULL #random score
      for(ks in 1:nB){ # random select same number of genes within input geneset to calculate score
        tg1<-as.vector(sample(allt,np))
        stmp<-c(stmp,sum(abs(tg1)*wg)/np)
      }

      if(sobs[i]==0){ # if observed t-score =0
        if(all(stmp==0)){    # and all random sampling produce 0 t-score
          ph[i]<-NA         # probability for Functional class scoring (FCS) then is NA
        }else{
          ph[i]<-1   # otherwise probability =1
        }

      }
      if(sobs[i]>0){ # if observed t-score >0
        numb<-sum(stmp>=sobs[i],na.rm = TRUE)  # how many time random sampling produce higher t-score
        ph[i]<-numb/length(stmp)*2   # probability
        if(ph[i]<=0){ph[i]<-1/nB/100}
        if(ph[i]>1){ph[i]<-1}
      }

      pfstmp<-NULL # pertubation P

      for (k in 1:nB){
        x<-rep(0,length(X));names(x)<-rownames(M);
        x[ok][sample(1:sum(ok),noMy)]<-as.vector(sample(de,noMy))
        x<-b*x
        tt<-solve(M,-x)
        pfstmp<-c(pfstmp,sum(tt-x))#
      }

      mnn<-median(pfstmp)
      pfstmp<-pfstmp-mnn # normalize pfstmp
      ob<-smPFS[i]-mnn
      tA[i]<-ob

      if(ob>0){
        pb[i]<-sum(pfstmp>=ob)/length(pfstmp)*2
        if(pb[i]<=0){pb[i]<-1/nB/100}
        if(pb[i]>1){pb[i]<-1}
      }
      if(ob<0){
        pb[i]<-sum(pfstmp<=ob)/length(pfstmp)*2
        if(pb[i]<=0){pb[i]<-1/nB/100}
        if(pb[i]>1){pb[i]<-1}
      }
      if(ob==0){
        if(all(pfstmp==0)){    #there is nothing to learn from perturbations
          pb[i]<-NA
        }else{
          pb[i]<-1
        }

      }

      if(plots){
        bwidth = sd(pfstmp)/4
        if (bwidth > 0) {
          plot(density(pfstmp,bw=bwidth),cex.lab=1.2,col="black",lwd=2,main=paste("pathway ID=",names(datp)[i],"  P PERT=",round(pb[i],5),sep=""),
               xlim=c(min(c(tA[i]-0.5,pfstmp)),max(c(tA[i]+0.5,pfstmp))),cex.main=0.8,xlab="Total Perturbation Accumulation (TA)")
        } else {
          pfsTab = table(pfstmp)
          plot(as.numeric(names(pfsTab)), as.numeric(pfsTab), cex.lab=1.2,col="black",main=paste("pathway ID=",names(datp)[i],"  P PERT=",round(pb[i],5),sep=""),
               xlim=c(min(c(tA[i]-0.5,pfstmp)),max(c(tA[i]+0.5,pfstmp))),cex.main=0.8,xlab="Total Perturbation Accumulation (TA)", ylab="frequency")
        }
        abline(v=0,col="grey",lwd=2)
        abline(v=tA[i],col="red",lwd=3)

      }
      pcomb[i]<-combfunc(pb[i],ph[i],combine)

    }else{
      pb[i]<-ph[i]<-smPFS[i]<-pcomb[i]<-tAraw[i]<-tA[i]<-KEGGLINK[i]<-NA}

    if(verbose){
      cat("\n");
      cat(paste("Done pathway ",i," : ",substr(path.names[names(datp)[i]],1,30),"..",sep=""))
    }

  }#end for each pathway

  if(plots){
    par(mfrow=c(1,1))
    dev.off()
  }

  pcombFDR=p.adjust(pcomb,"fdr");phFdr=p.adjust(ph,"fdr")
  pcombfwer=p.adjust(pcomb,"bonferroni")
  Name=path.names[names(datp)]
  Status=ifelse(tA>0,"Activated","Inhibited")
  res<-data.frame(Name,ID=names(datp),pSize,tScore=sobs,pgsa=ph,tA,pPERT=pb,pG=pcomb,pGFdr=pcombFDR,pGFWER=pcombfwer,Status,KEGGLINK,stringsAsFactors=FALSE)

  res<-res[!is.na(res$pgsa),]
  res<-res[order(res$pG),]

  rownames(res)<-NULL;
  return(res)
}

## download KEGG xml + list the names ####
for(i in keys(KEGGPATHID2NAME)){
  download.kegg(pathway.id = paste(i),species = "hsa",kegg.dir = "Kegg/",file.type = "xml")
}
p = list.files(paste0(dir,"Kegg/"),pattern = ".xml",full.names = FALSE)
makeSPIAdata(kgml.path="~/Documents/analysis/V2G/Kegg",organism="hsa",out.path="~/Documents/analysis/V2G/Kegg")

## comput BC and SP #####
betweennesslist = BC(mdir = "Kegg",pathwaynames = p,org = "hsa")
wfg = SP(mdir = "Kegg",pathwaynames = p,org = "hsa")
