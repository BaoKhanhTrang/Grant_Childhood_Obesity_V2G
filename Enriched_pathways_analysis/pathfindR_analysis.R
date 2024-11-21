args = commandArgs(trailingOnly=TRUE)
cell=args[1]

library(pathfindR)

dir="pathway/"
allv2g = readRDS("allv2g.rds")

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
a$ensembl_gene_id=rownames(a)

# remove NA p-values
if(any(is.na(a$pvalue))){a=a[-which(is.na(a$pvalue)),]}

# run pathfindR 
genes = unique(allv2g[which(allv2g$Celltype==cell ),c('ensembl_gene_id','gene_name')])
if(nrow(genes)>=2){
    dir.create(paste0(dir,"findpathR/",cell),recursive=TRUE)
    setwd(paste0(dir,"findpathR/",cell))

    # retrieve V2G of cell type as input
    input = a[which(a$gene_name %in% genes),c('gene_name','logFC','pvalue')]

    # save processed-input for visualization
    input_processed=NULL
    res <- tryCatch({inp1 = input_processing(input,p_val_threshold = 1,pin_name_path = "Biogrid")},
                    error = function(e) {return(NULL)},
                    warning = function(w) {return(NULL)})
    if(exists("inp1")){input_processed=rbind(input_processed,inp1)}

    res <- tryCatch({inp2 = input_processing(input,p_val_threshold = 1,pin_name_path = "KEGG")},
                    error = function(e) {return(NULL)},
                    warning = function(w) {return(NULL)})
    if(exists("inp2")){input_processed=rbind(input_processed,inp2)}

    res <- tryCatch({inp3 = input_processing(input,p_val_threshold = 1,pin_name_path = "STRING")},
                    error = function(e) {return(NULL)},
                    warning = function(w) {return(NULL)})
    if(exists("inp3")){input_processed=rbind(input_processed,inp3)}


     res <- tryCatch({inp4 = input_processing(input,p_val_threshold = 1,pin_name_path = "GeneMania")},
                    error = function(e) {return(NULL)},
                    warning = function(w) {return(NULL)})
    if(exists("inp4")){input_processed=rbind(input_processed,inp4)}


    res <- tryCatch({inp5 = input_processing(input,p_val_threshold = 1,pin_name_path = "IntAct")},
                    error = function(e) {return(NULL)},
                    warning = function(w) {return(NULL)})
    if(exists("inp5")){input_processed=rbind(input_processed,inp5)}
    saveRDS(input_processed,"input_processed.rds")

    # run with different PIP db

    combine <- run_pathfindR(input,min_gset_size = 1,p_val_threshold=1,output_dir="./Biogrid",pin_name_path="Biogrid",visualize_enriched_terms =F,plot_enrichment_chart=F)

    res2 <- run_pathfindR(input,min_gset_size = 1,p_val_threshold=1,output_dir="./KEGG",pin_name_path="KEGG",visualize_enriched_terms =F,plot_enrichment_chart=F)
    if(nrow(res2)!=0){
        if(nrow(combine)!=0){
            combine = combine_pathfindR_results(combine,res2,plot_common = FALSE)
        }else{
            combine = res2
        }
    }

    res3 <- run_pathfindR(input,min_gset_size = 1,p_val_threshold=1,output_dir="./STRING",pin_name_path="STRING",visualize_enriched_terms =F,plot_enrichment_chart=F)
    if(nrow(res3)!=0){
        if(nrow(combine)!=0){
            combine = combine_pathfindR_results(combine,res3,plot_common = FALSE)
        }else{
            combine = res3
        }
    }

    res4 <- run_pathfindR(input,min_gset_size = 1,p_val_threshold=1,output_dir="./GeneMania",pin_name_path="GeneMania",visualize_enriched_terms =F,plot_enrichment_chart=F)
    if(nrow(res4)!=0){
        if(nrow(combine)!=0){
            combine = combine_pathfindR_results(combine,res4,plot_common = FALSE)
        }else{
            combine = res4
        }
    }

    res5 <- run_pathfindR(input,min_gset_size = 1,p_val_threshold=1,output_dir="./IntAct",pin_name_path="IntAct",visualize_enriched_terms =F,plot_enrichment_chart=F)
    if(nrow(res5)!=0){
        if(nrow(combine)!=0){
            combine = combine_pathfindR_results(combine,res5,plot_common = FALSE)
        }else{
            combine = res5
        }
    }

    # combine the results from all PIP db

    if(nrow(combine)>0){
        final = combine[,1:2]
        final$Fold_Enrichment=0
        final$lowest_p=0
        final$highest_p=0
        final$Up_regulated=""
        final$Down_regulated=""
        for(r in 1:nrow(combine)){
            final$Fold_Enrichment[r] = max(unlist(combine[r,grep("Fold_Enrichment",colnames(combine))]),na.rm = T)
            final$lowest_p[r] = min(unlist(combine[r,grep("lowest_p",colnames(combine))]),na.rm = T)
            final$highest_p[r] = min(unlist(combine[r,grep("highest_p",colnames(combine))]),na.rm = T)
            gene = unique(unlist(strsplit(na.exclude(unlist(combine[r,grep("Up_regulated",colnames(combine))])),split=", ")))
            if(length(gene)!=0){
                final$Up_regulated[r] = paste(gene,collapse = ", ")
            }
            gene = unique(unlist(strsplit(na.exclude(unlist(combine[r,grep("Down_regulated",colnames(combine))])),split=", ")))
            if(length(gene)!=0){
                final$Down_regulated[r] = paste(gene,collapse = ", ")
            }
        }
        out=which(final$Up_regulated=="" & final$Down_regulated=="")
        if(length(out)>0){
            final = final[-out,]
        }

        # save result of each cell type
        saveRDS(final,file=paste0(cell,"_pathfindR.rds"))
    }
}