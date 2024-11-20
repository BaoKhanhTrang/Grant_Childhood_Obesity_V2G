############################# eQTL Configuration #############################
# Download GRCh38 dbSNP BED files from here: https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/BED/
# Run dbsnp_hash_table_maker.py to create hash tables from the GRCh38 dbSNP BED files in the directory containing the files.
# https://github.com/bvoightlab/ColocQuiaL/blob/main/dbsnp_hash_table_maker.py
cd GWAS_reference/dbSNP/hg38/bed/
python dbsnp_hash_table_maker.py

#provide path to GRCh38 variant hash table directory that contain
hash_table_dir = "GWAS_reference/dbSNP/hg38/bed/"

#provide path to eqtl tissue table
eQTL_tissue_table = "GTEx_v8_Tissue_Summary_with_filenames.csv"

#provide path to significant eqtl data tabix directory NOTE: don't include last slash
eQTL_sig_qtl_tabix_dir = "GTEx_Analysis_v8_eQTL"
#column number of column in significant eqtl files that contains geneID
eQTL_sig_geneID_col = 7

#provide path to the all eqtl data tabix directory
eQTL_all_qtl_tabix_dir = "GTEx_Analysis_v8_eQTL_all_associations/"
#provide the header to the eqtl data
eQTL_all_header = c("chrom_b38", "chromStart_b38", "chromEnd_b38", "eGeneID", "A1_eqtl", "A2_eqtl", "build", "tss_distance", "ma_samples", "ma_count", "maf", "pvalue_eQTL", "slope", "slope_se")
#name of column in header representing geneID
eQTL_all_geneID = "eGeneID"
#name of column in header representing chromosome
eQTL_all_chrom = "chrom_b38"
#name of column in header representing end coordinate 
eQTL_all_chromEnd = "chromEnd_b38"
#name of column in header representing p-value
eQTL_all_pvalue = "pvalue_eQTL"
