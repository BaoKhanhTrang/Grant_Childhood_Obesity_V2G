############################# TRAIT Configuration #############################
#ID for user's trait of interest. (Can be any string)
trait="Childhood_obesity"
#path to the input GWAS summary stats files
traitFilePath="ldsc/sumstats/ChObes/CHILDHOOD_OBESITY.TRANS_ANCESTRAL.RESULTS.hg38"
#column IDs from trait file
trait_A1col="EA"
trait_A2col="OA"
trait_SNPcol="RSID"
trait_CHRcol="CHR"
trait_BPcol="POS"
trait_Pcol="EUR_P"
trait_Ncol="EUR_N"
trait_MAFcol="EUR_FRQ"

traitType="cc" # case-control
traitProp=".65"  # proportion of samples that are cases in a case control GWAS = cases / case + controls

#locus information for running coloc. 
chrom = CHROMOSOME
colocStart = STARTBP
colocStop = STOPBP
#reference genome build: "hg19" or "hg38"
build = "hg38"
lead_SNP = "SNPNUMBER"
#"eqtl" or "sqtl"
qtlType = "eqtl"
#colc window size in bps, this is only required when no leadSNPsFilePath is set. If you are providing a leadSNPsFile leave the empty string field for window
window="500000"

#config file with paths to the qtl data
setup_config_R = $colocdir/setup_config.R

