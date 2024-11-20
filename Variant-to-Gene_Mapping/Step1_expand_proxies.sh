dir=v2g/SNPset

# 19 sentinel signals that achieved genome-wide significance in the trans-ancestral meta-analysis study
#Bradfield, Jonathan P et al. “A trans-ancestral meta-analysis of genome-wide association studies reveals loci associated with childhood obesity.” Human molecular genetics vol. 28,19 (2019): 3327-3338. doi:10.1093/hmg/ddz161
# 19 rsids (hg19) saved in CB20snps_hg19 file
sed "s/rs//g" CB20snps_hg19 file  > tmp
python LiftRsNumber.py tmp > tmp2
paste -d'\t' <(cut -f2 tmp2 | sed "s/^ /rs/g" ) > leads.txt

# use topLD to expand proxies
topld_api -inFile leads.txt -outputInfo topLD_info.txt -outputLD topLD_output.txt -pop "EUR"  -thres 0.8


