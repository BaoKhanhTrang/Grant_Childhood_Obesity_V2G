#!/bin/bash

cdir=ldsc
scriptdir=~/tools
cd $cdir/sumstats/ChObes

# Download sumstat
wget http://egg-consortium.org/Childhood_Obesity_2019/CHILDHOOD_OBESITY.TRANS_ANCESTRAL.RESULTS.txt.gz

# Lift from hg19 to hg38
zcat CHILDHOOD_OBESITY.TRANS_ANCESTRAL.RESULTS.txt.gz | sed '1d' | awk '{print "chr" $2 "\t" $3-1 "\t" $3 "\t" $1}' > CBtransancestral.hg19
$scriptdir/liftOver CBtransancestral.hg19 $scriptdir/hg19ToHg38.over.chain.gz CBtransancestral.hg38 CBtransancestral.unlift
for i in {1..22} "X" ; do
        chr="chr"$i"_"
        for pat in $(grep "^$chr" CBtransancestral.hg38 | cut -f1 | sort -u) ; do
                sed  -i "s/$pat/chr$i/g"  CBtransancestral.hg38
        done
done
for chr in $(cut -f1 tmp.lift | sort -u) ; do
        i=$(echo $chr | sed "s/chr//g")
        grep -w "^$chr" CBtransancestral.hg38 | cut -f3 > pos
        grep -w -f pos <(zcat $scriptdir/annotationFiles/GWAS_reference/dbSNP/hg38/bed/bed_chr_$i.bed.gz | cut -f1,3,4 )  >> tmp2
done
join -1 2 -2 4 -o 1.1,1.3,1.4,2.3 <(awk '{$2=$1"_"$3 ; print $0}' CBtransancestral.hg38 | sort -k2,2) <(awk '{$4=$1"_"$2; print $0}' tmp2 | sort -k4,4) | sort -u > genoPos2rsID_hg38
echo "CHR START POS RSID alleles EA OA AFR_N AFR_FRQ AFR_BETA AFR_SE AFR_P ASN_N ASN_FRQ ASN_BETA ASN_SE ASN_P EUR_N EUR_FRQ EUR_BETA EUR_SE EUR_P AMR_N AMR_FRQ AMR_BETA AMR_SE AMR_P BAYES" > 
join -1 6 -2 1 -t $'\t' \
        -o 1.1,1.2,1.3,1.4,1.5,1.6,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14,2.15,2.16,2.17,2.18,2.19,2.20,2.21,2.22,2.23,2.24,2.25,2.26 \
        <(sort -t $'\t' -k 6b,6 genoPos2rsID_hg38) \
        <(sed '1d' CHILDHOOD_OBESITY.TRANS_ANCESTRAL.RESULTS.txt | sort -t $'\t' -k 1b,1) > \
        CHILDHOOD_OBESITY.TRANS_ANCESTRAL.RESULTS.hg38

# munge the sumstat 
conda activate ldsc
python $scriptdir/ldsc/munge_sumstats.py \
        --sumstats CHILDHOOD_OBESITY.TRANS_ANCESTRAL.RESULTS.hg38 \
        --a1 EA \
        --a2 OA \
        --snp RSID \
        --signed-sumstats EUR_BETA,0 \
        --p EUR_P \
        --N-col EUR_N \
        --frq EUR_FRQ \
        --out CB_EGG_hg38_EUR.munged
