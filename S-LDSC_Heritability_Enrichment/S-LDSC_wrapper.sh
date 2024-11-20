#!/bin/bash

scriptdir=~/tools
cdir=ldsc
plink=$cdir/plink_files_hg38/1000G.EUR.hg38
frq=$cdir/plink_files_hg38/1000G.EUR.hg38
baseline=$cdir/baselineLD_v2.2_hg38/baselineLD
weight=$cdir/weights_hg38/weights.hm3_noMHC

# Run partitioned LD score regression #########################
gwas_sumstat=$cdir/sumstats/ChObes/CB_EGG_hg38_EUR.munged.sumstats.gz
gwas_base=($(basename $gwas_sumstat ".sumstats.gz"))

for base in $(cat celltypes.txt ); do
    dir="$base""_ldscore"    

    for cate in "narrow" "atac" "Prom" "expand" ; do
        annodir=$cdir/annotationfile_$cate
        outdir=$cdir/out_$cate

		python $scriptdir/ldsc/ldsc.py --h2 $gwas_sumstat \
            --ref-ld-chr $annodir/$dir/$base.,$baseline.  \
            --w-ld-chr $weight.  \
            --overlap-annot  \
            --frqfile-chr $frq.  \
            --out $outdir/${base}_$gwas_base
    done
done

# Summarize results
for cate in "narrow" "atac" "Prom" "expand" ; do
    outdir=$cdir/out_$cate
	echo "Celltype	Prop._SNPs	Prop._h2	Prop._h2_std_error	Enrichment	Enrichment_std_error	Enrichment_p" > $outdir/summary_$gwas_base.$cate
	ls $cdir/$outdir/*_${gwas_base}.results | while read file ; do
        base=$(basename $file "_${gwas_base}.results" )
        res=$(grep "L2_0" $cdir/$outdir/${base}_${gwas_base}.results | cut -f2- )
        echo $base $res | sed "s/ /\t/g" >> $outdir/summary_$gwas_base.$cate
	done
done
