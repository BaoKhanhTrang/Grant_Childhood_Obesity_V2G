# Use ATAC-seq BAM files and peaks BED files for motif analysis using RGT-HINT tool

proxies=v2g/SNPset/CB_EGG_proxy_topLD.bed


for cell in $(cat celltypes.txt ); do
    dir="atac/"$cell

    rgt-hint footprinting --atac-seq --paired-end --organism=hg38 \
        --output-location=$dir \
        --output-prefix=$cell.footprint \
        $dir/$cell.bam \
        $dir/atac.bed

    rgt-motifanalysis matching --organism=hg38 \
        --filter "species:sapiens;database:hocomoco,jaspar_vertebrates" \
        --output-location=$dir/motif_analysis/MotifMatching \
        --input-files $dir/$cell.footprint.bed \
        --rand-proportion 10

    rgt-motifanalysis enrichment --organism=hg38 \
        --filter "species:sapiens;database:hocomoco,jaspar_vertebrates" \
        --matching-location=$dir/motif_analysis/MotifMatching \
        --output-location=$dir/motif_analysis/enrichment \
        $dir/motif_analysis/MotifMatching/random_regions.bed \
        $dir/$cell.footprint.bed
    
    # Overlap the TF footprint with all the proxies
    path=$dir/motif_analysis/enrichment/$cell.footprint/mpbs_ev.bed
    intersectBed -wo -a $path -b $proxies > $$dir/$cell.mpbs.txt
done


