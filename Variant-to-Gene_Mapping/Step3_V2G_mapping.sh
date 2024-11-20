DIR=v2g
proxy=$DIR/SNPset/CB_EGG_proxy_topLD.bed
prefix="CB_EGG.topLD"
outdir=$DIR/${prefix}_v2g

mkdir -p $outdir
mkdir $DIR/${prefix}_v2loop
mkdir $DIR/${prefix}_v2openloops
mkdir $DIR/${prefix}_v2ocr

for cell in $(cat $DIR/celltypes.txt ); do
    ANNO_DIR=$DIR/celltypes/$cell
    OUT_PREFIX=${prefix}_$cell
    cd $DIR/celltypes/$cell
    Rscript v2g_proxybed_cisanno.R $proxy $ANNO_DIR $OUT_PREFIX $outdir

    cut -f1-3 $DIR/celltypes/$cell/$cell.cis.anno.bedpe | sort -u > $DIR/celltypes/$cell/$cell.tmp
    cut -f4-6 $DIR/celltypes/$cell/$cell.cis.anno.bedpe | sort -u >> $DIR/celltypes/$cell/$cell.tmp
    intersectBed -wa -u -a $DIR/celltypes/$cell/$cell.tmp -b $DIR/celltypes/$cell/atac.bed | sort -u > $DIR/celltypes/$cell/$cell.loopbed
    intersectBed -wa -u -a $proxy -b $DIR/celltypes/$cell/$cell.loopbed > $DIR/${prefix}_v2openloops/$cell.snps
    
    intersectBed -wa -u -a $proxy -b $DIR/celltypes/$cell/atac.bed > $DIR/${prefix}_v2ocr/$cell.snps;

    pairToBed -a $DIR/celltypes/$cell/${cell}_cis.bedpe -b $proxy | cut -f11 | sort -u > $DIR/${prefix}_v2loop/$cell.snps

done

