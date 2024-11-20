leadSNPsFilePath="v2g/CB_EGG.13sentinels_loci.txt"
colocdir=coloc
source $colocdir/trait_config_template.R
source $colocdir/setup_config.R

############################# run coloc for each locus #############################
cat $leadSNPsFilePath | while read line ; do
        echo $line
        #parse the SNP
        SNP=`echo $line | cut -f 1 -d " "`
        echo $SNP
        #parse chr start stop for coloc
        CHR=`echo $line | cut -f 2 -d " "`
        echo $CHR
        Start=`echo $line | cut -f 4 -d " "`
        echo $Start
        Stop=`echo $line | cut -f 5 -d " "`
        echo $Stop

        #make the analysis directory
        mkdir -p $colocdir/$trait/$SNP

        #cd into the dir
        cd $colocdir/$trait/$SNP

        #copy the colocquial.R file into the dir
        cp $colocdir/$trait/coloc.R ./

        #add an edited version of the QTL_config.R file to the dir (edit gtex file path, chr, start, and stop with sed)
        sed "s/SNPNUMBER/$SNP/" $colocdir/trait_config_template.R | sed "s/CHROMOSOME/$CHR/" | sed "s/STARTBP/$Start/" | sed "s/STOPBP/$Stop/" > $colocdir/$trait/$SNP/QTL_config.R

        # Run
        Rscript coloc.R

        #cd back into the main directory to go to the next SNP 
        cd $colocdir

done

echo "all lead SNP jobs have been submitted"

#set variables for date and time
today=`date +%F`
now=`date +%T`

#write header to the summary results file
printf "SNP\tGene\tGeneID-Tissue\tTrait\tPPID\tPP\n" > $colocdir/$trait"_"$qtlType"_coloc_results_all_summary.txt"
#go through each lead SNP directory
cat $leadSNPsFilePath | while read line ; do
    lead=`echo $line | cut -f 1 -d " "`
    dir=$colocdir/$trait/$lead
    cd  $dir

    #for each GTEx coloc output in the directory
    for colocOut in ./*coloc_results_summary.txt
    do
        #echo $colocOut
        colocfilename=`basename $colocOut`
        #echo $colocfilename

        #save a "trait" string for formatting
        trait_str="_"$trait"_"

        #grab the information we need from the locus coloc results file
        sed "s/^/$colocfilename /" $colocOut | sed "s/^/$lead /" | sed "s,$trait_str, $trait," | sed "s/coloc_results_summary.txt//" | sed "s/_ENSG/ ENSG/"| tr " " "\t" | tail -n+3 >> $colocdir/$trait"_"$qtlType"_coloc_results_all_summary.txt"
    done

    cd $colocdir

done

Rscript $colocdir"/summarize_qtl_coloc_PP3_PP4_results.R" $colocdir/$trait"_"$qtlType"_coloc_results_all_summary.txt"
rm $colocdir/$trait"_"$qtlType"_coloc_results_all_summary.txt"
echo "Your summary QTL file is ready!" 
