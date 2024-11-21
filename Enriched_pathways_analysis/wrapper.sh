dir="pathway"

# run DE analysis consider each cell type versus all other
for cell in $(cat celltypes.txt ); do
    Rscript DE_analysis.R $cell
done

# Use DE results as input for pathfindR
# run for each cell type
for cell in $(cat celltypes.txt ); do
    Rscript pathfindR_anaalysis.R $cell
done

# Use DE results as input for mSPIA
# run for each cell type
for cell in $(cat celltypes.txt ); do
    Rscript mSPIA_analysis.R $cell
done

#combine the results
Rscrript combine_results.R

# create WGCNA modules
Rscript WGCNA_analysis.R
