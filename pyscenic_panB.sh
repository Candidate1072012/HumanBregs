#!/bin/bash

# project name
#SBATCH -A immune-rep.prj

# job name
#SBATCH -J pyscenic_bcell_job

#SBATCH -o pyscenic_bcell_job-%j.out
#SBATCH -e pyscenic_bcell_job-%j.err

#SBATCH -p long
#SBATCH -c 10

# set working directory
#SBATCH -D /gpfs3/well/immune-rep/users/tma392/non_immune_PDAC/python


# Print some useful data about the job to help with debugging
echo "------------------------------------------------"
echo "SGE Job ID: $JOB_ID"
echo "SGE Job ID: $SGE_JOB_ID"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

# Finally, we can run our real computing job

module load Python/3.10.4-GCCcore-11.3.0
source SCENIC-skylake/bin/activate

# run pyscenic
dir=/gpfs3/well/immune-rep/users/tma392/non_immune_PDAC/python/cisTarget_databases

tfs=$dir/allTFs_hg38.txt
tbl=$dir/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
feather=$dir/hg38*rankings.feather

input_loom=./BCellCombined_sample.loom
ls $tfs  $feather  $tbl

#2.1 grn
pyscenic grn \
--num_workers 10 \
--output adj.sample.tsv \
--method grnboost2 \
BCellCombined_sample.loom \
$tfs

#2.2 cistarget
pyscenic ctx \
adj.sample.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom  \
--mode "dask_multiprocessing" \
--output reg.csv \
--num_workers 10  \
--mask_dropouts

#2.3 AUCell
pyscenic aucell \
$input_loom \
reg.csv \
--output BCellCombined_out_SCENIC.loom \
--num_workers 10
