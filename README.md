# scRNAsnakemake


Snakefile: defines where the samples lie, how the files are named, etc.

snakemakerules.STARrules.singleEnd: holds the most general rules for this

config.json: specifies for each rule how much computational power it gets on the cluster

countMatrix_standalone.R: called by the rule that produces a single countmatrix for all cells


run pipeline with:
snakemake -p -j 80 --cluster-config /home/wertek/scripts/snakemake/config.json --cluster
'qsub -l walltime={cluster.time},mem={cluster.mem},nodes=1:ppn={cluster.threads}'

