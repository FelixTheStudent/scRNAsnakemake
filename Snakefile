# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

"""
Author: Felix Wertek
Affiliation: DKFZ, G200
Aim: Apply my STAR Snakemake workflow ('pipeline') to the data of this project:
	/icgc/dkfzlsdf/project/G200_hsc_inflamm_stress/sequencing/rna_sequencing
Date: 23rd June 2016
"""


##----------------------------------------------------##
# How to use this script
#      - fill in the 'local folder structure' below, 
#      - fill in the 'define samples...' setion below, then
#      execute this Snakefile with the following call:
#      snakemake -p -j 80 --cluster-config /path/to/config.json --cluster 'qsub -l walltime={cluster.time},mem={cluster.mem},nodes=1:ppn={cluster.threads}'


##----------------------------------------------------##
##	local folder structure where the samples lie
#
# Note: the slashes matter. If below something starts/ends with a '/', whatever you put in
#       there should do so, too.
FASTQDIR = "/icgc/dkfzlsdf/project/G200_hsc_inflamm_stress/sequencing/rna_sequencing/core/run160531_SN7001427_0222_AC9C6UANXX/"  #    M I N D   T H E    L A S T    S L A S H
RESULTS = "results/"  #    M I N D   T H E    L A S T    S L A S H
FASTQSUFFIX="_R1.fastq.gz"


#-------------------------------------------------------------##
## Variables declaration                          
## Declaring some variables used by STAR and other tools... 
## (GTF file, INDEX, ...)

INDEX = "/abi/data/wertek/kyewski/data/brennecke/reference_with_SpikeIn/indexes/star_erccGTF"
GTF = "/abi/data/wertek/kyewski/data/brennecke/reference_with_SpikeIn/databases/ensembl/Mus_musculus.GRCm38.84_ERCC92.gtf"
STARVERSION="/ibios/tbi_cluster/13.1/x86_64/bin/STAR-2.5"   # star2.5 can output a sorted BAM directly (although coordinate-sorted, not name-sorted.
OUTBAM="BAM Unsorted"





###########################################################################
##
##		include rule file and the configfile it needs
##
#########################################################################
##  the following declaration is necessary, because my 'rule file' 
#	"snakemakerules.STARrules.singleEnd" 
##  accesses the value 'threads' from the config.json file. For some reason I have not understood yet,
##  snakemake requires me to load the configfile here in the Snakefile instead of in the rule file. The 
##  access happens via this call:		{config[star_map][threads]}
configfile: "/home/wertek/scripts/snakemake/config.json"

include: "/home/wertek/scripts/snakemake/snakemakerules.STARrules.singleEnd"
##       RULE FILE     and     CONFIGURATION FILE
#	the above rule file defines rules, and for each rules the cluster resources are
#	defined in the config.json file that is being loaded in the rule file.
#	for this to work you also need to specify them in the --cluster parameter like this:
# snakemake -p -j 24 --cluster 'qsub -l walltime={cluster.time},mem={cluster.mem},nodes=1:ppn={cluster.threads}'
# 	Let me stress once more that for some reason, you still need to specify this on the command line: --cluster-config /path/to/config.js, even though you have the 'configfile: /path/to/*.json' above.


###########################################################################
##
##		define the samples on which to run the snakemake workflow
##
#########################################################################

SAMPLES, = glob_wildcards(FASTQDIR + "{sample,(AS-)\d*(-LR-1680)\d}"+ FASTQSUFFIX)
CELLS,CHIPNO = glob_wildcards(FASTQDIR + "AS-{cell,\d*}-LR-16{batch,\d*}"+ FASTQSUFFIX)

CELLNAMES=["cell"+s for s in CELLS]
BATCHNAMES="160531_SN7001427_0222_AC9C6UANXX"

# comment on the name scheme of this data.
# Name scheme:					AS-11ssss-LR-16bbb  
#						where:	`ssss` is unique sample number
#							 `bbb` is batch number
#		notes:	all samples (individual cells) in one batch (C1 chip) are either even or odd numbers.
#			Sample numbers (11`ssss`) go from (11)7271 until (11)8357 (discontinuously) and 
#			batches go from (LR-16)802 continuously up to 808. From one batch to another, the
#			sample numbers jump from even numbers to odd numbers.






rule all:
    # use this to get countmatrix with all cells: input: RESULTS + "countMatrices/160531_SN7001427_0222_AC9C6UANXXcountmatrix.R"
    # use the following line to index all BAM files. You can then look at them in IGV:
    input: expand(RESULTS + "mapping/{batch}/{sample}.coordsorted.bam.bai", sample=SAMPLES,batch=BATCHNAMES)

##########################
#         this rule can not be outsourced to the rule file, as otherwise SAMPLES and BATCHNAMES
#	  are not recognized. 
rule countmatrix:
    input: expand(RESULTS + "mapping/{batch}/{sample}.rawCounts.txt", sample=SAMPLES,batch=BATCHNAMES)
    output: protected(RESULTS + "countMatrices/{batch}countmatrix.R")
    # script: "/home/wertek/scripts/snakemake/countMatrix.R" # does not work until you manage to install rpy2 package.
    shell: "Rscript /home/wertek/scripts/snakemake/countMatrix_standalone.R {output} {input}"






