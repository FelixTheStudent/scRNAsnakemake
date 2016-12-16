#!/usr/bin/env Rscript
## After my snakemake workflow for mapping with STAR and obtaining read counts with htseq-count
## has run through, there will be many output-files *.rawCounts.txt. This script is called by
## a further snakemake rule and passed the following arguments:
#	- first argument is for the output: the filepath and -name (path is relative to the Snakefile)
#	- all following arguments are filepaths and -names to the *.rawCounts.txt files
## This script will take all input files, read them in as dataframes, merge them into an expression
## matrix (rows = genes [ENSEMBL IDs], columns = cells [names from file type]) and save it as RData.


# as of 20th June 2016, the snakemake rule calling this script is placed within the Snakefile in the
# project directories. Not sure if I'll leave it like this or whether I put it into the rule file
# in /home/wertek/scripts/snakemake

####
#
#	note that as a side-effect of the merging below, the genes are sorted alphabetically by their gene names (special cases first, then ensembl IDs, then ERCCs)
#	As of 20th June 2016, I see no problem with this.
###


library(stringr)    # this is for the regEx below, to assign filename as column names 
# read in all arguments (1st is output file, all following are input files)
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error. Also remember whether there is
# more than one input file or not.
justOneInputFile=FALSE
if (length(args)<2) {
  stop("At least two file paths/names (one output file and one input file) must be supplied (input file).n", call.=FALSE)
##########################################################3
} else if (length(args)==2) {
	justOneInputFile=TRUE
}


# get output file names:
outputfile <- args[1]
# get input file names:
indexLast = length(args)
#inputfiles <- args[2:indexLast]


## 
## create countMatrix as dataframe by reading in the first cell:
##
countMatrix <- read.table(args[2],header=FALSE,sep="\t",row.names=1)
numberOfGenes <- nrow(countMatrix)  # use this below to verify merging went well (i.e. that you didn't loose any rows).
# get the filename and put it as column name:
filename <- basename(args[2])
filename <- str_match(filename,"(.*).rawCounts.txt")[,2]  # FYI: [,1] gives entire string in which it found the match, [,2] is then whatever is within the (brackets).
colnames(countMatrix)<-filename


##
## for loop reading in and merging by rownames the DFs for all cells
##
if (justOneInputFile) { 
	print("R-Script read out rawCounts for one input file.")
}else{
remainingFiles=args[3:indexLast]
print("Reading out rawCounts for this many files: ")
print(length(remainingFiles)+1)

for (file in remainingFiles){
	DFonecell <- read.table(file,header=FALSE,sep="\t",row.names=1)
	filename <- basename(file)
	filename <- str_match(filename,"(.*).rawCounts.txt")[,2]  # FYI: [,1] gives entire string in which it found the match, [,2] is then whatever is within the (brackets).
	colnames(DFonecell)<-filename

##################################
#### 		merge with countMatrix
#		Note that this step sorts the rownames alphabetically
#################################
countMatrix <- merge(countMatrix,DFonecell,by="row.names",all=FALSE)   # gives shortened rownumbers if cells have different gene names. So below, I can control everything went well by checking rowlength.
# the above command adds an additional column called 'Row.names', which we want to revert:
row.names(countMatrix) <- countMatrix$Row.names	
countMatrix$Row.names <- NULL



}
# do a very simple check if the countMatrix looks ok (initially, I did not trust merging too much)
if(nrow(countMatrix)==numberOfGenes) { print("Merging seems to have gone well. final CountMatrix has as many rows as first cell from first input file")}
}

# write out the countMatrix to a file that can be loaded and worked with.
save(countMatrix,file=outputfile)

