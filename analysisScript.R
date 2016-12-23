## figuresprint_header.R will be incorporated into the 'figuresprint.R' script onces it's done.

# Steps:
# - read inHTSeq count matrix. Remove HTSeq-count lines and batch 802
# - extract bulk samples from it and save it as bulkmatrix
# - remove bulks, use scater to filter genes and cells
# - compare filtered cells to 28 badCells from my own QC back in the day
# - compute HVGs on cleaned matrix and save it as single-cell matrix

# the bulkmatrix and single-cell matrix created as described above will be uploaded as "processed files" to GEO,
# together with the 'raw data files' (fastq.gz).

library(Biobase)
# Load scater and required packages
packagelist<- c("data.table", "ggplot2", "knitr", "matrixStats", "MASS",
                "plyr", "reshape2", "rjson", "testthat", "viridis",
                "rhdf5", 
                "scater",
                "scran")
sapply(packagelist, require, character.only = TRUE)
packageVersion('scater'); packageVersion('scran') # bioconductor-3.4 versions for both are 1.2, released 19th October 2016. Get devel versions in ~ November!




#' Load count matrix as outputted by HTseq-count:
load("/icgc/dkfzlsdf/analysis/G200/wertek/milsom/G200_hsc_inflamm_stress/results/countMatrices/160531_SN7001427_0222_AC9C6UANXXcountmatrix.R")
countMatrix <- countMatrix


#' remove batch 802 and HTSeqCount-specific rows
countMatrix <- countMatrix[, !substr(colnames(countMatrix), 16,18)=="802"]
countMatrix <- countMatrix[-grep("__",rownames(countMatrix)),]


# separate bulk and single-cell data (I know the bulkSamples from the Metafile of our core facility)
bulkSamples <- c("AS-117588-LR-16803","AS-117747-LR-16804","AS-117900-LR-16805","AS-118051-LR-16806","AS-118210-LR-16807","AS-118357-LR-16808")
 # FYI: for batch 802, the bulk sample is "AS-117341-LR-16802" 
bulkMatrix <- countMatrix[,bulkSamples]
countMatrix <- countMatrix[,! colnames(countMatrix) %in% bulkSamples]



#' # Scater / Scran

batchNames <- substr(colnames(countMatrix),16,18)
pData <- as.data.frame(cbind( Batch=batchNames,
                              Treatment=revalue(as.character(batchNames),
                                                c("803"="PBS","804"="PBS","805"="PBS","806"="treated","807"="treated","808"="treated")),
                              CellType = rep("LT-HSC",length(colnames(countMatrix)))
))
rownames(pData) <- colnames(countMatrix)
pMeta <- data.frame(labelDescription=
                      c( " One Batch equals one 96-well plate in C1 machine, as well as one sequencing run. ",  "Batches 3, 4, 5 are PBS controls, batches 6,7,8 come from pI:pC injected mice.", "All cells were sorted as Long-term HSCs, see methods."),
                    row.names=c("Batch","Treatment","CellType"))
milsomphenoData <- new("AnnotatedDataFrame", data=pData, varMetadata= pMeta)
# Create ExpressionSet and SCESet for this:
milsomSet <- newSCESet(countData = countMatrix, phenoData = milsomphenoData) 







#' before the next step, you should filter genes and cells as desribed below (section 'pasted from another script, untestetd').
#' If you are very impatient, do this right away to compute differentially expressed genes:

#' # use this cut-off:
fdr=0.05
# ERCC-fit: 
 var.fit.all <- trendVar(milsomSet, trend="loess", use.spikes=T, span=0.2)
 
# compute HVGs:
 var.out.all <- decomposeVar(milsomSet, var.fit.all)
 hvg.out.all <- var.out.all[which(var.out.all$FDR <= fdr & var.out.all$bio >= 0.5),]
 hvg.out.all <- hvg.out.all[order(hvg.out.all$bio, decreasing=TRUE),]
# get significantly variable genes:
hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >= 0.5),]
hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),]
nrow(hvg.out)
hist(hvg.out$FDR)


 plotfitGenes <- function(var.out){
 # fitplotGenes plots your variance fit as in the tutorial. It adds a blue line
   # to indicate that genes were fitted (black dots), not the spike-ins (red dots)
   plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression",
      ylab="Variance of log-expression", main=deparse(substitute(var.out)))
 o <- order(var.out$mean)
 lines(var.out$mean[o], var.out$tech[o], col="dodgerblue", lwd=2)
 spike.fit <- trendVar(milsomSet, use.spikes=TRUE) # To compute spike-in variances.
 points(spike.fit$mean, spike.fit$var, col="red", pch=16)
 }
 # using spikeins. Differences to plotfitGenes-function above:
 # you don't have to do a trendVar on spikeins, and the line is red ;)
 plotfitSpikes <- function(var.out, var.fit){
   plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression",
        ylab="Variance of log-expression", main=deparse(substitute(var.out)))
   points(var.fit$mean, var.fit$var, col="red", pch=16)
   o <- order(var.out$mean)
   lines(var.out$mean[o], var.out$tech[o], col="red", lwd=2)
 }


# fit model to spike-ins with different values for 'span' (plots are further below)
 fitSpikes <- trendVar(milsomSet, trend="loess", use.spikes=T)
 fitSpikes0.2 <- trendVar(milsomSet, trend="loess", use.spikes=T, span=0.2)
 fitSpikes0.3 <- trendVar(milsomSet, trend="loess", use.spikes=T, span=0.3)
 fitSpikes0.4 <- trendVar(milsomSet, trend="loess", use.spikes=T, span=0.4)
 fitSpikes1.2 <- trendVar(milsomSet, trend="loess", use.spikes=T, span=1.2)
 varSpikes <- decomposeVar(milsomSet, fitSpikes)
 varSpikes0.2 <- decomposeVar(milsomSet, fitSpikes0.2)
 varSpikes0.3 <- decomposeVar(milsomSet, fitSpikes0.3)
 varSpikes0.4 <- decomposeVar(milsomSet, fitSpikes0.4)
 varSpikes1.2 <- decomposeVar(milsomSet, fitSpikes1.2)


# fit model to genomic genes (NOT spike-ins):
fitGenes <- trendVar(milsomSet, trend="loess", use.spikes=F)
fitGenes0.2 <- trendVar(milsomSet, trend="loess", use.spikes=F, span=0.2)
fitGenes0.3 <- trendVar(milsomSet, trend="loess", use.spikes=F, span=0.3)
fitGenes0.4 <- trendVar(milsomSet, trend="loess", use.spikes=F, span=0.4)
fitGenes1.2 <- trendVar(milsomSet, trend="loess", use.spikes=F, span=1.2)
varGenes <- decomposeVar(milsomSet, fitGenes)
varGenes0.2 <- decomposeVar(milsomSet, fitGenes0.2)
varGenes0.3 <- decomposeVar(milsomSet, fitGenes0.3)
varGenes0.4 <- decomposeVar(milsomSet, fitGenes0.4)
varGenes1.2 <- decomposeVar(milsomSet, fitGenes1.2)

# plot the fits to spike-ins and endogenous genes, then pick best one:
par(mfrow=c(5,2))
plotfitGenes(varGenes)
plotfitSpikes(varSpikes, fitSpikes)
plotfitGenes(varGenes0.2)
plotfitSpikes(varSpikes0.2, fitSpikes0.2)
plotfitGenes(varGenes0.3)
plotfitSpikes(varSpikes0.3, fitSpikes0.3)
plotfitGenes(varGenes0.4)
plotfitSpikes(varSpikes0.4, fitSpikes0.4)
plotfitGenes(varGenes1.2)
plotfitSpikes(varSpikes1.2, fitSpikes1.2)
par(mfrow=c(1,1))
#' From this I conclude that fitting **to the Spike-ins, with span = 0.2** is the best and I should use it.
var.out <- varSpikes0.2

# get significantly variable genes:
hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >= 0.5),]
hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),]
nrow(hvg.out)
hist(hvg.out$FDR)










####################################################################################
##     pasted from another script, untestetd:
#

#' ## Filter Cells (ERCCs, mitoch, lib size)
# for ERCCs:
is.spike <- grepl("^ERCC", rownames(milsomSet))
# for mitochondrial genes:
library(biomaRt)
ensembl=useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org")  # uses version 85 at 23rd July, I used version 84 for mapping.
# use the mouse Mart:
ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
# retrieve all 37 genes annotated on the murine mitochondrial chromosome:
#  (in the past, this only worked for me connected to eduroam, not ethernet-cable - if so leave out or find better solution)
mtGenes<-getBM(attributes= "ensembl_gene_id",
               filters="chromosome_name",
               values="MT", mart=ensembl)
is.mito <- rownames(milsomSet) %in% mtGenes

milsomSet <- calculateQCMetrics(milsomSet, feature_controls=list(ERCC=is.spike, Mt=is.mito))
isSpike(milsomSet) <- "ERCC"

par(mfrow=c(1,2))
hist(milsomSet$total_counts/1e6, xlab="Library sizes (millions)", main="",
     breaks=20, col="grey80", ylab="Number of cells")
hist(milsomSet$total_features, xlab="Number of expressed genes", main="",
     breaks=20, col="grey80", ylab="Number of cells")
par(mfrow=c(1,1))

libsize.drop <- isOutlier(milsomSet$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(milsomSet$total_features, nmads=3, type="lower", log=TRUE)
mito.drop <- isOutlier(milsomSet$pct_counts_feature_controls_Mt, nmads=3, type="higher")
spike.drop <- isOutlier(milsomSet$pct_counts_feature_controls_ERCC, nmads=3, type="higher")
milsomSet <- milsomSet[,!(libsize.drop | feature.drop | mito.drop | spike.drop)]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
           ByMito=sum(mito.drop), BySpike=sum(spike.drop), Remaining=ncol(milsomSet))

#' ## filter Genes 
ave.counts <- rowMeans(counts(milsomSet))
AvExprthreshold <- 0.5
keep <- ave.counts >= AvExprthreshold # this is applied below, after more plots with all genes
table(keep)

hist(log10(ave.counts), breaks=100, main="", col="grey80",
     xlab=expression(Log[10]~"average count"))
abline(v=log10(AvExprthreshold), col="blue", lwd=2, lty=2)

plotQC(milsomSet, type = "highest-expression", n=50)

numcells <- nexprs(milsomSet, byrow=TRUE)
smoothScatter(log10(ave.counts), numcells, xlab=expression(Log[10]~"average count"),
              ylab="Number of expressing cells")
is.ercc <- isSpike(milsomSet, type="ERCC")
points(log10(ave.counts[is.ercc]), numcells[is.ercc], col="red", pch=16, cex=0.5)

milsomSet <- milsomSet[keep,]

















