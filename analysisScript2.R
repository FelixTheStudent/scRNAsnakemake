# Scripts steps:
#  incorporate hierarchical clustering into milsomSet
#  plot PCA and t-SNE with it
#  find ideal number of clusters
#  ... there's so much more you can do ...!

#  after filtering genes in previous script, you can save it out.
#  Then load it here:


  milsomSet <- readRDS("....rds")





#' ## For further analysis: normalized, transformed matrices
#' This is useful to compute clusters, PCA, etc., independently from the scater / scran packages 

   milsomexprs <- milsomSet[rownames(hvg.out.all)] # normalized with scran's deconvolution (Lun, Genome Biol 2016)
   milsomexprs <- t(exprs(milsomexprs))
   # create treatment-vectors for above matrices to plot with color=treatment   (does not apply to your data, charles, just leave it out...)
   treatment <- revalue(substr(rownames(milsomexprs),18,18), 
                        c("3"="pbs", "4"="pbs", "5"="pbs", 
                          "6"="treated", "7"="treated", "8"="treated") )



#' ## Compute clusters of cells:
 HC.euclidean <- dist(milsomexprs, method="euclidean")
 HC.euclideantree <- hclust(HC.euclidean, method="ward.D2")
 library(dynamicTreeCut)
 HC.euclidean <- unname(cutreeDynamic(HC.euclideantree, distM=as.matrix(HC.euclidean), verbose=0))
 HC.euclidean <- as.factor(HC.euclidean)
# note that below, I have a section on how to find 'ideal number of clusters' with NBClust - this is cool but perhaps for later timepoints.
 

#############################################################################3
#' ## Add new information to the milsomSet 
#' Here I add the NbClust cluster definitions (if you computed them, see below) to the milsomSet. Also, I put the hvg library size, etc.,
#' to later be able to display them in PCA / DiffusionMaps or to compute correlation with the first components.
 
 # It's good to check library size as a potential confounder interfering with PCA. To be rigorous, I on top
 # also want to use the library size of the 1764 HVGs 
    hvg1764_counts <- rowSums(milsomexprs)
    # number of detected HVGs in each cell:
    hvg1764_features <- apply(milsomexprs, 1, function(cell) sum(cell >0))

 # add information to milsomSet:   (NOTE: you can use the 'mutate' function instead of the following complicated code, see https://rdrr.io/bioc/scater/man/mutate.html)
   milsomPheno <- pData(milsomSet)
    milsomPheno <- cbind( 
		# add the clustering from above:
			HC.euclidean,
               # add library sizes, in different versions:
                         hvg1764_counts, hvg1764_features, 
		# do the following commented lines only if you followed section 'ideal number of clusters' below...
			   # # add NbClust's ideal cluster definitions, according to Silhouette index:    
                           # clusters.complete.eucl.s, clusters.ward.eucl.s, clusters.kmeans.eucl.s,
                           # # add NbClust's ideal cluster definitions, according to Dunn index:
                           # clusters.complete.eucl.d, clusters.ward.eucl.d, clusters.kmeans.eucl.d, clusters.kmeans.spear.d,clusters.kmeans.spear.s,
               # append the previous milsomPheno information as well, of course:			
                         milsomPheno)
    pd <- new("AnnotatedDataFrame", data = milsomPheno)
    pData(milsomSet) <- pd

    # confirm it worked:
    head(varLabels(milsomSet), n=10)  # should show you "clusters...." 
  


  
 # now you can color anything inside varlabels(milsomSet) in PCAs and t-SNE: 
plotPCA(milsomSet, feature_set=rownames(hvg.out.all), color=HC.euclidean)
# plotTSNE( ...)   # see ? plotTSNE, it's awesome!


# Now felix is curious, once you're here send him a screenshot ;)

 
#############################################################################3
#' #     ideal number of clusters with NbClust 
  # this is perhaps for later analysis steps, not for the first quick look on the data. That's why I commented it out for now.

   #   library(NbClust) 
   # 
   # nbwrapper <- function( dataset=milsomexprs, index, min.nc=2, max.nc=8){
   #     # nbwrapper computes ideal number of clusters with euclidean matrix (4 different methods) and spearman distance
   #     # (one method, as the other two do not work for me in this package implementation). The 4 methods are the following:
   #            # average: "UPGMA"
   #            # ward.D2: "Ward's method"
   #            # kmeans
   #            # complete: "complete-linkage clustering"
   #    # computation based on euclidean distance, 3 methods
   #     c<- NbClust(dataset, distance = "euclidean", min.nc=min.nc, max.nc=max.nc, method = "complete", index = index)
   #     av<- NbClust(dataset, distance = "euclidean", min.nc=min.nc, max.nc=max.nc, method = "average", index = index)
   #     wa<- NbClust(dataset, distance = "euclidean", min.nc=min.nc, max.nc=max.nc, method = "ward.D2", index = index)
   #     k<- NbClust(dataset, distance =  "euclidean", min.nc=min.nc, max.nc=max.nc, method = "kmeans", index = index)
   #    # computation based on spearman distance, 1 method: 
   #    HC.spearman <-cor(t(dataset), method='spearman')
   #    HC.spearman <- (1-HC.spearman)/2 
   #    ##   the following 2 options (s.av, s.wa) do not work with datasets this large. I have tried to do it only on the 217 treated cells,
   #     #   but still get this error: 
   #     #         Error in if (is.na(n) || n > 65536L) stop("size cannot be NA nor exceed 65536") : 
   #     #         missing value where TRUE/FALSE needed
   #    # s.av<- NbClust(dataset, diss=HC.spearman, distance = NULL, min.nc=min.nc, max.nc=max.nc, method = "average", index = index)
   #    # s.wa<- NbClust(dataset, diss=HC.spearman, distance = NULL, min.nc=min.nc, max.nc=max.nc, method = "ward.D2", index = index)
   #    s.k<- NbClust(dataset, diss=HC.spearman, distance = NULL, min.nc=min.nc, max.nc=max.nc, method = "kmeans", index = index)
   #   # extract information you want to return: index holds score for each cluster size, clusters holds ideal cluster assignments for each cell 
   #    index.all <- rbind(c$All.index ,av$All.index,wa$All.index,k$All.index,s.k$All.index)# s.av$All.index,s.wa$All.index,
   #    rownames(index.all) <- c("eucli.complete" ,"eucli.average","eucli.ward","eucli.kmeans",  "spear.kmeans")#"spear.av","spear.wa",
   #    clusters.all <- rbind(c$Best.partition, av$Best.partition,wa$Best.partition,k$Best.partition,s.k$Best.partition)#s.av$Best.partition,s.wa$Best.partition,
   #    rownames(clusters.all) <- c( "eucli.complete", "eucli.average","eucli.ward","eucli.kmeans",  "spear.kmeans")# "spear.av","spear.wa",
   #    return(list(index.all, clusters.all))
   # } 
   # 
   # 
   # # Compute ideal clusters with Dunn and Silhouette indices, or load the pre-computed ones:
   #   # compute and save:
   #      # dunn.list <- nbwrapper(milsomexprs, "dunn", max.nc=12)
   #      # silhouette.list <- nbwrapper(milsomexprs, "silhouette", max.nc=12)
   #      # saveRDS(dunn.list, file=paste(DIR, "nbclust_dunn-list",format(Sys.time(), "_%Y%m%d-%H:%M."), "rds", sep='')  )
   #      # saveRDS(silhouette.list, file=paste(DIR, "nbclust_silhouette-list",format(Sys.time(), "_%Y%m%d-%H:%M."), "rds", sep='')  )
   #   # load the above:
   #   dunn.list <- readRDS('/home/felix/PhDother/milsomData/stressedCells/octobersprint/resultsanddata/20161217_Draft1_firstFigs/nbclust_dunn-list_20161217-16:18.rds')
   #   silhouette.list <- readRDS('/home/felix/PhDother/milsomData/stressedCells/octobersprint/resultsanddata/20161217_Draft1_firstFigs/nbclust_silhouette-list_20161217-16:18.rds')
   #   
   # # plot index scores for each number of clusters:  
   #   # massage score data for dunn index:
   #   dunn.index <- dunn.list[[1]]
   #   dunn.index <- scale( t(dunn.index), scale=F)
   #   dunn.index <- dunn.index[,-2] # exclude average linkage, as all cells except 1 in same cluster (not shown here)
   #   dunn.index <- melt(dunn.index)
   #   dunn.index <- cbind(dunn.index, rep("dunn",dim(dunn.index)[1]))
   #   colnames(dunn.index) <- c("Number_of_Clusters","distance_and_method",  "IndexScore", "Index")
   #   # massage score data for silhouette index:
   #   silh.index<- silhouette.list[[1]]
   #   silh.index <- scale( t(silh.index), scale=F)
   #   silh.index <- silh.index[,-2]# exclude average linkage, as all cells except 1 in same cluster (not shown here)
   #   silh.index<- melt(silh.index)
   #   silh.index <- cbind(silh.index, rep("Silhouette", dim(silh.index)[1]))
   #   colnames(silh.index) <- c("Number_of_Clusters","distance_and_method",  "IndexScore", "Index")
   #   # bring both into same dataframe:
   #   indices <- rbind(dunn.index, silh.index)
   #   # plot ideal scores of all: 
   #   indices <- ggplot(indices, aes(x=Number_of_Clusters, y=IndexScore, color=distance_and_method, linetype=Index)) + 
   #     geom_line() + scale_x_continuous(breaks=as.numeric(2:12),) + ylab("Index Score (centered)")
   #   #plotit(indices, DIR, "NBclust_1764HVG", knit= knitall)
   #  
   #    dunn.clusters <- as.data.frame(t(dunn.list[[2]]))
   #    silhouette.clusters <- as.data.frame(t(silhouette.list[[2]]))
   #    
   #    
   #    #' I do not show plots for 'average linkage', as the ideal clusters were not meaningful (all cells except 1 in same cluster),
   #    #' for both Dunn and Silhouette indices. Therefore, they are also left out in all following analysis. The following clusters,
   #    #' however, are interesting and will be investigated further. To this end, we extract the ideal clusterings and, further below,
   #    #' add them to the milsomSet.
   # 
   # # Silhouette:
   #         #   1 eucli.complete
   #         #   3 eucli.ward
   #         #   4 eucli.kmeans
   #         #   5 spear kmeans
   # clusters.complete.eucl.s <- as.factor(silhouette.clusters[,1])
   # clusters.ward.eucl.s <- as.factor(silhouette.clusters[,3])
   # clusters.kmeans.eucl.s <- as.factor(silhouette.clusters[,4])
   # clusters.kmeans.spear.s <- as.factor(silhouette.clusters[,5])
   # # Dunn:
   #         #   1 eucli.complete
   #         #   3. eucli.ward
   #         #   4. eucli.kmeans
   #         #   5. spear.kmeans
   # clusters.complete.eucl.d <- as.factor(dunn.clusters[,1])
   # clusters.ward.eucl.d <- as.factor(dunn.clusters[,3])
   # clusters.kmeans.eucl.d <- as.factor(dunn.clusters[,4])
   # clusters.kmeans.spear.d <- as.factor(dunn.clusters[,5])
 
 
 



   
  
  
   
 
 
 
#######################################################
#' ## Compare different cluster definitions

 
 #' 
 pc.1764x424 <- prcomp(milsomexprs,center=T,scale=F)
 # extract scores (cell's values in PC1/2) and loads (best genes have high scores)
 scores.1764x424 <- as.data.frame(pc.1764x424$x)
 loads.1764x424 <- as.data.frame(pc.1764x424$rotation) # coefficients of linear combinations of variables
 
 #' Compute t-SNE           
  seed=100
  tsne.hvg <- plotTSNE(milsomSet, feature_set=rownames(hvg.out.all), colour_by="hvg1764_features")
  tsne.scores <- tsne.hvg$data 
  stopifnot(all( rownames(tsne.scores) == colnames(milsomSet) )) # makes sure rows of tsne values and in milsomSet match
 
# ... 
