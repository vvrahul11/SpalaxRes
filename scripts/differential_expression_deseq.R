# This scripts take RNAseq count data of Spalax and Mouse and
# performs differential gene expression, and enrichment analysis#

# Install these libraries
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Mm.eg.db")
## BiocManager::install("edgeR") ##this for another version of pipeline
# BiocManager::install("DESeq2") ##this for current version of pipeline
# BiocManager::install("ReactomePA")
# BiocManager::install("topGO")
# BiocManager::install(c("gage", "FGNet", "pathview"))
# BiocManager::install("RCurl")
# BiocManager::install(c("RGtk2",
#            "RDAVIDWebService", "KEGGprofile"))
# BiocManager::install(("enrichplot"))
# install.packages("ggplot2")


#### Remove all objects ####
rm(list = ls())

#### Load libraries ####
#library(edgeR) ##this for another version of pipeline
library(DESeq2)
library(KEGGREST)
library(KEGGgraph)
library(org.Mm.eg.db)
library(ReactomePA)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
require(pathview)
library(ggplot2)


data(korg, package="pathview")
showWarnings = FALSE
load("data/kegg.gs.mmu.RData")

#----------------------------

#### DESeq2 differential expression analysis ####
calculate_diff <- function(d1, 
                           metadata,
                           groupNames,
                           mainDir, 
                           IDtype){
  #make DESeq2 object
  d <- DESeqDataSetFromMatrix(countData = d1,
                              colData = metadata[colnames(d1),],
                              design = ~ group)
  #normalize counts and find differential expressed genes
  dds <- DESeq(d)
  
  #find p-values and correct for multiple testing
  res <- results(dds,
                 pAdjustMethod = "fdr",
                 alpha = 0.01,
                 lfcThreshold = 0,
                 contrast=c("group", groupNames[2], groupNames[1]))
  
  summary(res)
  
  res_df = data.frame(res)
  res_df$padj[is.na(res_df$padj)] = 1 # replace NA in p.adj with 1
  res_df_filt = res_df[res_df$padj < 0.01, ] #threshold on p.adj
  res_df_filt = res_df_filt[abs(res_df_filt$log2FoldChange) > 2.5, ] #threshold on LFC
  
  #DESeq2::plotMA(res,  alpha = 0.01)
  #abline(h=c(-2.5, 2.5), col="red")
  
  genes = row.names(res_df_filt)
  genesUp = row.names(res_df_filt[res_df_filt$log2FoldChange > 0,])
  genesDown = row.names(res_df_filt[res_df_filt$log2FoldChange < 0,])
  
  backgroundIDs = row.names(res_df)
  foldchange.dge.genes = res_df_filt$log2FoldChange
  foldchange.dge.genesUp = res_df_filt[res_df_filt$log2FoldChange > 0,]$log2FoldChange
  foldchange.dge.genesDown = res_df_filt[res_df_filt$log2FoldChange < 0,]$log2FoldChange
  
  subDir <- paste(groupNames, 
                  collapse = "-", 
                  sep ="")
  
  if (file.exists(subDir)){
    setwd(file.path(mainDir, subDir))
    SaveCSV(res_df_filt, 
            res_df[!is.na(res_df$log2FoldChange),], #do not save genes with NA in LFC
            groupNames)
    
    Cluster.ProfileR(genes,
                     backgroundIDs,
                     foldchange.dge.genes,
                     IDtype,
                     "")
    Cluster.ProfileR(genesUp,
                     backgroundIDs,
                     foldchange.dge.genesUp,
                     IDtype,
                     "Up")
    Cluster.ProfileR(genesDown,
                     backgroundIDs,
                     foldchange.dge.genesDown,
                     IDtype,
                     "Down")
    
    write.csv(rownames(res_df), 
              "Allgenes.csv", 
              row.names = F)
    write.table(rownames(res_df), 
                "Allgenes.txt", 
                row.names = F, 
                quote= F)
    
  } else {
    dir.create(file.path(mainDir, subDir))
    setwd(file.path(mainDir, subDir))
    SaveCSV(res_df_filt, 
            res_df[!is.na(res_df$log2FoldChange),],
            groupNames)
    
    Cluster.ProfileR(genes,
                     backgroundIDs,
                     foldchange.dge.genes,
                     IDtype,
                     "")
    Cluster.ProfileR(genesUp,
                     backgroundIDs,
                     foldchange.dge.genesUp,
                     IDtype,
                     "Up")
    Cluster.ProfileR(genesDown,
                     backgroundIDs,
                     foldchange.dge.genesDown,
                     IDtype,
                     "Down")
    
    write.csv(rownames(res_df), 
              "Allgenes.csv", 
              row.names = F)
    write.table(rownames(res_df), 
                "Allgenes.txt", 
                row.names = F, 
                quote= F)
  }
  
  setwd(file.path(mainDir))
  
  mydiffgenelist = list(genesUp,
                        genesDown,
                        genes)
  return(mydiffgenelist)
}



#### Draw KEGG Maps ####
pathviewMake <- function(path.ids, 
                         exp.fc, 
                         out.suffix){
  load("C:/projects/Final Rat Trans Project Result2023/FinalProjectWorkAugust/kegg.gs.mmu.RData")
  data(korg, package="pathview")
  print("Going to plot pathview maps")
  print(path.ids)
  pv.out.list <- sapply(path.ids, 
                        function(pid) pathview(gene.data =  exp.fc, 
                                               pathway.id = pid,
                                               species = "mmu",
                                               out.suffix="edger",
                                               #limit      = list(gene=max(abs(exp.fc)), cpd=1)
                        )
  )
}

#### Draw venn diagram ####
makeVenn <-function(v1, v2, v3, v4){
  library(VennDiagram)
  png("venn_diagram.png", height=1500, width=2000)
  venn.plot <- venn.diagram(list(v1, v2, v3, v4), NULL,
                            fill=c("red", "green", "blue", "brown"),
                            height = 3000, 
                            width = 3000, 
                            resolution = 500,
                            alpha=c(0.5,0.5, 0.5, 0.5), 
                            cex = 4,
                            cat.cex=4, 
                            cat.fontface=2,
                            category.names=c("M.Wound", "M.Scar", "S.Wound", "S.Scar"))
  grid.draw(venn.plot)
  dev.off()
}



#### Save major foldchange files ####
SaveCSV <- function(filt_df, all_df, groupNames){
  write.csv(file = paste(paste(groupNames, collapse = "-", "-UP", sep ="") ### Upregulated
                         ,".csv", sep ="")
            , filt_df[filt_df$log2FoldChange > 0, c("log2FoldChange", "padj")])
  write.csv(file = paste(paste(groupNames, collapse = "-", "-DOWN", sep ="") ### Downregulated
                         ,".csv", sep ="")
            , filt_df[filt_df$log2FoldChange < 0, c("log2FoldChange", "padj")])
  write.csv(file = paste(paste(groupNames, collapse = "-", sep ="") ### All gene list
                         ,".csv", sep ="")
            , filt_df[, c("log2FoldChange", "padj")])
  myfold = data.frame("names.exp.fc" = rownames(all_df), 
                      "exp.fc" = all_df$log2FoldChange)
  rownames(myfold) = myfold$names.exp.fc
  write.csv(myfold, paste("foldchange", groupNames
                          , collapse = "-", sep ="", ".csv"),  )
}

### plots clusterProfiler ####
plotEnrichment <- function(ego, enrichmentType, geneList){
  
  #### barplot ####
  png(paste("barplot.", enrichmentType, ".png", sep = ""), height= 2500, width = 2000)
  mybar = barplot(ego, showCategory=20)
  plot(mybar)
  dev.off()
  
  #### cnet plot ####
  png(paste("Cnetplot.", enrichmentType, ".png", sep = ""), height= 2500, width = 2000)
  #ego_cnet <- pairwise_termsim(ego)
  # if you want to run the new version of clusterprofile use ego_cnet as parameter inside cnetplot
  cnetplot_ego = cnetplot(ego, categorySize="pvalue", 
                          foldChange=geneList,  
                          font.size = 8)#,
                          #layout = "star")
  plot(cnetplot_ego)
  dev.off()
  
  ### dot plot ####
  png(paste("dotplot.", enrichmentType, ".png", sep = ""), height= 2500, width = 2000)
  dotplot_ego = dotplot(ego, showCategory=10)
  plot(dotplot_ego)
  dev.off()
  
  
  #emapplot(ego, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
  #### EnrichMap ####
  png(paste("EnrichMap.", enrichmentType, ".png", sep = ""), height= 2500, width = 2000)
  #emapplot_ego = emapplot(ego_cnet) # for new version of enrichmap
  #plot(emapplot_ego)
  
  Enrich = enrichMap(ego) # for old version of enrichmap
  plot(Enrich)
  dev.off()
  
  #### Gorilla ####
  #png(paste("GorillaPlot.", enrichmentType, ".png", sep = ""), height= 2500, width = 2000)
  #plotGOgraph(ego)
  #dev.off()

  #### heatmap plot ####
  #png(paste("heatmapplot.", enrichmentType, ".png", sep = ""), height= 2500, width = 2000)
  #heatplot(ego, showCategory=10)
  #dev.off()
}


#### ClusterProfiler enrichment Analysis ####
Cluster.ProfileR <- function(genes, 
                             backgroundIDs, 
                             foldchange,
                             IDtype,
                             extension){
  if(IDtype == 2){
    # backgroundIDs1 = bitr(backgroundIDs, 
    #                      fromType="SYMBOL", 
    #                      toType="ENTREZID", 
    #                      OrgDb="org.Mm.eg.db")
    # eg = bitr(genes, 
    #           fromType="SYMBOL", 
    #           toType="ENTREZID", 
    #           OrgDb="org.Mm.eg.db")
    names(foldchange) = genes
    #### enrichGO ####
    ego <- enrichGO(gene          = genes,
                    universe      = backgroundIDs,
                    OrgDb         = org.Mm.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "fdr",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
    write.csv(ego, paste("GeneOntologyBP", extension, ".csv", sep=""))
    tryCatch({
      plotEnrichment(ego, "EnrichGO", foldchange)  
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    
    #### Reactome ####
    reactomePathways <- enrichPathway(gene= genes,
                                      organism = "mouse",
                                      pvalueCutoff = 0.05,
                                      pAdjustMethod = "fdr", 
                                      qvalueCutoff = 0.05,
                                      universe=backgroundIDs,
                                      minGSSize = 10,
                                      maxGSSize = 500)
    write.csv(reactomePathways, paste("reactome_pathways.csv", extension, ".csv", sep=""))
    
    tryCatch({
      plotEnrichment(reactomePathways, "Reactome", foldchange) 
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    
    #### KEGG ####
    keggpathways <- enrichKEGG(gene         = genes,
                               organism     = 'mouse',
                               pvalueCutoff = 0.05,  # This value is set hight to generate figures
                               pAdjustMethod     = "fdr",
                               universe = backgroundIDs,
                               minGSSize         = 10,
                               maxGSSize         = 500,
                               qvalueCutoff      = 0.05)
    write.csv(keggpathways, paste("Kegg_pathway.csv", extension, ".csv", sep=""))
    tryCatch({
      plotEnrichment(keggpathways, "KEGG", foldchange) 
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    #kegg.reactome.pathways= cbind(keggpathways, reactomePathways)
    #write.csv(kegg-reactome.pathways, "kegg-reactomePathways.csv")    
    
    kegg.pids = keggpathways[which(keggpathways$qvalue < 0.05), ]$ID
    if(length(kegg.pids > 1)){
      tryCatch({
        pathviewMake(kegg.pids, foldchange, "deseq")  
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
    
  }
}

#### merges d.f by rownames ####
merge.all <- function(x, ..., by = "row.names") {
  L <- list(...)
  for (i in seq_along(L)) {
    x <- merge(x, L[[i]], by = by)
    rownames(x) <- x$Row.names
    x$Row.names <- NULL
  }
  return(x)
}



cbind.fill<-function(...){
  nm <- list(...) 
  nm<-lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}


plot_pca <-function(GenewiseCounts, samples_metadata){
  d <- DESeqDataSetFromMatrix(countData = GenewiseCounts,
                               colData = samples_metadata,
                               design = ~ group)
  dds <- DESeq(d)
  vsd <- vst(dds, blind=TRUE)
  pca_fig = plotPCA(vsd, intgroup=c("group"), ntop = 15000) + theme_classic() + 
                ggtitle("PCA of RNA-seq Data") + 
                theme(plot.title = element_text(size=18, face = "bold", hjust = 0.5))
  print(pca_fig) #imediatly print the pca
  
  # Hierarchical Clustering Heatmap
  library(pheatmap)
  sampleDist <- dist(t(assay(vsd)))
  pheatmap(as.matrix(sampleDist), clustering_distance_rows=sampleDist,
           clustering_distance_cols=sampleDist, main="Sample Clustering Heatmap",
            labels_row = samples_metadata$group)
}



main <-function(){
  setwd("data")  
  expression = read.csv("DMBA-Mouse-Spalax-expressionGENESYMBOL.txt", row.names =1 )[, -c(31, 32)] #remove entrez ID and symbol
  expression <- expression[, -c(17, 23, 25)] #Remove RNASeq scar ovation sample TRUESEQ
  expressionID = read.csv("DMBA-Mouse-Spalax-expressionENTREZID.txt", row.names =1 )[, -c(31, 32)] #remove entrez ID and symbol
  expressionID <- expressionID[, -c(17, 23, 25)] #Remove RNASeq scar ovation sample TRUESEQ
  setwd("../DMBA")
  dim(expression)
  
  
  # Sample names
  # Mouse 1185, and Spalax 6042 looks like outliers
  S.names = c("m1", "m2","m3", "m4",  "m5", 
              #"m1185", "M20", "mM50", "mF60", # with outliers
              "M20", "mM50", "mF60", 
              "M1184","M1190", "M1188", 
              "Sp2866", "Sp6049","Sp6058", "Sp6053",  
              #"Sp6042", "Sp6026","Sp6029",  # with outliers
              "Sp6026","Sp6029", 
              "Sp6038_O", "Sp6028_O","Sp6032_O", 
              "Sp6038LD", "Sp6039I2", 
              "Sp6029H",  "Sp6038H", "Sp6039H")
  
  Counts  = expression[, S.names]
  CountsID  = expressionID[, S.names]
  group = c( rep(c("MM.control"), 5), 
             rep(c("M.wound"), 3),  #4 if including outliers 3 otherwise
             rep(c("M.scar"), 3), 
             rep(c("S.control"), 4), 
             rep(c("S.wound"), 2),  #3 if including outliers 2 otherwise
             rep(c("S.scar"), 3), 
             rep(c("S.LD"), 2), 
             rep(c("S.healing"), 3) )
  
  ##### Make metadata object for DESeq2 #########
  
  samples_metadata = data.frame("group" = group,
                                "sample_name" = S.names)
  rownames(samples_metadata) = samples_metadata$sample_name
  
  ## to be sure about the same order of rows in metadata and columns in counts matrix
  samples_metadata = samples_metadata[colnames(Counts), ]
  
  
  ###### Group Names #########
  ## choose pairs of groups to compare
  mygroup1 <- c("MM.control", "M.wound")
  mygroup2 <- c("MM.control", "M.scar")
  mygroup3 <- c("S.control", "S.wound")
  mygroup4 <- c("S.control", "S.scar")
  mygroup5 <- c("S.control", "S.healing")
  
  mygroups <- list(mygroup1, 
                   mygroup2,
                   mygroup3, 
                   mygroup4, 
                   mygroup5
  )
  
  d1Sym <- plot_pca(Counts, samples_metadata)
  
  ###### Normalize the count matrix, fit Negative Binomial GLM, estimate diff genes ############
  #d1Sym <- normalizeDataAndFindDiffGenes(Counts, group)
  #d1ID  <- normalizeDataAndFindDiffGenes(CountsID, group)
  
  for(i in c(1, 2)){
    IDtype = i
    if(i ==1){
      counts_table = Counts
      ## create directory Results_DESEQ2
      dir.create(paste(getwd(), "/Results/genesymbol/", sep = ""), recursive = TRUE)
      setwd(paste(getwd(), "/Results/genesymbol/", sep = ""))
      mainDir <- getwd()  
    }else{
      counts_table = CountsID
      ## create subdirectory in Results_DESEQ2
      dir.create("../entrezID/")
      setwd("../entrezID/")
      mainDir <- getwd()   
    }
    
    #### Create empty list and dataframe ####
    TotalDGEList = list()
    DGElist4save = data.frame(Empty = character()
                              ,stringsAsFactors=FALSE)
    for(groupNames in mygroups){
      #samples = combinations[,i]
      subset_metadata = samples_metadata[samples_metadata$group %in% groupNames,]
      subset_counts = counts_table[,rownames(subset_metadata)]
      mydiffgenelist = calculate_diff(subset_counts, 
                                      subset_metadata, 
                                      groupNames,
                                      mainDir, 
                                      i)
      TotalDGEList[[length(TotalDGEList)+1]] <- mydiffgenelist[[3]]
      DGElist4save <- cbind.fill(DGElist4save,
                                 mydiffgenelist[[1]],
                                 mydiffgenelist[[2]],
                                 mydiffgenelist[[3]],
                                 fill = NA)
      print(sampleNames)
      print(head(DGElist4save))
      #print(head(MYfoldchange))
      
    }  
    DGElist4save = DGElist4save[,-1]
    names(DGElist4save) <- rep(unique(unlist(mygroups))[-c(1,4)], each = 3)
    write.csv(DGElist4save, "DGEinEachCondition.csv")
    
    makeVenn(TotalDGEList[[1]], 
             TotalDGEList[[2]], 
             TotalDGEList[[3]], 
             TotalDGEList[[4]])
    
  }
}

main()

