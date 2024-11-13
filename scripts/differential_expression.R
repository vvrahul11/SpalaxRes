# This scripts take RNAseq count data of Spalax and Mouse and
# performs differential gene expression, and enrichment analysis#

# Install these libraries
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Mm.eg.db")
# BiocManager::install("edgeR")
# BiocManager::install("ReactomePA")
# BiocManager::install("topGO")
# BiocManager::install(c("gage", "FGNet", "pathview"))
# BiocManager::install(c("RGtk2", "RCurl",
#            "RDAVIDWebService","topGO", "KEGGprofile"))
# BiocManager::install(("enrichplot"))


#### Remove all objects ####
rm(list = ls())

#### Load libraries ####
library(edgeR)
library(KEGGREST)
library(KEGGgraph)
library(org.Mm.eg.db)
library(ReactomePA)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
require(pathview)

data(korg, package="pathview")
showWarnings = FALSE

#### edgeR TMM normalization ####
normalizeData <- function(GenewiseCounts, group){
  d <- DGEList(counts=GenewiseCounts, group=factor(group))
  countsPerMillion <- cpm(d)
  summary(countsPerMillion)
  # From biological point of view, a gene must be expressed at some minimal level before 
  # it is likely to be translated into a protein or to be considered biologically important. 
  # From a statistical point of view, genes with consistently low counts are very unlikely be 
  #assessed as significantly DE because low counts do not provide enough statistical 
  # evidence for a reliable judgement to be made.
  # As a rule of thumb, we require that a gene have a count of at least 10–15 in at 
  # least some libraries before it is considered to be expressed in the study.
  ## countsPerMillion>10 -> '1' is roughly equal to 10/L where L is the minimum library size in millions
  
  keep <- rowSums(countsPerMillion>1.7) >= 3 # rows with counts more than 100 and present in at least 2 samples
  # 10/6 = 1.7
  #The requirement of ≥3 libraries is because each group contains >= 3 replicates. 
  # This ensures that a gene will be retained if it is expressed in both the libraries belonging to any of the 5 groups.
  d <- d[keep,]
  dim(d)
  d$samples$lib.size <- colSums(d$counts)
  d$samples
  
  write.csv(rownames(d$counts), 
            "Allgenes-spalax-mouse.csv", 
            row.names = F)
  write.table(rownames(d$counts), 
              "Allgenes-spalax-mouse.txt", 
              row.names = F, 
              quote= F)
  d <- calcNormFactors( d, method="TMM")
  d$samples
  #apply(d$counts, 2, sum) #Total gene count per sample
  #d$samples$lib.size * d$samples$norm.factors
  
  #png("sample.comparison.png", width = 1000, height= 500)
  #plotMDS(d, method="bcv", col=as.numeric(d$samples$group), , labels = colnames( d$counts ))
  #legend("topleft", as.character(unique(d$samples$group)), pch=22)
  #dev.off()
  
  d1 <- estimateCommonDisp(d, verbose=T)
  ## Disp = 0.07776 , BCV = 0.2789
  names( d1 )
  d1$common.dispersion
  
  d1 <- estimateTagwiseDisp(d1)
  names(d1)
  summary(d1$tagwise.dispersion )
  #meanVarPlot <- plotMeanVar( d1 ,show.raw.vars=TRUE , show.tagwise.vars=TRUE , show.binned.common.disp.vars=FALSE , show.ave.raw.vars=FALSE , dispersion.method = "qcml" , NBline = TRUE , nbins = 100 , #these are arguments about what is plotted
  #                            pch = 16 , xlab ="Mean Expression (Log10 Scale)" , ylab = "Variance (Log10 Scale)" , main = "Mean-Variance Plot" ) #these arguments are to make it look prettier
  
  #png("plotBCV-dispersion1.png", width = 1000, height= 500)
  #plotBCV(d1)
  #dev.off()
  return(d1)
}


#### edgeR differential expression analysis ####
calculate_diff <- function(d1, 
                           sampleNames, 
                           kegg.gs.mmu, 
                           mainDir, 
                           IDtype){
  # Exact tests are used to identify differences in the means between two groups of gene counts. 
  # the default normalization method of edgeR is the Trimmed Mean of M-values (TMM), and TMM normalization factors do not take into account library sizes.                         
  et <- exactTest(d1, pair=c(sampleNames)) # compare groups 2 and 3
  samplegroup = d1$samples$group
  #print(topTags(et, n=10))
  edger.fc=et$table$logFC
  names(edger.fc)=rownames(et$table)
  exp.fc=edger.fc
  out.suffix="edger"
  FGNetCountData = d1$counts
  
  de1 <- decideTestsDGE(et, 
                        adjust.method="fdr", 
                        p.value=0.05, 
                        lfc=1.5)
  print(summary(de1))
  detags <- rownames(d1)[as.logical(de1)]

  plotSmear(et, de.tags=detags)
  abline(h=c(-1, 1), col="blue")
  
  #png(paste(paste(et$comparison, collapse = "-"), ".png", sep= ""), width = 500, height= 500)
  #plotMD(et, status=is.de, values=c(1,-1), col=c("red","blue"),
  #       legend="topright")
  #dev.off()
  
  genes = row.names(et$table[de1==1 | de1==-1,])
  genesUp = row.names(et$table[de1==1,])
  genesDown = row.names(et$table[de1==-1,])
  backgroundIDs = row.names(et$table)
  foldchange.dge.genes = et$table[de1==1 | de1==-1,]$logFC
  foldchange.dge.genesUp = et$table[de1==1,]$logFC
  foldchange.dge.genesDown = et$table[de1==-1,]$logFC
  
  subDir <- paste(et$comparison, 
                  collapse = "-", 
                  sep ="")
  
  if (file.exists(subDir)){
    setwd(file.path(mainDir, subDir))
    SaveCSV(et, 
            exp.fc, 
            de1)
    Cluster.ProfileR(genes,
                     backgroundIDs,
                     foldchange.dge.genes,
                     IDtype)
    
  } else {
    dir.create(file.path(mainDir, subDir))
    setwd(file.path(mainDir, subDir))
    SaveCSV(et, 
            exp.fc, 
            de1)
    
    Cluster.ProfileR(genes,
                     backgroundIDs,
                     foldchange.dge.genes,
                     IDtype)
  }
  setwd(file.path(mainDir))
  
  mydiffgenelist = list(row.names(et$table[de1==1,]),
                        row.names(et$table[de1==-1,]),
                        row.names(et$table[de1==1 | de1==-1,]))
  #geneNames = et$table[de1==1 | de1==-1,]
  #return(row.names(geneNames))
  return(mydiffgenelist)
}

#### Draw KEGG Maps ####
pathviewMake <- function(path.ids, 
                         exp.fc, 
                         out.suffix){
  load("C:/projects/Final Rat Trans Project Result2023/FinalProjectWorkAugust/kegg.gs.mmu.RData")
  data(korg, package="pathview")
  print("Going to plot pathview maps")
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
SaveCSV <- function(et, exp.fc, de1){
  write.csv(file = paste(paste(et$comparison, collapse = "-", "-UP", sep ="") ### Upregulated
                         ,".csv", sep ="")
            , et$table[de1==1,])
  write.csv(file = paste(paste(et$comparison, collapse = "-", "-DOWN", sep ="") ### Downregulated
                         ,".csv", sep ="")
            , et$table[de1==-1,])
  write.csv(file = paste(paste(et$comparison, collapse = "-", sep ="") ### All gene list
                         ,".csv", sep ="")
            , et$table[de1==1 | de1==-1,])
  myfold = data.frame(names(exp.fc), exp.fc)
  write.csv(myfold, paste("foldchange", et$comparison
                          , collapse = "-", sep ="", ".csv"),  )
  #png(paste(paste(et$comparison, collapse = "-"), ".png", sep= ""), width = 500, height= 500)
  #de1tags12 <- rownames(d1)[as.logical(de1)] 
  #plotSmear(et, de.tags=de1tags12)
  #abline(h = c(-2, 2), col = "blue")
  #dev.off()
  
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
                          font.size = 8,
                          layout = "star")
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
                             IDtype){
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
    write.csv(ego, "GeneOntologyBP.csv")
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
    write.csv(reactomePathways, "reactome_pathways.csv")
    tryCatch({
      plotEnrichment(reactomePathways, "Reactome", foldchange) 
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    
    #### KEGG ####
    keggpathways <- enrichKEGG(gene         = genes,
                               organism     = 'mouse',
                               pvalueCutoff = 0.5,  # This value is set hight to generate figures
                               pAdjustMethod     = "fdr",
                               universe = backgroundIDs,
                               minGSSize         = 10,
                               maxGSSize         = 500,
                               qvalueCutoff      = 0.5)
    write.csv(keggpathways, "Kegg_pathway.csv")
    tryCatch({
      plotEnrichment(keggpathways, "KEGG", foldchange) 
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    #kegg.reactome.pathways= cbind(keggpathways, reactomePathways)
    #write.csv(kegg-reactome.pathways, "kegg-reactomePathways.csv")    
    
    kegg.pids = keggpathways[which(keggpathways$qvalue < 0.5), ]$ID
    if(length(kegg.pids > 1)){
      tryCatch({
        pathviewMake(kegg.pids, foldchange, "edger")  
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
#### Data processing step 1 ####  
convert2entrezID <- function(expression){
  ### for pathview get the gene name id in entrez id 
  ### and get the fold change of the genes
  genesymbols <- row.names(expression)
  Gsymbol2ID = bitr(genesymbols, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  Gsymbol2ID = Gsymbol2ID[-8765,]
  row.names(Gsymbol2ID) = Gsymbol2ID$SYMBOL
  de_entrez = merge.all(expression, Gsymbol2ID)
  de_genesymbol = de_entrez
  row.names(de_entrez) = de_entrez$ENTREZID
  row.names(de_genesymbol) = de_entrez$SYMBOL
  write.csv(de_entrez, "DMBA-Mouse-Spalax-expressionENTREZID.txt")
  write.csv(de_genesymbol, "DMBA-Mouse-Spalax-expressionGENESYMBOL.txt")
}

compare_pathwaysbetweenTimePoints <- function(genes_list){
  #require(clusterProfiler)
  data(gcSample)
  res <- compareCluster(gcSample, fun="enrichPathway")
  png("CompareCluster_timePoints.png")
  plot(res)
  dev.off()
}

cbind.fill<-function(...){
  nm <- list(...) 
  nm<-lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

plot_pca <-function(GenewiseCounts, group){
  d <- DGEList(counts=GenewiseCounts, group=factor(group))
  d <- calcNormFactors(d, method="TMM")
  keep <- rowSums(cpm(d) > 1.7) >= 3
  d <- d[keep,]

  # Perform PCA
  logCPM <- cpm(d, log=TRUE)
  pca <- prcomp(t(logCPM), scale. = TRUE)

  # Plot PCA
  plot(pca$x[,1], pca$x[,2], col=as.numeric(factor(group)), pch=19,
      xlab="PC1", ylab="PC2", main="PCA of RNA-seq Data")
  legend("topright", legend=levels(factor(group)), col=1:length(unique(group)), pch=19)

  # Hierarchical Clustering Heatmap
  library(pheatmap)
  sampleDist <- dist(t(logCPM))
  pheatmap(as.matrix(sampleDist), clustering_distance_rows=sampleDist,
          clustering_distance_cols=sampleDist, main="Sample Clustering Heatmap")
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

  ###### Group Names #########
  # How the groups should be ordered. See here: https://morphoscape.wordpress.com/2022/08/09/downstream-bioinformatics-analysis-of-omics-data-with-edger/
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

  d1Sym <- plot_pca(Counts, group)

  ###### Normalize the count matrix  ############
  d1Sym <- normalizeData(Counts, group)
  d1ID <- normalizeData(CountsID, group)


  for(i in c(1, 2)){
    if(i ==1){
      d1 = d1Sym
      ## create directory Results
      setwd(paste(getwd(), "/Results/genesymbol/", sep = ""))
      mainDir <- getwd()  
    }else{
      d1 = d1ID
      ## create directory Results
      setwd("../entrezID/")
      mainDir <- getwd()   
    }
    
    #### Create empty list and dataframe ####
    TotalDGEList = list()
    DGElist4save = data.frame(Empty = character()
                              ,stringsAsFactors=FALSE)
    for(sampleNames in mygroups){
      #samples = combinations[,i]
      mydiffgenelist = calculate_diff(d1, 
                                      sampleNames, 
                                      kegg.gs.mmu, 
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
