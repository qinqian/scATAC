library(ArchR)
library('Cairo')
library(GetoptLong)
library(circlize)
library(DOSE)

addArchRThreads(threads = 10)
addArchRGenome("hg19")

setwd("~/langenau/stjude/scATAC-seq/")

inputFiles = Sys.glob("*/*fragments.tsv.gz")

ArrowFiles <- createArrowFiles(inputFiles = inputFiles,
                               sampleNames = dirname(inputFiles),
                               filterTSS = 4, #Dont set this too high because you can always increase later
                               filterFrags = 1000,
                               addTileMat = TRUE,
                               ## bcTag='RG',
                               addGeneScoreMat = TRUE)

doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1)


proj <- ArchRProject(
    ArrowFiles = ArrowFiles,
    outputDirectory = "StJude_RMS",
    copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

paste0("Memory Size = ", round(object.size(proj) / 10^6, 3), " MB")


getAvailableMatrices(proj)

head(proj$cellNames)

head(proj$Sample)

quantile(proj$TSSEnrichment)

proj[1:100, ]

proj[proj$cellNames[1:100], ]

idxSample <- BiocGenerics::which(proj$Sample %in% "SJRHB010927_X1")
cellsSample <- proj$cellNames[idxSample]
proj[cellsSample, ]

pdf("RMS_ATAC_QC2.pdf")
p1 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "ridges")
p2 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
)
p3 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "ridges"
   )
p4 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
p5 <- plotFragmentSizes(ArchRProj = proj)
p6 <- plotTSSEnrichment(ArchRProj = proj)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
dev.off()


saveArchRProject(ArchRProj = proj, outputDirectory = "Save-Proj", load = FALSE)


proj2 <- filterDoublets(proj)

proj2 <- addIterativeLSI(
    ArchRProj = proj2,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)

proj2 <- addHarmony(
    ArchRProj = proj2,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)

proj2 <- addClusters(
    input = proj2,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8
)

table(proj2$Clusters)

cM <- confusionMatrix(paste0(proj2$Clusters), paste0(proj2$Sample))

library(pheatmap)
cM <- cM / Matrix::rowSums(cM)

pdf("RMS_ATAC_clusters.pdf")
p <- pheatmap::pheatmap(
                   mat = as.matrix(cM), 
                   color = paletteContinuous("whiteBlue"), 
                   border_color = "black")
print(p)
dev.off()


proj2 <- addUMAP(
    ArchRProj = proj2,
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")


proj2 <- addUMAP(
    ArchRProj = proj2,
    reducedDims = "Harmony", 
    name = "UMAPHarmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

p3 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
ggAlignPlots(p3, p4, type = "h")
plotPDF(p1, p2, p3, p4, name = "RMS_ATAC_Plot-UMAP-Clusters.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)

proj2 <- addClusters(
    input = proj2,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "ClustersHarmony",
    resolution = 0.8
)

cM <- confusionMatrix(paste0(proj2$ClustersHarmony), paste0(proj2$Sample))
p3 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = "ClustersHarmony", embedding = "UMAPHarmony")
plotPDF(p3, p4, name = "RMS_ATAC_Plot-UMAP-HarmonyClusters.pdf", ArchRProj = proj2, addDOC = FALSE, width = 5, height = 5)

pdf("RMS_ATAC_Harmonyclusters.pdf")
p <- pheatmap::pheatmap(
                   mat = as.matrix(cM), 
                   color = paletteContinuous("whiteBlue"),
                   display_numbers=T,
                   border_color = "black")
print(p)
dev.off()


proj3 = proj2[proj2$Clusters!='C1' | proj2$ClustersHarmony!='C1', ]

## proj3 <- addHarmony(
##     ArchRProj = proj3,
##     reducedDims = "IterativeLSI",
##     name = "Harmony",
##     groupBy = "Sample",
##     force=T
## )

proj3 <- addUMAP(
    ArchRProj = proj3,
    reducedDims = "Harmony", 
    name = "UMAPHarmony_subset", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    ## force=T
)

proj3 <- addClusters(
    input = proj3,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "ClustersHarmony_subset",
    resolution = 1.0 #, force=T
)

cM <- confusionMatrix(paste0(proj3$ClustersHarmony_subset), paste0(proj3$Sample))

pdf("RMS_ATAC_Proj3_Harmonyclusters.pdf")
p <- pheatmap::pheatmap(
                   mat = as.matrix(cM), 
                   color = paletteContinuous("whiteBlue"),
                   display_numbers=T,
                   border_color = "black")
print(p)
dev.off()

p3 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony_subset")
p4 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "ClustersHarmony_subset", embedding = "UMAPHarmony_subset")
plotPDF(p3, p4, name = "RMS_ATAC_Plot-UMAP-HarmonyClusters2.pdf", ArchRProj = proj3, addDOC = FALSE, width = 5, height = 5)

markersGS <- getMarkerFeatures(
    ArchRProj = proj3, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "ClustersHarmony_subset",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

gene.modules <- Sys.glob('~/langenau/projects/01_sc_rms/final_annotations/gene_modules/*txt')
gene.list <- lapply(gene.modules, scan, what='')

rms.markers = unlist(gene.list)

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  ## labelMarkers=c("ADAMTS1", "TNNT2", "CD44", "BUB3"),
  labelMarkers = intersect(rms.markers, getFeatures(proj2)),
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapGS, name = "RMS-GeneScores-Marker-Heatmap", width = 24, height = 6, ArchRProj = proj3, addDOC = FALSE)

p <- plotEmbedding(
    ArchRProj = proj3, 
    colorBy = "GeneScoreMatrix", 
    ## name = intersect(rms.markers, getFeatures(proj3)), #c("MYOG", "MYL4", "TNNT3", "TNNT2"),
    name = c("CD44", "EGFR", "MEOX2", "OGN", "POSTN", # EMT
             "MYOG", "MYL4", "TNNT3", "TNNT2", "MEF2C", "DES", # Muscle
             "BUB3", "BUB1", "E2F1", "CDK1", "MKI67", "TOP2A",
             "NOTCH1", "NOTCH3"),
    embedding = "UMAPHarmony_subset",
    quantCut = c(0.01, 0.99),
    imputeWeights = NULL
)

plotPDF(plotList = p, 
    name = "Plot-UMAP-Marker-Genes-WO-Imputation.pdf", 
    ArchRProj = proj3,
    addDOC = FALSE, width = 5, height = 5)

p <- plotBrowserTrack(
    ArchRProj = proj3, 
    groupBy = "ClustersHarmony_subset", 
    geneSymbol = c("CD44", "EGFR", "MEOX2", "OGN", "POSTN", # EMT
                   "MYOG", "MYL4", "TNNT3", "TNNT2", "MEF2C", "DES", # Muscle
                   "BUB3", "BUB1", "E2F1", "CDK1", "MKI67", "TOP2A", # Prolif
                   "NOTCH1", "NOTCH3"),
    upstream = 100000,
    downstream = 100000
)

plotPDF(plotList = p, 
        name = "RMS-Plot-Tracks-Marker-Genes.pdf",
        ArchRProj = proj3,
        addDOC = FALSE, width = 12, height = 6)

allgenes = readRDS("allgenes_signatures.rds")

library(clusterProfiler)
## integrative clustering cannot annotate cells well
clusters = list()
for (clust in names(markerList)) {
    print(dim(markerList[[clust]]))
    genes = markerList[[clust]]$name
    y1 = enricher(genes, TERM2GENE=allgenes[,-2], maxGSSize=1000, minGSSize = 5)
    if (is.null(y1))
        clusters[[as.character(clust)]] <- y1
    else{
        y1 = y1@result
        rownames(y1) = paste0(clust, ".", rownames(y1))
        y1$clusters = clust
        clusters[[as.character(clust)]] <- y1
    }
}

clusters = do.call('rbind', clusters)
clusters = clusters %>% mutate(FoldEnrichment = parse_ratio(GeneRatio)/parse_ratio(BgRatio))
clusters = subset(clusters, p.adjust <= 0.05 & qvalue <= 0.05 & Count >= 3 & FoldEnrichment >= 1.1)
## clusters = subset(clusters, p.adjust <= 0.01 & qvalue <= 0.01 & Count >= 3 & FoldEnrichment >= 1.2)
readr::write_csv(clusters, path='RMS_StJude_integrative_cluster_enrichr_filter_annotation.csv')

## C1: Muscle
## C2: Ground
## C3: Ground
## C4: Muscle
## C5: Ground
## C6: Ground
## C7: Ground
## C8: Muscle
## C9: Ground
## C10: Ground
## C11: EMT
## C12: Ground
## C13: Ground
## C14: Ground
## C15: Ground
## C16: Ground
## C17: Ground
## C18: Ground
## C19: Ground
## C20: Ground
## C21: Ground
proj3$Cellstate = proj3$ClustersHarmony_subset

cell_names = factor(proj3$Cellstate)
integrative_annotation = c("Muscle", "Ground", "Ground", "Muscle", "Ground", "Ground",
                           "Ground", "Muscle", "Ground", "Ground", "EMT", "Ground",
                           "Ground", "Ground", "Ground", "Ground", "Ground", "Ground",
                           "Ground", "Ground", "Ground")
names(integrative_annotation) = paste0("C", seq(1, 21))
levels(cell_names) = integrative_annotation[levels(cell_names)]

proj3$Cellstate = cell_names

p4 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "Cellstate", embedding = "UMAPHarmony_subset")
plotPDF(p4, name = "RMS_ATAC_Plot-UMAP-IntegrativeClustersAnnotation.pdf", ArchRProj = proj3, addDOC = FALSE, width = 5, height = 5)

p <- plotBrowserTrack(
    ArchRProj = proj3, 
    groupBy = "Cellstate", 
    geneSymbol = c("CD44", "EGFR", "MEOX2", "OGN", "POSTN", # EMT
                   "MYOG", "MYL4", "TNNT3", "TNNT2", "MEF2C", "DES", # Muscle
                   "BUB3", "BUB1", "E2F1", "CDK1", "MKI67", "TOP2A", # Prolif
                   "NOTCH1", "NOTCH3"),
    upstream = 100000,
    downstream = 100000
)

plotPDF(plotList = p, 
        name = "RMS-Plot-Tracks-Marker-Genes-IntegrativeClustersAnnotation.pdf",
        ArchRProj = proj3,
        addDOC = FALSE, width = 12, height = 6)

## remove_doublet_or_not = c()
## names(remove_doublet_or_not) = unique(proj3$Sample)

## split samples for annotation
## to see if Prolif cell state exists
sample_list = list()
for (samp in unique(proj3$Sample)) {
    proj_tmp = proj3[proj3$Sample==samp, ]
    proj_tmp <- addIterativeLSI(
        ArchRProj = proj_tmp,
        useMatrix = "TileMatrix", 
        name = "IterativeLSI", 
        iterations = 4, 
        clusterParams = list( #See Seurat::FindClusters
            resolution = c(0.2), 
            sampleCells = 10000, 
            n.start = 10
        ), 
        varFeatures = 25000, 
        dimsToUse = 1:30,
        force=T,
    )
    proj_tmp <- addUMAP(
        ArchRProj = proj_tmp,
        reducedDims = "IterativeLSI", 
        name = "UMAP", 
        nNeighbors = 30, 
        minDist = 0.5, 
        metric = "cosine",
        force=T
    )
    proj_tmp <- addClusters(
        input = proj_tmp,
        reducedDims = "IterativeLSI",
        method = "Seurat",
        name = "ClustersHarmony_individualsample",
        resolution = 0.8,
        force=T
    )
    sample_list[[samp]] = proj_tmp
}


for (samp in unique(proj3$Sample)) {
    proj_tmp = sample_list[[samp]]
    markersGS <- getMarkerFeatures(
        ArchRProj = proj_tmp, 
        useMatrix = "GeneScoreMatrix", 
        groupBy = "ClustersHarmony_individualsample",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        testMethod = "wilcoxon"
    )
    markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
    ## integrative clustering cannot annotate cells well
    clusters = list()
    for (clust in names(markerList)) {
        print(dim(markerList[[clust]]))
        genes = markerList[[clust]]$name
        y1 = enricher(genes, TERM2GENE=allgenes[,-2], maxGSSize=1000, minGSSize = 5)
        if (is.null(y1))
            clusters[[as.character(clust)]] <- y1
        else{
            y1 = y1@result
            rownames(y1) = paste0(clust, ".", rownames(y1))
            y1$clusters = clust
            clusters[[as.character(clust)]] <- y1
        }
    }
    clusters = do.call('rbind', clusters)
    clusters = clusters %>% mutate(FoldEnrichment = parse_ratio(GeneRatio)/parse_ratio(BgRatio))
    clusters = subset(clusters, p.adjust <= 0.05 & qvalue <= 0.05 & Count >= 3 & FoldEnrichment >= 1.1)
    readr::write_csv(clusters, path=paste0(samp, 'RMS_StJude_cluster_enrichr_filter_annotation.csv'))
    heatmapGS <- markerHeatmap(
        seMarker = markersGS, 
        cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
        ## labelMarkers=c("ADAMTS1", "TNNT2", "CD44", "BUB3"),
        labelMarkers = intersect(rms.markers, getFeatures(proj_tmp)),
        transpose = TRUE
    )
    ## ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
    plotPDF(heatmapGS, name = paste0(samp, "RMS-GeneScores-Marker-Heatmap"), width = 24, height = 6, ArchRProj = proj_tmp, addDOC = FALSE)
    p <- plotEmbedding(
        ArchRProj = proj_tmp, 
        colorBy = "GeneScoreMatrix", 
        name = c("CD44", "EGFR", "MEOX2", "OGN", "POSTN", # EMT
                 "MYOG", "MYL4", "TNNT3", "TNNT2", "MEF2C", "DES", # Muscle
                 "BUB3", "BUB1", "E2F1", "CDK1", "MKI67", "TOP2A",
                 "NOTCH1", "NOTCH3"),
        embedding = "UMAP",
        quantCut = c(0.01, 0.99),
        imputeWeights = NULL
    )
    plotPDF(plotList = p, 
            name = paste0(samp, "Plot-UMAP-Marker-Genes-WO-Imputation.pdf"),
            ArchRProj = proj_tmp,
            addDOC = FALSE, width = 5, height = 5)
    p <- plotBrowserTrack(
        ArchRProj = proj_tmp, 
        groupBy = "ClustersHarmony_individualsample",
        geneSymbol = c("CD44", "EGFR", "MEOX2", "OGN", "POSTN", # EMT
                       "MYOG", "MYL4", "TNNT3", "TNNT2", "MEF2C", "DES", # Muscle
                       "BUB3", "BUB1", "E2F1", "CDK1", "MKI67", "TOP2A", # Prolif
                       "NOTCH1", "NOTCH3"),
        upstream = 100000,
        downstream = 100000
    )
    plotPDF(plotList = p, 
            name = paste0(samp, "RMS-Plot-Tracks-Marker-Genes-IntegrativeClustersAnnotation.pdf"),
            ArchRProj = proj_tmp,
            addDOC = FALSE, width = 12, height = 6)
}

## Use macs2
## pathToMacs2 <- findMacs2()
## proj3 <- addGroupCoverages(ArchRProj = proj3, groupBy = "ClustersHarmony_subset")
## proj3 <- addReproduciblePeakSet(
##     ArchRProj = proj3, 
##     groupBy = "ClustersHarmony_subset", 
##     pathToMacs2 = pathToMacs2
## )
## proj3 <- addMotifAnnotations(ArchRProj = proj3,
##                              motifSet = "cisbp", name = "Motif")
